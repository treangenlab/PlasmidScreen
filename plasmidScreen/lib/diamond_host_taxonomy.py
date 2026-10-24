"""DIAMOND blastx ORF intervals and host taxonomy for codon CAI.

DIAMOND finds translated ORFs on nucleotide reads (--min-orf) and reports
query coordinates (qstart/qend) plus subject taxonomy (staxids).
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
import subprocess
from pathlib import Path
from typing import Iterable, Optional, Sequence

from plasmidScreen.lib.exceptions import DiamondError, InputError
from plasmidScreen.lib.models import MEMORY_MODE, MEMORY_CONFIG

# outfmt 6: qseqid sseqid qstart qend pident length evalue bitscore staxids sscinames
DIAMOND_OUTFMT = (
    "6",
    "qseqid",
    "sseqid",
    "qstart",
    "qend",
    "pident",
    "length",
    "evalue",
    "bitscore",
    "staxids",
    "sscinames",
)
DIAMOND_MIN_COLUMNS = 9  # through staxids (0-based index 8)


@dataclass(frozen=True)
class DiamondHit:
    read_id: str
    subject_id: str
    qstart: int
    qend: int
    pident: float
    aln_length: int
    evalue: float
    bitscore: float
    staxids: tuple[str, ...]
    sscinames: tuple[str, ...] = ()

    @property
    def query_start(self) -> int:
        """0-based start on read (half-open interval with query_end)."""
        return min(self.qstart, self.qend) - 1

    @property
    def query_end(self) -> int:
        """0-based exclusive end on read."""
        return max(self.qstart, self.qend)

    @property
    def strand(self) -> str:
        return "-" if self.qstart > self.qend else "+"


@dataclass(frozen=True)
class OrfInterval:
    read_id: str
    start: int
    end: int
    strand: str
    hits: tuple[DiamondHit, ...]

    @property
    def length_bp(self) -> int:
        return max(0, self.end - self.start)


def _parse_staxids_field(raw: str) -> tuple[str, ...]:
    s = raw.strip()
    if not s or s == "0":
        return tuple()
    parts = [p.strip() for p in s.replace(",", ";").split(";")]
    return tuple(p for p in parts if p and p != "0")


def _parse_sscinames_field(raw: str) -> tuple[str, ...]:
    s = raw.strip()
    if not s:
        return tuple()
    parts = [p.strip() for p in s.replace(",", ";").split(";")]
    return tuple(p for p in parts if p)


def read_diamond_output(path: str | Path) -> list[str]:
    """Load DIAMOND outfmt 6 TSV lines from a saved file."""
    path = Path(path)
    try:
        return path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise InputError(f"Could not read DIAMOND output {path}: {exc}") from exc


def diamond_technology_flags(sequencing_technology: int) -> list[str]:
    """
    DIAMOND blastx flags for long-read (0) vs short-read (1) inputs.

    Long-read / ONT: frameshift + range-culling (error-prone, multi-ORF queries).
    Short-read: more-sensitive search without frameshift alignment.
    """
    if sequencing_technology == 0:
        return [
            "--fast",
            "--range-culling",
            "--culling-overlap",
            "50",
            "--range-cover",
            "50",
            "--frameshift",
            "15",
            "--min-orf",
            "30",
            "--masking",
            "seg",
            "-k",
            "1",
        ]
    if sequencing_technology == 1:
        return [
            "--more-sensitive",
            "--min-orf",
            "10",
            "--masking",
            "0",
            "--top",
            "5",
        ]
    raise ValueError(
        f"Expecting sequence technology to be either 0 or 1. Instead received "
        f"{sequencing_technology}. Where 0 is ONT/long-read and 1 is short-read."
    )


def resolve_diamond_lines(
        reads_path: str | Path,
        diamond_db: str | Path | None,
        *,
        run_diamond: bool = True,
        mem_mode: int = 0,
        output_path: str | Path | None = None,
        debug_write_output: bool = False,
        threads: int = 4,
        sequencing_technology: int = 1,
) -> tuple[list[str], Path | None]:
    """
    Run DIAMOND or load a saved TSV.

    Returns (lines, path) where path is set when reading/writing ``output_path``.
    """
    if debug_write_output and output_path is None:
        raise ValueError(
            "output_path is required when debug_write_output=True."
        )
    if not run_diamond:
        if output_path is None:
            raise ValueError(
                "output_path is required when run_diamond=False (precomputed DIAMOND TSV)."
            )
        path = Path(output_path)
        if not path.is_file():
            raise InputError(f"DIAMOND output not found: {path}")
        logging.info("Loading DIAMOND output from %s", path)
        return read_diamond_output(path), path

    if diamond_db is None:
        raise ValueError("diamond_db is required when run_diamond=True.")

    write_path = Path(output_path) if debug_write_output else None
    lines = run_diamond_blastx(
        reads_path,
        diamond_db,
        threads=threads,
        mem_mode=mem_mode,
        output_path=write_path,
        sequencing_technology=sequencing_technology,
    )
    return lines, write_path


def run_diamond_blastx(
        reads_path: str | Path,
        diamond_db: str | Path,
        threads: int = 4,
        mem_mode: int = 0,
        output_path: str | Path | None = None,
        sequencing_technology: int = 1,
) -> list[str]:
    """Run DIAMOND blastx; ORF detection is built in (--min-orf)."""
    reads_path = Path(reads_path)
    diamond_db = Path(diamond_db)

    if mem_mode == 0:  # MEMORY_CONFIG.get(mem_mode) == MEMORY_MODE.LOW.value:
        block_size = "4"
        index_chunks = "3"
    elif mem_mode == 1:  # MEMORY_CONFIG.get(mem_mode) == MEMORY_MODE.MEDIUM.value:
        block_size = "6"
        index_chunks = "1"
    elif mem_mode == 2:  # MEMORY_CONFIG.get(mem_mode) == MEMORY_MODE.HIGH.value:
        block_size = "12"
        index_chunks = "1"
    else:
        raise ValueError(
            "memory_mode was expecting a value between 0-2 where:"
            "\nlow(0) <= 32GB  RAM "
            "\nmedium(1) <= 64GB  RAM "
            "\nhigh(2) <= 128GB RAM"
        )

    technology_specific_commands = diamond_technology_flags(sequencing_technology)
    if sequencing_technology == 0:
        logging.info(
            "Running DIAMOND in long-read/ONT mode (--seq-type 0)."
        )
    else:
        logging.info(
            "Running DIAMOND in short-read mode (--seq-type 1)."
        )

    if not reads_path.is_file():
        raise InputError(f"DIAMOND query file not found: {reads_path}")
    if not diamond_db.exists():
        raise InputError(f"DIAMOND database not found: {diamond_db}")
    cmd: list[str] = [
        "diamond",
        "blastx",
        "-q",
        str(reads_path),
        "-d",
        str(diamond_db),
        "--evalue",
        "10",
        "--threads",
        str(threads),
        "--block-size",
        block_size,
        "--index-chunks",
        index_chunks,
        "--salltitles",
        "-f",
        *DIAMOND_OUTFMT,
    ]
    cmd.extend(technology_specific_commands)

    try:
        if output_path is not None:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            cmd.extend(["-o", str(out)])
            logging.info("Running DIAMOND: %s", " ".join(cmd))
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info("DIAMOND output written to %s", out)
            return read_diamond_output(out)

        logging.info("Running DIAMOND: %s", " ".join(cmd))
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return proc.stdout.splitlines()
    except FileNotFoundError as exc:
        logging.error("DIAMOND executable not found on PATH")
        raise DiamondError(
            127, "diamond executable not found; is DIAMOND installed and on PATH?"
        ) from exc
    except subprocess.CalledProcessError as exc:
        logging.error(
            "DIAMOND failed (code %s). stderr: %s",
            exc.returncode,
            (exc.stderr or "").strip(),
        )
        raise DiamondError(exc.returncode, exc.stderr) from exc


def parse_diamond_tsv(lines: Iterable[str]) -> dict[str, list[DiamondHit]]:
    """Parse DIAMOND outfmt 6 TSV into read_id -> hits (sorted by bitscore)."""
    out: dict[str, list[DiamondHit]] = {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < DIAMOND_MIN_COLUMNS:
            continue
        read_id = parts[0]
        try:
            qstart = int(parts[2])
            qend = int(parts[3])
            pident = float(parts[4])
            aln_length = int(parts[5])
            evalue = float(parts[6])
            bitscore = float(parts[7])
        except ValueError:
            logging.debug("Skipping malformed DIAMOND line: %s", line[:120])
            continue
        staxids = _parse_staxids_field(parts[8])
        sscinames = _parse_sscinames_field(parts[9]) if len(parts) > 9 else tuple()
        hit = DiamondHit(
            read_id=read_id,
            subject_id=parts[1],
            qstart=qstart,
            qend=qend,
            pident=pident,
            aln_length=aln_length,
            evalue=evalue,
            bitscore=bitscore,
            staxids=staxids,
            sscinames=sscinames,
        )
        out.setdefault(read_id, []).append(hit)
    for rid in out:
        out[rid].sort(key=lambda h: (-h.bitscore, h.query_start, h.query_end))
    return out


def hits_to_orf_intervals(
        read_id: str,
        hits: Sequence[DiamondHit],
        *,
        merge_gap_bp: int = 30,
        min_orf_len_bp: int = 30,
) -> list[OrfInterval]:
    """Merge DIAMOND hit coordinates into ORF intervals on the read (coordinates only)."""
    if not hits:
        return []
    sorted_hits = sorted(hits, key=lambda h: (h.query_start, h.query_end))
    intervals: list[OrfInterval] = []
    cur_start = sorted_hits[0].query_start
    cur_end = sorted_hits[0].query_end
    cur_strand = sorted_hits[0].strand
    cur_hits: list[DiamondHit] = [sorted_hits[0]]
    for h in sorted_hits[1:]:
        if h.strand == cur_strand and h.query_start <= cur_end + merge_gap_bp:
            cur_end = max(cur_end, h.query_end)
            cur_hits.append(h)
        else:
            if cur_end - cur_start >= min_orf_len_bp:
                intervals.append(
                    OrfInterval(
                        read_id=read_id,
                        start=cur_start,
                        end=cur_end,
                        strand=cur_strand,
                        hits=tuple(cur_hits),
                    )
                )
            cur_start, cur_end, cur_strand, cur_hits = (
                h.query_start,
                h.query_end,
                h.strand,
                [h],
            )
    if cur_end - cur_start >= min_orf_len_bp:
        intervals.append(
            OrfInterval(
                read_id=read_id,
                start=cur_start,
                end=cur_end,
                strand=cur_strand,
                hits=tuple(cur_hits),
            )
        )
    return intervals


def majority_taxid(taxids: Sequence[str], ) -> Optional[str]:
    """
    Gets the majority taxid along a subsequence from diamond.
    """
    clean = [t for t in (str(x) for x in taxids) if t not in ("0", "")]
    if not clean:
        return None
    counts: dict[str, int] = {}
    for t in clean:
        counts[t] = counts.get(t, 0) + 1
    best = max(counts.items(), key=lambda kv: kv[1])
    return best[0]


def staxids_from_hits(hits: Sequence[DiamondHit]) -> list[str]:
    taxids: list[str] = []
    for h in hits:
        taxids.extend(h.staxids)
    return taxids


def infer_orfs_and_host_taxids(
        diamond_lines: Iterable[str],
        merge_gap_bp: int = 30,
        min_orf_len_bp: int = 30,
) -> tuple[dict[str, list[OrfInterval]], dict[str, Optional[str]]]:
    """
    Parse DIAMOND outfmt 6 into per-read ORF intervals and host taxids.

    ORFs are merged hit intervals on each read; host taxid is the majority
    ``staxids`` value across all hits for that read.
    """
    hits_by_read = parse_diamond_tsv(diamond_lines)
    orfs_by_read: dict[str, list[OrfInterval]] = {}
    host_by_read: dict[str, Optional[str]] = {}

    for read_id, hits in hits_by_read.items():
        orfs_by_read[read_id] = hits_to_orf_intervals(
            read_id,
            hits,
            merge_gap_bp=merge_gap_bp,
            min_orf_len_bp=min_orf_len_bp,
        )
        host_by_read[read_id] = majority_taxid(staxids_from_hits(hits))

    return orfs_by_read, host_by_read
