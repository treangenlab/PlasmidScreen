"""DIAMOND-based ORF inference and host taxonomy assignment (SeqScreen-Nano style).

This module is build/runtime offline-friendly: it runs local DIAMOND and uses local
taxonomy parents (from codon reference build) to resolve LCA.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
import subprocess
from pathlib import Path
from typing import Iterable, Optional, Sequence


@dataclass(frozen=True)
class DiamondHit:
    read_id: str
    qstart: int
    qend: int
    bitscore: float
    staxids: tuple[str, ...]


@dataclass(frozen=True)
class OrfInterval:
    read_id: str
    start: int
    end: int
    hits: tuple[DiamondHit, ...]
    taxid: Optional[str]

    @property
    def length_bp(self) -> int:
        return max(0, self.end - self.start)


def _parse_staxids_field(raw: str) -> tuple[str, ...]:
    # DIAMOND staxids can be like "9606" or "9606;10090" or empty.
    s = raw.strip()
    if not s or s == "0":
        return tuple()
    parts = [p.strip() for p in s.replace(",", ";").split(";")]
    return tuple(p for p in parts if p and p != "0")


def run_diamond_blastx(
    reads_path: str | Path,
    diamond_db: str | Path,
    *,
    threads: int = 4,
    extra_args: Sequence[str] | None = None,
) -> list[str]:
    """Run DIAMOND blastx and return output lines (TSV).

    We require an outfmt that includes taxonomy IDs via `staxids`.
    """
    reads_path = Path(reads_path)
    diamond_db = Path(diamond_db)
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
        "200",
        "--index-chunks",
        "1",
        "--salltitles",
        "--more-sensitive",
        "--min-orf",
        "10",
        "--masking",
        "0",
        "--top",
        "5",
        "-f",
        "100"
    ]
    if extra_args:
        cmd.extend(list(extra_args))

    logging.info("Running DIAMOND: %s", " ".join(cmd))
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return proc.stdout.splitlines()


def parse_diamond_tsv(lines: Iterable[str]) -> dict[str, list[DiamondHit]]:
    """Parse DIAMOND outfmt 6 lines into read_id -> hits list."""
    out: dict[str, list[DiamondHit]] = {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        read_id = parts[0]
        try:
            qstart = int(parts[1])
            qend = int(parts[2])
            bitscore = float(parts[3])
        except ValueError:
            continue
        staxids = _parse_staxids_field(parts[4])
        hit = DiamondHit(
            read_id=read_id,
            qstart=min(qstart, qend),
            qend=max(qstart, qend),
            bitscore=bitscore,
            staxids=staxids,
        )
        out.setdefault(read_id, []).append(hit)
    # Deterministic ordering: higher bitscore first, then coordinate.
    for rid in out:
        out[rid].sort(key=lambda h: (-h.bitscore, h.qstart, h.qend))
    return out


def hits_to_orf_intervals(
    read_id: str,
    hits: Sequence[DiamondHit],
    *,
    merge_gap_bp: int = 30,
    min_orf_len_bp: int = 90,
) -> list[tuple[int, int, list[DiamondHit]]]:
    """Merge hits into ORF intervals (approximation of SeqScreen-Nano candidate ORFs)."""
    if not hits:
        return []
    sorted_hits = sorted(hits, key=lambda h: (h.qstart, h.qend))
    intervals: list[tuple[int, int, list[DiamondHit]]] = []
    cur_start = sorted_hits[0].qstart
    cur_end = sorted_hits[0].qend
    cur_hits: list[DiamondHit] = [sorted_hits[0]]
    for h in sorted_hits[1:]:
        if h.qstart <= cur_end + merge_gap_bp:
            cur_end = max(cur_end, h.qend)
            cur_hits.append(h)
        else:
            if cur_end - cur_start >= min_orf_len_bp:
                intervals.append((cur_start, cur_end, cur_hits))
            cur_start, cur_end, cur_hits = h.qstart, h.qend, [h]
    if cur_end - cur_start >= min_orf_len_bp:
        intervals.append((cur_start, cur_end, cur_hits))
    return intervals


def _lineage(taxid: str, parents: dict[str, str]) -> list[str]:
    out = []
    cur = str(taxid)
    seen: set[str] = set()
    while cur and cur not in seen and cur not in ("0", ""):
        out.append(cur)
        seen.add(cur)
        parent = parents.get(cur)
        if parent is None or parent == cur:
            break
        cur = parent
    return out


def lca_taxid(taxids: Sequence[str], parents: dict[str, str]) -> Optional[str]:
    """Lowest common ancestor using taxonomy_parents mapping."""
    clean = [t for t in (str(x) for x in taxids) if t not in ("0", "")]
    if not clean:
        return None
    lineages = [_lineage(t, parents) for t in clean]
    if any(not lin for lin in lineages):
        return None
    common = set(lineages[0])
    for lin in lineages[1:]:
        common.intersection_update(lin)
        if not common:
            return None
    # choose deepest: max depth in first lineage order
    for t in lineages[0]:
        if t in common:
            return t
    return None


def majority_taxid(
    taxids: Sequence[str],
    parents: dict[str, str],
) -> Optional[str]:
    """Majority vote; if ties/no majority, fall back to LCA."""
    clean = [t for t in (str(x) for x in taxids) if t not in ("0", "")]
    if not clean:
        return None
    counts: dict[str, int] = {}
    for t in clean:
        counts[t] = counts.get(t, 0) + 1
    best = max(counts.items(), key=lambda kv: kv[1])
    if best[1] > len(clean) / 2:
        return best[0]
    return lca_taxid(clean, parents)


def orf_taxid_from_hits(hits: Sequence[DiamondHit], parents: dict[str, str]) -> Optional[str]:
    """Assign a taxid to an ORF interval from its hits (SeqScreen-Nano majority/LCA)."""
    taxids: list[str] = []
    for h in hits:
        taxids.extend(list(h.staxids))
    return majority_taxid(taxids, parents)


def read_host_taxid_from_orfs(orfs: Sequence[OrfInterval], parents: dict[str, str]) -> Optional[str]:
    """Assign a host taxid to a read from ORF taxids (majority/LCA)."""
    taxids = [o.taxid for o in orfs if o.taxid]
    return majority_taxid([t for t in taxids if t is not None], parents)


def infer_orfs_and_host_taxids(
    diamond_lines: Iterable[str],
    *,
    taxonomy_parents: dict[str, str],
    merge_gap_bp: int = 30,
    min_orf_len_bp: int = 90,
) -> tuple[dict[str, list[OrfInterval]], dict[str, Optional[str]]]:
    """End-to-end inference: DIAMOND lines -> ORF intervals + per-read host taxid."""
    hits_by_read = parse_diamond_tsv(diamond_lines)
    orfs_by_read: dict[str, list[OrfInterval]] = {}
    host_by_read: dict[str, Optional[str]] = {}

    for read_id, hits in hits_by_read.items():
        merged = hits_to_orf_intervals(
            read_id, hits, merge_gap_bp=merge_gap_bp, min_orf_len_bp=min_orf_len_bp
        )
        orfs: list[OrfInterval] = []
        for start, end, interval_hits in merged:
            taxid = orf_taxid_from_hits(interval_hits, taxonomy_parents)
            orfs.append(
                OrfInterval(
                    read_id=read_id,
                    start=start,
                    end=end,
                    hits=tuple(interval_hits),
                    taxid=taxid,
                )
            )
        orfs_by_read[read_id] = orfs
        host_by_read[read_id] = read_host_taxid_from_orfs(orfs, taxonomy_parents)

    return orfs_by_read, host_by_read

