"""
Reference database lookup for ScreenResult enrichment.

Decoupled from ML / k-mer / codon modeling. Engineered reads are searched with
``blastn`` against a nucleotide reference (FASTA subject or BLAST DB). Hits
above a bitscore threshold are retained; the sample's top-N reference sequences
(by # reads covered, bitscore tiebreak) are kept for reporting and visualization.

Codon CAI (DIAMOND blastx) remains a separate stage and is not used here.
"""
from __future__ import annotations

import logging
import os
import subprocess
import tempfile
import time
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

from Bio import SeqIO

from plasmidScreen.lib.models import (
    NO_DATABASE_MATCH_NOTE,
    QueryMatchConfidence,
    SampleReferenceRank,
    ScreenResult,
)

logger = logging.getLogger(__name__)

# BLAST outfmt 6 columns.
_BLAST_QSEQID = 0
_BLAST_SSEQID = 1
_BLAST_PIDENT = 2
_BLAST_LENGTH = 3
_BLAST_QSTART = 6
_BLAST_QEND = 7
_BLAST_SSTART = 8
_BLAST_SEND = 9
_BLAST_EVALUE = 10
_BLAST_BITSCORE = 11
_BLAST_MIN_COLS = 12

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn", ".fas", ".faa"}


@dataclass(frozen=True)
class HitFilterConfig:
    """Thresholds for BLAST hit retention and sample-level reference ranking."""

    min_bitscore: float = 50.0
    top_n_references: int = 5
    max_evalue: Optional[float] = None
    min_identity: float = 0.0


@dataclass
class ReferenceSearchConfig:
    """Runtime options for reference search engines."""

    reference_db: Path
    database_source: Optional[str] = None
    threads: int = 4
    filters: HitFilterConfig = field(default_factory=HitFilterConfig)
    output_path: Optional[Path] = None
    run_blast: bool = True
    blastn_binary: str = "blastn"
    blast_task: str = "blastn"


@dataclass(frozen=True)
class QueryReferenceMatch:
    """Retained alignments plus a primary (best-bitscore) hit for one query."""

    hits: tuple[QueryMatchConfidence, ...] = ()
    summary: Optional[QueryMatchConfidence] = None


@dataclass(frozen=True)
class SampleReferenceSearchResult:
    """Per-query matches plus sample-wide top reference ranking."""

    matches: dict[str, QueryReferenceMatch]
    top_references: tuple[SampleReferenceRank, ...] = ()


def _bitscore_of(hit: QueryMatchConfidence) -> float:
    if hit.bitscore is not None:
        return float(hit.bitscore)
    return float(hit.match_score)


def parse_blast_outfmt6_line(
        line: str,
        *,
        database_source: str,
) -> QueryMatchConfidence | None:
    """Parse one BLAST tabular (outfmt 6) line into a hit row."""
    if not line or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < _BLAST_MIN_COLS:
        return None
    try:
        pident = float(parts[_BLAST_PIDENT])
        aln_len = int(parts[_BLAST_LENGTH])
        q_start_1 = int(parts[_BLAST_QSTART])
        q_end_1 = int(parts[_BLAST_QEND])
        s_start_1 = int(parts[_BLAST_SSTART])
        s_end_1 = int(parts[_BLAST_SEND])
        evalue = float(parts[_BLAST_EVALUE])
        bitscore = float(parts[_BLAST_BITSCORE])
    except ValueError:
        return None
    if aln_len <= 0:
        return None

    # BLAST coords are 1-based inclusive; convert to 0-based half-open.
    q0 = min(q_start_1, q_end_1) - 1
    q1 = max(q_start_1, q_end_1)
    t0 = min(s_start_1, s_end_1) - 1
    t1 = max(s_start_1, s_end_1)
    strand = "+" if s_start_1 <= s_end_1 else "-"
    return QueryMatchConfidence(
        accession_id=parts[_BLAST_SSEQID],
        database_source=database_source,
        match_score=bitscore,
        sequence_identity=pident,
        query_id=parts[_BLAST_QSEQID],
        evalue=evalue,
        bitscore=bitscore,
        alignment_length=aln_len,
        query_start=q0,
        query_end=q1,
        target_start=t0,
        target_end=t1,
        functional_metadata={
            "strand": strand,
            "alignment_tool": "blastn",
        },
    )


def parse_blast_outfmt6(
        lines: Iterable[str],
        *,
        database_source: str,
) -> dict[str, list[QueryMatchConfidence]]:
    """Parse BLAST outfmt 6 lines into query_id -> hits (sorted by bitscore)."""
    out: dict[str, list[QueryMatchConfidence]] = {}
    for line in lines:
        hit = parse_blast_outfmt6_line(line, database_source=database_source)
        if hit is None or hit.query_id is None:
            continue
        out.setdefault(hit.query_id, []).append(hit)
    for rid in out:
        out[rid].sort(key=lambda h: (-_bitscore_of(h), -(h.sequence_identity or 0.0)))
    return out


def _pass_significance(hit: QueryMatchConfidence, cfg: HitFilterConfig) -> bool:
    if cfg.max_evalue is not None and hit.evalue is not None and hit.evalue > cfg.max_evalue:
        return False
    if hit.sequence_identity is not None and hit.sequence_identity < cfg.min_identity:
        return False
    if _bitscore_of(hit) < cfg.min_bitscore:
        return False
    return True


def filter_significant_hits(
        hits: Sequence[QueryMatchConfidence],
        config: HitFilterConfig | None = None,
) -> list[QueryMatchConfidence]:
    """Return hits that pass bitscore / optional identity / E-value thresholds."""
    cfg = config or HitFilterConfig()
    return [h for h in hits if _pass_significance(h, cfg)]


def _top_n_references(cfg: HitFilterConfig) -> int:
    return max(0, cfg.top_n_references)


def rank_sample_references(
        hits_by_query: Mapping[str, Sequence[QueryMatchConfidence]],
        *,
        top_n: int = 5,
) -> tuple[SampleReferenceRank, ...]:
    """
    Rank reference accessions by how many queries they cover.

    Primary key: number of distinct query reads with ≥1 hit.
    Tiebreak: sum of per-read best bitscores to that accession.
    """
    if top_n <= 0:
        return ()

    reads_by_acc: dict[str, set[str]] = defaultdict(set)
    best_bits_by_acc_read: dict[tuple[str, str], float] = {}
    for query_id, hits in hits_by_query.items():
        for hit in hits:
            acc = hit.accession_id
            reads_by_acc[acc].add(query_id)
            key = (acc, query_id)
            bits = _bitscore_of(hit)
            prev = best_bits_by_acc_read.get(key)
            if prev is None or bits > prev:
                best_bits_by_acc_read[key] = bits

    ranks: list[SampleReferenceRank] = []
    for acc, read_ids in reads_by_acc.items():
        sum_bits = sum(best_bits_by_acc_read[(acc, rid)] for rid in read_ids)
        ranks.append(
            SampleReferenceRank(
                accession_id=acc,
                n_reads=len(read_ids),
                sum_bitscore=sum_bits,
            )
        )
    ranks.sort(key=lambda r: (-r.n_reads, -r.sum_bitscore, r.accession_id))
    return tuple(ranks[:top_n])


def make_primary_summary(
        selected_hits: Sequence[QueryMatchConfidence],
) -> QueryMatchConfidence | None:
    """Return the best-bitscore retained hit (primary alignment)."""
    if not selected_hits:
        return None
    return max(
        selected_hits,
        key=lambda h: (_bitscore_of(h), h.sequence_identity or 0.0),
    )


def build_query_matches(
        hits_by_query: Mapping[str, Sequence[QueryMatchConfidence]],
        config: HitFilterConfig | None = None,
) -> SampleReferenceSearchResult:
    """
    Bitscore-filter hits, pick sample top-N references, attach per-query hits.

    Ranking: references covering the most distinct reads win; bitscore sum
    breaks ties. Only hits to those top references are retained per read.
    """
    cfg = config or HitFilterConfig()

    filtered: dict[str, list[QueryMatchConfidence]] = {}
    for query_id, hits in hits_by_query.items():
        kept = filter_significant_hits(hits, cfg)
        filtered[query_id] = kept

    top_refs = rank_sample_references(
        filtered, top_n=_top_n_references(cfg)
    )
    keep_acc = {r.accession_id for r in top_refs}

    out: dict[str, QueryReferenceMatch] = {}
    for query_id, hits in filtered.items():
        selected = [h for h in hits if h.accession_id in keep_acc]
        selected.sort(key=lambda h: (-_bitscore_of(h), -(h.sequence_identity or 0.0)))
        out[query_id] = QueryReferenceMatch(
            hits=tuple(selected),
            summary=make_primary_summary(selected),
        )
    return SampleReferenceSearchResult(matches=out, top_references=top_refs)


def _seq_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".fq", ".fastq"}:
        return "fastq"
    return "fasta"


def write_engineered_query_fasta(
        query_ids: set[str],
        source_reads: str | Path,
        *,
        output_path: str | Path | None = None,
) -> Path | None:
    """
    Write a temporary FASTA of engineered read IDs for blastn.

    Returns ``None`` when ``query_ids`` is empty.
    """
    if not query_ids:
        return None

    source = Path(source_reads)
    fmt = _seq_format(source)
    wanted = set(query_ids)
    records = []
    for record in SeqIO.parse(str(source), fmt):
        if record.id in wanted:
            records.append(record)
            wanted.discard(record.id)
            if not wanted:
                break

    if not records:
        logger.warning(
            "No engineered read IDs found in %s (%d requested)",
            source,
            len(query_ids),
        )
        return None

    if output_path is None:
        fd, name = tempfile.mkstemp(prefix="plasmidScreen_eng_", suffix=".fasta")
        os.close(fd)
        out = Path(name)
    else:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)

    SeqIO.write(records, str(out), "fasta")
    logger.info(
        "Wrote %d engineered quer(ies) for blastn to %s",
        len(records),
        out,
    )
    return out


class ReferenceSearchEngine(ABC):
    """Abstract interface for batch reference-database lookup."""

    @abstractmethod
    def search_batch(
            self,
            query_fasta: str | Path,
            *,
            query_ids: Sequence[str] | None = None,
    ) -> SampleReferenceSearchResult:
        """Query the reference database for sequences in ``query_fasta``."""


class BlastnReferenceSearchEngine(ReferenceSearchEngine):
    """
    Nucleotide search via local ``blastn``.

    Aligns engineered query reads to a FASTA subject or BLAST DB, keeps hits
    above ``min_bitscore``, then retains only the sample's top-N references
    (by # reads covered, bitscore sum as tiebreak).
    """

    def __init__(self, config: ReferenceSearchConfig) -> None:
        self.config = config
        self.database_source = (
                config.database_source
                or f"blastn:{Path(config.reference_db).name}"
        )

    def search_batch(
            self,
            query_fasta: str | Path,
            *,
            query_ids: Sequence[str] | None = None,
    ) -> SampleReferenceSearchResult:
        t0 = time.perf_counter()
        lines = self._resolve_alignment_lines(query_fasta)
        result = self._matches_from_lines(lines)
        matches = result.matches

        if query_ids is not None:
            ordered = {
                qid: matches.get(qid, QueryReferenceMatch(hits=()))
                for qid in query_ids
            }
        else:
            ordered = matches

        elapsed = time.perf_counter() - t0
        n_with_hits = sum(1 for m in ordered.values() if m.hits)
        n_hits = sum(len(m.hits) for m in ordered.values())
        logger.info(
            "blastn reference search finished in %.3fs: %d queries, "
            "%d with ≥1 retained hit, %d total hits, top_refs=%s "
            "(database=%s, min_bitscore=%.1f, top_n=%d)",
            elapsed,
            len(ordered),
            n_with_hits,
            n_hits,
            [r.accession_id for r in result.top_references],
            self.database_source,
            self.config.filters.min_bitscore,
            _top_n_references(self.config.filters),
        )
        return SampleReferenceSearchResult(
            matches=ordered, top_references=result.top_references
        )

    def search_from_blast_lines(
            self,
            lines: Iterable[str],
            *,
            query_ids: Sequence[str] | None = None,
    ) -> SampleReferenceSearchResult:
        """Parse/filter precomputed BLAST outfmt 6 lines (no subprocess)."""
        result = self._matches_from_lines(lines)
        if query_ids is None:
            return result
        ordered = {
            qid: result.matches.get(qid, QueryReferenceMatch(hits=()))
            for qid in query_ids
        }
        return SampleReferenceSearchResult(
            matches=ordered, top_references=result.top_references
        )

    def _resolve_alignment_lines(self, query_fasta: str | Path) -> list[str]:
        cfg = self.config
        if not cfg.run_blast:
            if cfg.output_path is None:
                raise ValueError(
                    "output_path is required when run_blast=False "
                    "(precomputed blastn outfmt 6 for reference search)."
                )
            path = Path(cfg.output_path)
            if not path.is_file():
                raise FileNotFoundError(f"blastn output not found: {path}")
            logger.info("Loading reference-search blastn output from %s", path)
            return path.read_text(encoding="utf-8").splitlines()

        return run_blastn_alignment(
            Path(query_fasta),
            Path(cfg.reference_db),
            threads=cfg.threads,
            blastn_binary=cfg.blastn_binary,
            blast_task=cfg.blast_task,
            max_evalue=cfg.filters.max_evalue,
            output_path=cfg.output_path,
        )

    def _matches_from_lines(
            self,
            lines: Iterable[str],
    ) -> SampleReferenceSearchResult:
        parsed = parse_blast_outfmt6(lines, database_source=self.database_source)
        return build_query_matches(parsed, self.config.filters)


def _reference_is_fasta(path: Path) -> bool:
    return path.suffix.lower() in _FASTA_SUFFIXES


def run_blastn_alignment(
        reads_path: Path,
        reference_db: Path,
        *,
        threads: int = 4,
        blastn_binary: str = "blastn",
        blast_task: str = "blastn",
        max_evalue: float | None = None,
        output_path: str | Path | None = None,
) -> list[str]:
    """
    Run ``blastn`` and return outfmt 6 tabular lines.

    ``reference_db`` may be a nucleotide FASTA (``-subject``) or a BLAST DB
    prefix (``-db``).
    """
    reads_path = Path(reads_path)
    reference_db = Path(reference_db)
    cmd: list[str] = [
        blastn_binary,
        "-task",
        blast_task,
        "-query",
        str(reads_path),
        "-outfmt",
        "6",
        "-num_threads",
        str(threads),
        "-dust",
        "no",
    ]
    if max_evalue is not None:
        cmd.extend(["-evalue", str(max_evalue)])
    if _reference_is_fasta(reference_db):
        cmd.extend(["-subject", str(reference_db)])
    else:
        cmd.extend(["-db", str(reference_db)])

    logger.info("Running blastn reference search: %s", " ".join(cmd))
    t0 = time.perf_counter()
    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd.extend(["-out", str(out)])
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.perf_counter() - t0
        logger.info("blastn wrote %s in %.3fs", out, elapsed)
        return out.read_text(encoding="utf-8").splitlines()

    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    lines = proc.stdout.splitlines()
    logger.info("blastn returned %d hit lines in %.3fs", len(lines), elapsed)
    return lines


def enrich_screen_result(
        screen_result: ScreenResult,
        matches_by_read: Mapping[str, QueryReferenceMatch | Sequence[QueryMatchConfidence]]
                         | SampleReferenceSearchResult,
        *,
        database_source: str | None = None,
        searched: bool = True,
        annotate_only_ids: set[str] | None = None,
        top_references: tuple[SampleReferenceRank, ...] | None = None,
) -> ScreenResult:
    """
    Attach reference-match hits onto a finished :class:`ScreenResult`.

    When ``annotate_only_ids`` is set (engineered read IDs), only those reads
    receive match notes / summaries; other reads are left unchanged.

    Mapping is informational only: engineered / natural overall labels are never
    rewritten based on database match success or failure.
    """
    ranks = top_references
    if isinstance(matches_by_read, SampleReferenceSearchResult):
        ranks = matches_by_read.top_references if ranks is None else ranks
        matches_by_read = matches_by_read.matches

    annotate = annotate_only_ids
    new_per_read = []
    for detail in screen_result.per_read:
        if annotate is not None and detail.read_id not in annotate:
            new_per_read.append(detail)
            continue

        raw = matches_by_read.get(detail.read_id)
        if isinstance(raw, QueryReferenceMatch):
            hits = tuple(raw.hits)
            summary = raw.summary
        elif raw is None:
            hits = ()
            summary = None
        else:
            hits = tuple(raw)
            summary = make_primary_summary(hits)

        has_match = bool(hits) and summary is not None
        note: str | None = None
        if searched and not has_match:
            note = NO_DATABASE_MATCH_NOTE

        new_per_read.append(
            replace(
                detail,
                database_match_note=note,
                database_match_summary=summary if has_match else None,
                database_match_hits=hits if has_match else (),
            )
        )

    source = database_source or screen_result.reference_database_source
    if source is None:
        for detail in new_per_read:
            if detail.database_match_summary is not None:
                source = detail.database_match_summary.database_source
                break

    if ranks is None:
        hits_by_query = {
            d.read_id: d.database_match_hits
            for d in new_per_read
            if d.database_match_hits
        }
        ranks = rank_sample_references(hits_by_query, top_n=5)

    # Mapping pileup summaries are viz-only; see visualize_engineered_mappings.
    return replace(
        screen_result,
        per_read=new_per_read,
        reference_database_source=source,
        top_engineered_references=ranks,
        engineered_mapping_hits=(),
    )


def enrich_screen_result_with_engine(
        screen_result: ScreenResult,
        engine: ReferenceSearchEngine,
        query_fasta: str | Path,
) -> ScreenResult | None:
    """
    Align engineered reads with blastn (or legacy engine) and join hits onto
    ``screen_result``.

    Natural-only reads are left without a database-match annotation.
    """
    engineered_ids = set(screen_result.engineered_read_ids)
    if not engineered_ids:
        logger.info("No engineered reads; skipping reference search.")
        return None

    subset = write_engineered_query_fasta(engineered_ids, query_fasta)
    if subset is None:
        logger.warning(
            "Could not fetch engineered sequences from %s; skipping reference search.",
            query_fasta,
        )
        return None

    try:
        search = engine.search_batch(
            subset, query_ids=sorted(engineered_ids)
        )

        source = getattr(engine, "database_source", None)
        return enrich_screen_result(
            screen_result,
            search,
            database_source=source,
            searched=True,
            annotate_only_ids=engineered_ids,
        )
    finally:
        try:
            subset.unlink(missing_ok=True)
        except OSError:
            logger.warning("Could not unlink temporary query FASTA %s", subset)
