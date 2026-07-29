"""Structured results returned by the PlasmidScreen library API."""
from __future__ import annotations
from enum import Enum
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional

# Provenance tag when a queried sequence has no significant reference hits.
NO_DATABASE_MATCH_NOTE = "No database match found"


class SupportDataTypes(Enum):
    """minimap2 ``-x`` presets for reference search on engineered reads."""

    LONG_READ_ONT = "map-ont"
    LONG_READ_PB = "map-pb"
    SHORT_READ = "sr"
    ASM = "asm5"


@dataclass(frozen=True)
class ReferenceHit:
    """
    A single verified hit from a reference sequence database search.

    ``match_score`` is backend-specific: for DIAMOND/BLAST it is typically the
    bit-score (higher is better); for vector search it may be cosine similarity
    or negative L2 distance. Alignment E-value and coverage details belong in
    ``functional_metadata`` (and optional dedicated fields) for full provenance.

    ``query_start`` / ``query_end`` are 0-based half-open coordinates on the
    query so mosaic / multi-accession tiling can cover distinct query regions.
    """

    accession_id: str
    database_source: str
    match_score: float
    sequence_identity: Optional[float] = None
    functional_metadata: dict[str, Any] = field(default_factory=dict)
    query_id: Optional[str] = None
    evalue: Optional[float] = None
    bitscore: Optional[float] = None
    alignment_length: Optional[int] = None
    query_start: Optional[int] = None
    query_end: Optional[int] = None
    mapq: Optional[int] = None

    @property
    def query_span_bp(self) -> int:
        """Aligned query bases (0 if coordinates are missing)."""
        if self.query_start is None or self.query_end is None:
            return 0
        return max(0, self.query_end - self.query_start)


@dataclass(frozen=True)
class QueryMatchConfidence:
    """
    Aggregate confidence that a query is explained by its retained reference hits.

    Inspired by PlasmidHawk CORRECT-mode uniqueness weighting: fragments
    (hits) shared by many competing accessions contribute less. The published
    PlasmidHawk analysis related uniqueness to scores via linear regression; we
    mirror that strategy as a linear combination of coverage, identity, and
    uniqueness features scaled to ``[0, 1]``.
    """

    confidence: float
    query_coverage: float
    covered_bp: int
    query_length_bp: int
    mean_identity: float
    uniqueness: float
    plasmidhawk_score: float
    n_accessions: int
    n_hits: int
    features: dict[str, float] = field(default_factory=dict)


@dataclass
class CodonAdaptationResult:
    labels: list[CodonAdaptationRead] | None = field(default=list)

    @property
    def natural_read_ids(self) -> set[str]:
        return {r.read_id for r in self.labels if r.label == "Natural"}

    @property
    def engineered_read_ids(self) -> set[str]:
        return {r.read_id for r in self.labels if r.label == "Synthetic"}


@dataclass(frozen=True)
class CodonAdaptationRead:
    """Per-read codon adaptation scores from DIAMOND ORFs and CSDB reference weights."""

    read_id: str
    #label: Literal["Natural", "Synthetic"]
    cds_strand: str  # "+" or "-" on the read
    cds_start: int  # 0-based start (half-open interval with cds_end)
    cds_end: int
    host_taxid: str  # NCBI taxid from DIAMOND staxids (majority over hits)
    reference_taxid: Optional[str]  # CSDB table used after lineage resolution
    cds_len_bp: int
    cai_vs_host: Optional[float]  # Sharp & Li CAI vs host reference (0–1)
    host_taxid_method: Optional[str] = None  # e.g. "diamond"


@dataclass(frozen=True)
class ReadEngineeringLabel:
    read_id: str
    label: Literal["Natural", "Synthetic"]


@dataclass
class EngineeredScanResult:
    labels: list[ReadEngineeringLabel] = field(default_factory=list)
    synthetic_count: int = 0
    natural_count: int = 0

    @property
    def natural_read_ids(self) -> set[str]:
        return {r.read_id for r in self.labels if r.label == "Natural"}

    @property
    def engineered_read_ids(self) -> set[str]:
        return {r.read_id for r in self.labels if r.label == "Synthetic"}

    @property
    def any_synthetic(self) -> bool:
        return self.synthetic_count > 0


@dataclass
class BuildCodonReferenceResult:
    """Outcome of an offline codon reference build."""

    data_dir: Path
    taxids_requested: list[str]
    taxids_added: list[str]
    taxids_skipped: list[str]
    taxids_failed: list[str]


def compute_engineered_overall(
        *,
        engineered_by_kmer_scan: bool,
        engineered_by_codon_cai: bool | None,
) -> bool:
    """
    Combined engineered call for a read using thresholds active in ``run_screen``.

    Engineered if the k-mer scan flagged Synthetic, or if codon CAI flagging is
    enabled (threshold set) and CAI is below that threshold.
    """
    if engineered_by_kmer_scan:
        return True
    return engineered_by_codon_cai is True


@dataclass
class ScreenResult:
    """
    Full screening run result (engineered k-mer scan + optional codon usage).

    Use ``per_read`` for per-read ``engineered_overall`` / ``overall_label`` and
    ``overall_synthetic_count`` for a run total that respects both k-mer and codon thresholds.

    When reference search is enabled, ``database_hits`` holds a flat aggregation of
    significant hits across all reads; per-read hits and match notes live on
    :class:`ReadFlagDetail`.
    """

    engineered_scan: EngineeredScanResult
    codon_adaptation: CodonAdaptationResult
    engineered_kmer_threshold: int
    engineered_kmer_window_size: int
    codon_cai_engineered_threshold: float
    per_read: list["ReadFlagDetail"] = field(default_factory=list)
    engineered_report_path: Optional[Path] = None
    codon_usage_report_path: Optional[Path] = None
    diamond_output_path: Optional[Path] = None
    database_hits: list[ReferenceHit] = field(default_factory=list)
    reference_database_source: Optional[str] = None

    @property
    def overall_engineered_read_count(self) -> int:
        """Reads classified as engineered under the combined k-mer + codon rules."""
        return sum(1 for r in self.per_read if r.engineered_overall)

    @property
    def overall_natural_read_count(self) -> int:
        return sum(1 for r in self.per_read if not r.engineered_overall)

    @property
    def engineered_read_ids(self) -> set[str]:
        return {r.read_id for r in self.per_read if r.engineered_overall}

    @property
    def natural_read_ids_overall(self) -> set[str]:
        return {r.read_id for r in self.per_read if not r.engineered_overall}


@dataclass(frozen=True)
class ReadFlagDetail:
    """Per-read summary of which method(s) flagged engineered and the overall call."""

    read_id: str
    kmer_label: Literal["Natural", "Synthetic"]
    engineered_by_kmer_scan: bool
    engineered_overall: bool
    overall_label: Literal["Natural", "Synthetic"]
    engineered_kmer_max_in_window: Optional[int] = None
    engineered_kmer_threshold: Optional[int] = None
    engineered_kmer_window_size: Optional[int] = None
    cai_vs_host: Optional[float] = None
    engineered_by_codon_cai: Optional[bool] = None
    codon_cai_threshold: Optional[float] = None
    database_hits: tuple[ReferenceHit, ...] = ()
    database_match_note: Optional[str] = None
    database_match_confidence: Optional[float] = None
    database_query_coverage: Optional[float] = None
    database_match_summary: Optional[QueryMatchConfidence] = None

    @property
    def engineered_methods(self) -> list[str]:
        methods: list[str] = []
        if self.engineered_by_kmer_scan:
            methods.append("engineered_kmer_scan")
        if self.engineered_by_codon_cai:
            methods.append("codon_cai")
        return methods

    @property
    def engineered_any(self) -> bool:
        """Alias for :attr:`engineered_overall` (combined threshold decision)."""
        return self.engineered_overall

    @property
    def has_database_match(self) -> bool:
        """True when at least one significant reference hit was retained."""
        return bool(self.database_hits)
