"""Structured results returned by the PlasmidScreen library API."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional
from enum import Enum

MEMORY_CONFIG = {"LOW": 0, "MEDIUM": 1, "HIGH": 2}

BACKGROUND_RATE = 0.0003

# Provenance tag when a queried sequence has no significant reference hits.
NO_DATABASE_MATCH_NOTE = "No database match found"


class MEMORY_MODE(Enum):
    LOW = "LOW"
    MEDIUM = "MEDIUM"
    HIGH = "HIGH"


class SupportDataTypes(Enum):
    """minimap2 ``-x`` presets for reference search on engineered reads."""

    LONG_READ_ONT = "map-ont"
    LONG_READ_PB = "map-pb"
    SHORT_READ = "sr"
    ASM = "asm5"


@dataclass(frozen=True)
class QueryMatchConfidence:
    """
    Unified reference-match record: one alignment *or* a scored per-query summary.

    Hit-level fields describe a single alignment; aggregate fields
    (``confidence``, ``query_coverage``, …) describe scoring over retained hits.
    """

    accession_id: str
    database_source: str
    match_score: float
    sequence_identity: Optional[float] = None
    query_id: Optional[str] = None
    evalue: Optional[float] = None
    bitscore: Optional[float] = None
    alignment_length: Optional[int] = None
    query_start: Optional[int] = None
    query_end: Optional[int] = None
    mapq: Optional[int] = None
    target_start: Optional[int] = None
    target_end: Optional[int] = None
    target_length: Optional[int] = None
    functional_metadata: dict[str, Any] = field(default_factory=dict)
    confidence: float = 0.0
    query_coverage: float = 0.0
    covered_bp: int = 0
    query_length_bp: int = 0
    mean_identity: float = 0.0
    uniqueness: float = 0.0
    plasmidhawk_score: float = 0.0
    n_accessions: int = 1
    n_hits: int = 1
    features: dict[str, float] = field(default_factory=dict)

    @property
    def query_span_bp(self) -> int:
        if self.query_start is None or self.query_end is None:
            return 0
        return max(0, self.query_end - self.query_start)

    @property
    def target_span_bp(self) -> int:
        if self.target_start is None or self.target_end is None:
            return 0
        return max(0, self.target_end - self.target_start)


@dataclass
class CodonAdaptationResult:
    """Per-read CAI scores (no engineered threshold applied here).

    For reads flagged engineered by CAI, use
    :attr:`ScreenResult.codon_cai_engineered_read_ids`.
    """

    labels: list[CodonAdaptationRead] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.labels)

    def __iter__(self):
        return iter(self.labels)

    def __getitem__(self, index: int) -> CodonAdaptationRead:
        return self.labels[index]

    @property
    def scored_read_ids(self) -> set[str]:
        """Read IDs that received a codon-adaptation result (host + ORF + CAI attempt)."""
        return {r.read_id for r in self.labels}


@dataclass(frozen=True)
class CodonAdaptationRead:
    """Per-read codon adaptation scores from DIAMOND ORFs and CSDB reference weights."""

    read_id: str
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
    eng_kmer_coverage: float

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
        engineered_by_perseus: bool | None = None,
) -> bool:
    """
    Combined engineered call for a read using thresholds active in ``run_screen``.

    Engineered if the k-mer scan flagged Synthetic, Perseus assigned taxid 32630
    on a Natural read, or codon CAI flagging is enabled and CAI is above threshold.
    """
    if engineered_by_kmer_scan:
        return True
    if engineered_by_perseus is True:
        return True
    return engineered_by_codon_cai is True


@dataclass
class ScreenResult:
    """
    Full screening run result (engineered k-mer scan + optional codon usage).

    Use ``per_read`` for per-read ``engineered_overall`` / ``overall_label`` and
    ``overall_synthetic_count`` for a run total that respects both k-mer and codon thresholds.
    """

    engineered_scan: EngineeredScanResult
    codon_adaptation: CodonAdaptationResult
    engineered_kmer_threshold: int
    engineered_kmer_window_size: int
    codon_cai_engineered_threshold: float | None
    per_read: list["ReadFlagDetail"] = field(default_factory=list)
    engineered_report_path: Optional[Path] = None
    codon_usage_report_path: Optional[Path] = None
    diamond_output_path: Optional[Path] = None
    perseus_output_path: Optional[Path] = None

    @property
    def overall_engineered_read_count(self) -> int:
        """Reads classified as engineered under the combined k-mer + codon rules."""
        return sum(1 for r in self.per_read if r.engineered_overall)

    @property
    def overall_natural_read_count(self) -> int:
        return sum(1 for r in self.per_read if not r.engineered_overall)

    # README / API aliases
    @property
    def overall_synthetic_count(self) -> int:
        return self.overall_engineered_read_count

    @property
    def overall_natural_count(self) -> int:
        return self.overall_natural_read_count

    @property
    def engineered_read_ids(self) -> set[str]:
        """Reads engineered overall (k-mer, Perseus, and/or CAI)."""
        return {r.read_id for r in self.per_read if r.engineered_overall}

    @property
    def natural_read_ids_overall(self) -> set[str]:
        return {r.read_id for r in self.per_read if not r.engineered_overall}

    @property
    def codon_cai_engineered_read_ids(self) -> set[str]:
        """Reads with CAI above ``codon_cai_engineered_threshold`` (empty if unset)."""
        return {r.read_id for r in self.per_read if r.engineered_by_codon_cai}

    @property
    def perseus_engineered_read_ids(self) -> set[str]:
        """Reads Perseus assigned to synthetic construct (taxid 32630)."""
        return {r.read_id for r in self.per_read if r.engineered_by_perseus}

    @property
    def per_read_engineered(self) -> list["ReadFlagDetail"]:
        return [r for r in self.per_read if r.overall_label == "Synthetic"]

    @property
    def per_read_natural(self) -> list["ReadFlagDetail"]:
        return [r for r in self.per_read if r.overall_label == "Natural"]


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
    engineered_by_perseus: Optional[bool] = None
    perseus_taxid: Optional[int] = None
    p_value: Optional[float] = None
    eng_kmer_coverage: Optional[float] = None
    database_match_note: Optional[str] = None
    database_match_confidence: Optional[float] = None
    database_query_coverage: Optional[float] = None
    database_match_summary: Optional[QueryMatchConfidence] = None
    database_match_hits: tuple[QueryMatchConfidence, ...] = ()

    @property
    def engineered_methods(self) -> list[str]:
        methods: list[str] = []
        if self.engineered_by_kmer_scan:
            methods.append("engineered_kmer_scan")
        if self.engineered_by_perseus:
            methods.append("perseus")
        if self.engineered_by_codon_cai:
            methods.append("codon_cai")
        return methods

    @property
    def has_database_match(self) -> bool:
        return self.database_match_summary is not None
