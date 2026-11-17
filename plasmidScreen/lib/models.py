"""Structured results returned by the PlasmidScreen library API."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional
from enum import Enum

MEMORY_CONFIG = {"LOW": 0, "MEDIUM": 1, "HIGH": 2}

BACKGROUND_RATE = 0.001
# Reject window-threshold Synthetic calls when the binomial p-value exceeds this.
PVALUE_ALPHA = 0.05

# Provenance tag when a queried sequence has no significant reference hits.
NO_DATABASE_MATCH_NOTE = "No database match found"


class MEMORY_MODE(Enum):
    LOW = "LOW"
    MEDIUM = "MEDIUM"
    HIGH = "HIGH"




@dataclass(frozen=True)
class QueryMatchConfidence:
    """One BLAST (or equivalent) alignment of a query read to a reference."""

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
    target_start: Optional[int] = None
    target_end: Optional[int] = None
    target_length: Optional[int] = None
    functional_metadata: dict[str, Any] = field(default_factory=dict)

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


@dataclass(frozen=True)
class MappedQueryInterval:
    """One query alignment interval on a reference plasmid (pileup row)."""

    query_id: str
    target_start: int  # 0-based inclusive
    target_end: int  # 0-based exclusive
    sequence_identity: Optional[float] = None
    match_score: Optional[float] = None
    overlapping_features: tuple[str, ...] = ()


@dataclass(frozen=True)
class ReferenceFeatureInterval:
    """One GenBank feature interval shown on the plasmid map."""

    label: str
    feature_type: str
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    strand: int = 0


@dataclass(frozen=True)
class SampleReferenceRank:
    """One reference sequence ranked by how many engineered reads it covers."""

    accession_id: str
    n_reads: int
    sum_bitscore: float


@dataclass(frozen=True)
class EngineeredMappingHit:
    """
    Per-plasmid mapping summary for engineered reads (mirrors the pileup figure).

    Built from retained BLAST hits (and optional GenBank features). Use
    :attr:`alignments` for where each query maps; coverage fields summarize
    the aggregated depth track used when many reads hit the same plasmid.
    """

    accession_id: str
    target_length: int
    n_queries: int
    n_alignments: int
    query_ids: tuple[str, ...]
    alignments: tuple[MappedQueryInterval, ...]
    features: tuple[ReferenceFeatureInterval, ...] = ()
    covered_bp: int = 0
    mean_coverage: float = 0.0
    max_coverage: int = 0


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
        # engineered_by_perseus: bool | None = None,  # Perseus OR disabled
) -> bool:
    """
    Combined engineered call for a read using thresholds active in ``run_screen``.

    Engineered if the k-mer scan flagged Synthetic, or codon CAI flagging is
    enabled and CAI is above threshold.
    """
    if engineered_by_kmer_scan:
        return True
    # if engineered_by_perseus is True:
    #     return True
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
    # Perseus disabled — field kept commented for re-enable.
    # perseus_output_path: Optional[Path] = None
    reference_database_source: Optional[str] = None
    # Sample-wide top reference sequences (BLAST bitscore filter + read coverage rank).
    top_engineered_references: tuple[SampleReferenceRank, ...] = ()
    # Filled only when visualization mode runs (:func:`visualize_engineered_mappings`).
    engineered_mapping_hits: tuple[EngineeredMappingHit, ...] = ()
    # Full viz artifact (pileup paths + mapping summaries); set when viz is enabled.
    visualization: Optional["VisualizationReport"] = None

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
        """Reads engineered overall (k-mer and/or CAI)."""
        return {r.read_id for r in self.per_read if r.engineered_overall}

    @property
    def natural_read_ids_overall(self) -> set[str]:
        return {r.read_id for r in self.per_read if not r.engineered_overall}

    @property
    def codon_cai_engineered_read_ids(self) -> set[str]:
        """Reads with CAI above ``codon_cai_engineered_threshold`` (empty if unset)."""
        return {r.read_id for r in self.per_read if r.engineered_by_codon_cai}

    # Perseus disabled
    # @property
    # def perseus_engineered_read_ids(self) -> set[str]:
    #     """Reads Perseus assigned to synthetic construct (taxid 32630)."""
    #     return {r.read_id for r in self.per_read if r.engineered_by_perseus}

    @property
    def per_read_engineered(self) -> list["ReadFlagDetail"]:
        return [r for r in self.per_read if r.overall_label == "Synthetic"]

    @property
    def per_read_natural(self) -> list["ReadFlagDetail"]:
        return [r for r in self.per_read if r.overall_label == "Natural"]

    @property
    def database_match_summaries(self) -> list[QueryMatchConfidence]:
        """Primary per-read database-match summaries (reads with a retained hit)."""
        return [
            r.database_match_summary
            for r in self.per_read
            if r.database_match_summary is not None
        ]


@dataclass(frozen=True)
class VisualizationReport:
    """
    Optional viz-mode artifact (pileups + mapping summaries).

    Produced by :func:`plasmidScreen.api.visualize_engineered_mappings` and attached
    to ``ScreenResult.visualization`` when visualization is enabled. Also returned
    directly from the viz API helpers. Primary BLAST / feature columns live on
    ``ReadFlagDetail`` and the engineered TSV.
    """

    output_dir: Path
    pileup_paths: tuple[Path, ...] = ()
    engineered_mapping_hits: tuple[EngineeredMappingHit, ...] = ()
    top_engineered_references: tuple[SampleReferenceRank, ...] = ()


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
    # Perseus disabled
    # engineered_by_perseus: Optional[bool] = None
    # perseus_taxid: Optional[int] = None
    p_value: Optional[float] = None
    eng_kmer_coverage: Optional[float] = None
    database_match_note: Optional[str] = None
    database_match_summary: Optional[QueryMatchConfidence] = None
    database_match_hits: tuple[QueryMatchConfidence, ...] = ()
    # GenBank labels overlapped by the primary hit (filled during map-hits).
    overlapping_features: tuple[str, ...] = ()

    @property
    def engineered_methods(self) -> list[str]:
        methods: list[str] = []
        if self.engineered_by_kmer_scan:
            methods.append("engineered_kmer_scan")
        # if self.engineered_by_perseus:
        #     methods.append("perseus")
        if self.engineered_by_codon_cai:
            methods.append("codon_cai")
        return methods

    @property
    def has_database_match(self) -> bool:
        return self.database_match_summary is not None

    @property
    def matched_plasmid(self) -> Optional[str]:
        """Primary BLAST accession (best bitscore), if mapped."""
        summary = self.database_match_summary
        return None if summary is None else summary.accession_id

    @property
    def match_target_start(self) -> Optional[int]:
        summary = self.database_match_summary
        return None if summary is None else summary.target_start

    @property
    def match_target_end(self) -> Optional[int]:
        summary = self.database_match_summary
        return None if summary is None else summary.target_end

    @property
    def match_sequence_identity(self) -> Optional[float]:
        summary = self.database_match_summary
        return None if summary is None else summary.sequence_identity

    @property
    def match_score(self) -> Optional[float]:
        summary = self.database_match_summary
        return None if summary is None else summary.match_score
