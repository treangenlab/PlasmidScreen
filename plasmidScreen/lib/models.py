"""Structured results returned by the PlasmidScreen library API."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional


@dataclass(frozen=True)
class CodonAdaptationResult:
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


@dataclass
class EngineeredScanResult:
    labels: list[ReadEngineeringLabel] = field(default_factory=list)
    synthetic_count: int = 0
    natural_count: int = 0

    @property
    def natural_read_ids(self) -> set[str]:
        return {r.read_id for r in self.labels if r.label == "Natural"}

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
    """

    engineered_scan: EngineeredScanResult
    codon_adaptation: list[CodonAdaptationResult] = field(default_factory=list)
    per_read: list["ReadFlagDetail"] = field(default_factory=list)
    engineered_report_path: Optional[Path] = None
    codon_usage_report_path: Optional[Path] = None
    diamond_output_path: Optional[Path] = None
    engineered_kmer_threshold: int = 25
    engineered_kmer_window_size: int = 200
    codon_cai_engineered_threshold: Optional[float] = None

    @property
    def overall_synthetic_count(self) -> int:
        """Reads classified as engineered under the combined k-mer + codon rules."""
        return sum(1 for r in self.per_read if r.engineered_overall)

    @property
    def overall_natural_count(self) -> int:
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
