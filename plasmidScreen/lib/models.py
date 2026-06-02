"""Structured results returned by the PlasmidScreen library API."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional


@dataclass(frozen=True)
class CodonAdaptationResult:
    read_id: str
    cds_strand: str
    cds_start: int
    cds_end: int
    host_taxid: str
    cds_len_bp: int
    cai_vs_host: Optional[float]
    host_taxid_method: Optional[str] = None


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


@dataclass
class ScreenResult:
    """Full screening run result (engineered k-mer scan + optional codon usage)."""

    engineered_scan: EngineeredScanResult
    codon_adaptation: list[CodonAdaptationResult] = field(default_factory=list)
    per_read: list["ReadFlagDetail"] = field(default_factory=list)
    engineered_report_path: Optional[Path] = None
    codon_usage_report_path: Optional[Path] = None
    diamond_output_path: Optional[Path] = None


@dataclass(frozen=True)
class ReadFlagDetail:
    """Per-read summary of which method(s) flagged engineered."""

    read_id: str
    kmer_label: Literal["Natural", "Synthetic"]
    engineered_by_kmer_scan: bool
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
        return bool(self.engineered_by_kmer_scan or self.engineered_by_codon_cai)
