"""
PlasmidScreen public library API.

All screening assumes codon reference data was built offline (airgapped-safe).
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

from plasmidScreen.lib.codon_usage_build import build_codon_reference, default_reference_taxids
from plasmidScreen.lib.codon_usage_db import (
    CodonUsageStore,
    default_codon_usage_dir,
    taxids_from_kraken_output,
)
from plasmidScreen.lib.models import (
    BuildCodonReferenceResult,
    CodonAdaptationResult,
    ScreenResult,
)
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
    write_codon_adaptation_tsv,
)
from plasmidScreen.src.plasmidScreen import Workflow

__all__ = [
    "analyze_codon_adaptation",
    "build_codon_reference",
    "default_reference_taxids",
    "default_codon_usage_dir",
    "run_screen",
    "taxids_from_kraken_output",
    "write_codon_adaptation_tsv",
    "BuildCodonReferenceResult",
    "CodonAdaptationResult",
    "CodonUsageStore",
    "ScreenResult",
]


def run_screen(
    fasta_file: str | Path,
    output_report_path: str | Path,
    kraken_raw_output: str | Path,
    kraken_db: str | Path,
    *,
    threads: int = 4,
    window_size: int = 200,
    engineered_kmer_threshold: int = 25,
    codon_usage_dir: str | Path | None = None,
    run_kraken: bool = True,
    run_codon_usage: bool = True,
    codon_usage_output_path: str | Path | None = None,
    kmer_len: int = 35,
    codon_cai_engineered_threshold: float | None = None,
) -> ScreenResult:
    """
    Run engineered k-mer screening and optional codon adaptation analysis.

    Parameters
    ----------
    codon_usage_dir
        Pre-built reference directory (codon_tables.json). Required when
        run_codon_usage is True.
    run_kraken
        If False, kraken_raw_output must already exist (airgapped Kraken step done).
    """
    workflow = Workflow(
        str(fasta_file),
        str(output_report_path),
        str(kraken_db),
        threads,
        str(kraken_raw_output),
        window_size=window_size,
        engineered_kmer_threshold=engineered_kmer_threshold,
        codon_usage_output_path=str(codon_usage_output_path) if codon_usage_output_path else None,
        codon_usage_dir=str(codon_usage_dir) if codon_usage_dir else None,
        run_kraken=run_kraken,
        run_codon_usage=run_codon_usage,
        kmer_len=kmer_len,
        codon_cai_engineered_threshold=codon_cai_engineered_threshold,
    )
    return workflow.run()
