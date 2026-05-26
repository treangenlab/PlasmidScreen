"""
PlasmidScreen: engineered DNA detection and codon adaptation analysis.

Library usage (airgapped screening requires a pre-built codon reference):

    from plasmidScreen import build_codon_reference, analyze_codon_adaptation, run_screen

    # 1. On a networked machine — build reference once
    build_codon_reference("codon_usage/", taxids=["9606", "511145"])

    # 2. On airgapped machine — analyze
    results = analyze_codon_adaptation("reads.fa", "kraken.out", codon_usage_dir="codon_usage/")
    for r in results:
        print(r.read_id, r.cai_vs_host)
"""
from plasmidScreen.api import (
    analyze_codon_adaptation,
    build_codon_reference,
    default_codon_usage_dir,
    run_screen,
    taxids_from_kraken_output,
    write_codon_adaptation_tsv,
)
from plasmidScreen.lib.codon_usage_db import CodonUsageStore
from plasmidScreen.lib.exceptions import CodonReferenceNotFoundError, MissingCodonReferenceError
from plasmidScreen.lib.models import (
    BuildCodonReferenceResult,
    CodonAdaptationResult,
    EngineeredScanResult,
    ReadEngineeringLabel,
    ScreenResult,
)

__all__ = [
    "analyze_codon_adaptation",
    "build_codon_reference",
    "default_codon_usage_dir",
    "run_screen",
    "taxids_from_kraken_output",
    "write_codon_adaptation_tsv",
    "BuildCodonReferenceResult",
    "CodonAdaptationResult",
    "CodonReferenceNotFoundError",
    "CodonUsageStore",
    "EngineeredScanResult",
    "MissingCodonReferenceError",
    "ReadEngineeringLabel",
    "ScreenResult",
]
