"""
PlasmidScreen: engineered DNA detection and codon adaptation analysis.

Library usage (airgapped screening requires a pre-built codon reference):

    from plasmidScreen import build_codon_reference, analyze_codon_adaptation, run_screen

    # Step 1 — networked machine: build reference JSON from Codon Statistics Database
    build_codon_reference("codon_usage/")  # all taxids in CSDB by default

    # 2. On airgapped machine — analyze (DIAMOND blastx + pre-built CSDB tables)
    results, _diamond_path = analyze_codon_adaptation(
        "reads.fa",
        diamond_db="/path/to/protein.dmnd",
        codon_usage_dir="codon_usage/",
    )
    for r in results:
        print(r.read_id, r.host_taxid, r.cai_vs_host)

    screen = run_screen("reads.fa", kraken_db, diamond_db="protein.dmnd",
                        codon_cai_engineered_threshold=0.7)
    for detail in screen.per_read:
        print(detail.read_id, detail.overall_label, detail.engineered_overall)
    print(screen.overall_synthetic_count)
"""
from plasmidScreen.api import (
    analyze_codon_adaptation,
    build_codon_database,
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
    ReadFlagDetail,
    ScreenResult,
    compute_engineered_overall,
)

__all__ = [
    "analyze_codon_adaptation",
    "build_codon_database",
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
    "ReadFlagDetail",
    "ScreenResult",
    "compute_engineered_overall",
]
