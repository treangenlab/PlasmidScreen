"""
plasmidScreen: engineered DNA detection and codon adaptation analysis.
"""
from plasmidScreen.api import (
    PrefetchResult,
    build_codon_database,
    default_genbank_cache_dir,
    map_and_visualize_engineered,
    map_engineered_reads,
    prefetch_genbank_from_fasta,
    run_screen,
    screen_result_from_engineered_report,
    visualize_engineered_mappings,
)
from plasmidScreen.src.analyze_codon_usage import analyze_codon_adaptation
from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference
from plasmidScreen.src.codon_usage.codon_usage_db import CodonUsageStore
from plasmidScreen.lib.exceptions import CodonReferenceNotFoundError, MissingCodonReferenceError
from plasmidScreen.lib.models import (
    BuildCodonReferenceResult,
    CodonAdaptationResult,
    EngineeredMappingHit,
    EngineeredScanResult,
    MappedQueryInterval,
    ReadEngineeringLabel,
    ReadFlagDetail,
    ReferenceFeatureInterval,
    SampleReferenceRank,
    ScreenResult,
    VisualizationReport,
    compute_engineered_overall,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "build_codon_database",
    "run_screen",
    "map_engineered_reads",
    "visualize_engineered_mappings",
    "map_and_visualize_engineered",
    "screen_result_from_engineered_report",
    "prefetch_genbank_from_fasta",
    "default_genbank_cache_dir",
    "PrefetchResult",
    "BuildCodonReferenceResult",
    "CodonAdaptationResult",
    "CodonReferenceNotFoundError",
    "CodonUsageStore",
    "EngineeredMappingHit",
    "EngineeredScanResult",
    "MappedQueryInterval",
    "MissingCodonReferenceError",
    "ReadEngineeringLabel",
    "ReadFlagDetail",
    "ReferenceFeatureInterval",
    "SampleReferenceRank",
    "ScreenResult",
    "VisualizationReport",
    "compute_engineered_overall",
    "build_codon_reference",
    "analyze_codon_adaptation",
]
