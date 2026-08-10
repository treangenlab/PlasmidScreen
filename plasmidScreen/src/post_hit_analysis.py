"""Post-screen reference-hit enrichment via minimap2 (engineered reads only)."""
from __future__ import annotations

import logging
from pathlib import Path

from plasmidScreen.lib.db_search import (
    HitFilterConfig,
    MinimapReferenceSearchEngine,
    ReferenceSearchConfig,
    enrich_screen_result_with_engine,
)
from plasmidScreen.lib.models import ScreenResult, SupportDataTypes

logger = logging.getLogger(__name__)


def post_analysis_report(
        screen_result: ScreenResult,
        *,
        query_fasta: str | Path,
        reference_db: str | Path,
        threads: int = 4,
        top_n: int = 5,
        min_identity: float = 0.0,
        min_mapq: int = 0,
        min_bitscore: float | None = None,
        database_source: str | None = None,
        output_path: str | Path | None = None,
        run_minimap: bool = True,
        data_type: SupportDataTypes = SupportDataTypes.LONG_READ_ONT,
) -> ScreenResult:
    """
    Align engineered reads to a nucleotide reference and attach tiled best hits.

    Uses minimap2 (PAF). Unengineered reads keep their ML / k-mer labels and are
    not annotated with database-match notes.
    """
    config = ReferenceSearchConfig(
        reference_db=Path(reference_db),
        database_source=database_source,
        threads=threads,
        data_type=data_type,
        filters=HitFilterConfig(
            top_n=top_n,
            min_identity=min_identity,
            min_mapq=min_mapq,
            min_bitscore=min_bitscore,
        ),
        output_path=Path(output_path) if output_path else None,
        run_minimap=run_minimap,
    )
    engine = MinimapReferenceSearchEngine(config)
    logger.info(
        "Enriching ScreenResult (%d engineered reads) with minimap2 hits from %s",
        len(screen_result.engineered_read_ids),
        config.database_source or reference_db,
    )
    return enrich_screen_result_with_engine(screen_result, engine, Path(query_fasta))
