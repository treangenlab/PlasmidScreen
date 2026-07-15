"""Post-screen reference-hit enrichment (SeqScreen-style database provenance)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from plasmidScreen.lib.db_search import (
    DiamondReferenceSearchEngine,
    HitFilterConfig,
    ReferenceSearchConfig,
    enrich_screen_result_with_engine,
)
from plasmidScreen.lib.models import ScreenResult

logger = logging.getLogger(__name__)


def post_analysis_report(
        screen_result: ScreenResult,
        *,
        query_fasta: str | Path,
        reference_db: str | Path,
        threads: int = 4,
        top_n: int = 5,
        max_evalue: float = 1e-5,
        min_identity: float = 0.0,
        min_bitscore: float | None = None,
        database_source: str | None = None,
        diamond_output_path: str | Path | None = None,
        run_diamond: bool = True,
) -> ScreenResult:
    """
    Intersect ML / k-mer screen predictions with a reference database.

    Runs a batch DIAMOND search, selects best hits under the given thresholds,
    and returns an enriched :class:`ScreenResult` with per-read
    ``database_hits`` and run-level ``database_hits`` flattened for reporters.
    """
    config = ReferenceSearchConfig(
        diamond_db=Path(reference_db),
        database_source=database_source,
        threads=threads,
        filters=HitFilterConfig(
            top_n=top_n,
            max_evalue=max_evalue,
            min_identity=min_identity,
            min_bitscore=min_bitscore,
        ),
        output_path=Path(diamond_output_path) if diamond_output_path else None,
        run_diamond=run_diamond,
    )
    engine = DiamondReferenceSearchEngine(config)
    logger.info(
        "Enriching ScreenResult (%d reads) with reference hits from %s",
        len(screen_result.per_read),
        config.database_source or reference_db,
    )
    return enrich_screen_result_with_engine(screen_result, engine, query_fasta)
