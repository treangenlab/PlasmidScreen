"""Post-screen reference-hit enrichment via blastn (engineered reads only)."""
from __future__ import annotations

import logging
from pathlib import Path

from plasmidScreen.lib.db_search import (
    BlastnReferenceSearchEngine,
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
        min_bitscore: float = 50.0,
        max_evalue: float | None = None,
        database_source: str | None = None,
        output_path: str | Path | None = None,
        run_blast: bool = True,
        blastn_binary: str = "blastn",
        blast_task: str = "blastn",
) -> ScreenResult:
    """
    BLAST engineered reads against a nucleotide reference and attach top hits.

    Keeps alignments with bitscore ≥ ``min_bitscore``, then retains the sample's
    top-``top_n`` reference sequences (ranked by # reads covered, bitscore sum
    as tiebreak). Mapping is informational only: overall engineered / natural
    labels are not rewritten. Unengineered reads are not annotated.
    """
    config = ReferenceSearchConfig(
        reference_db=Path(reference_db),
        database_source=database_source,
        threads=threads,
        filters=HitFilterConfig(
            min_bitscore=min_bitscore,
            top_n_references=top_n,
            max_evalue=max_evalue,
        ),
        output_path=Path(output_path) if output_path else None,
        run_blast=run_blast,
        blastn_binary=blastn_binary,
        blast_task=blast_task,
    )
    engine = BlastnReferenceSearchEngine(config)
    logger.info(
        "Enriching ScreenResult (%d engineered reads) with blastn hits from %s",
        len(screen_result.engineered_read_ids),
        config.database_source or reference_db,
    )
    enriched = enrich_screen_result_with_engine(
        screen_result, engine, Path(query_fasta)
    )
    return enriched if enriched is not None else screen_result
