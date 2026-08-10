"""
PlasmidScreen public library API.

Screening and codon scoring are airgapped-safe once ``codon_tables.json`` has been
built offline with :func:`build_codon_reference` or :func:`build_codon_database`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference
from plasmidScreen.src.codon_usage.codon_usage_db import (
    default_codon_usage_dir,
)
from plasmidScreen.lib.models import (
    BuildCodonReferenceResult,
    ReferenceHit,
    ScreenResult,
    SupportDataTypes,
)
from plasmidScreen.lib.types import GeneSet
from plasmidScreen.src.plasmidScreen import Workflow
from plasmidScreen.src.post_hit_analysis import post_analysis_report

__all__ = [
    "build_codon_database",
    "run_screen",
    "ScreenResult",
    "BuildCodonReferenceResult",
    "ReferenceHit",
    "SupportDataTypes",
    "post_analysis_report",
]


def run_screen(
    fasta_file: str | Path,
    kraken_db: str | Path,
    engineered_report_path: str | Path | None = None,
    kraken_output_path: str | Path | None = None,
    threads: int = 4,
    window_size: int = 200,
    engineered_kmer_threshold: int = 10,
    codon_usage_dir: str | Path | None = None,
    run_kraken: bool = True,
    debug_write_kraken_output: bool = False,
    debug_write_kraken_report: bool = False,
    run_codon_usage: bool = True,
    codon_usage_output_path: str | Path | None = None,
    codon_cai_engineered_threshold: float = 0.6,
    diamond_db: str | Path | None = None,
    diamond_output_path: str | Path | None = None,
    debug_write_diamond_output: bool = False,
    run_diamond: bool = True,
    filter_hits: bool = False,
    reference_db: str | Path | None = None,
    reference_top_n: int = 5,
    reference_min_identity: float = 0.0,
    reference_min_mapq: int = 0,
    reference_min_bitscore: float | None = None,
    reference_database_source: str | None = None,
    reference_output_path: str | Path | None = None,
    run_reference_minimap: bool = True,
    reference_data_type: SupportDataTypes = SupportDataTypes.LONG_READ_ONT,
    quiet_mode:bool = True,
    visualize: bool = False
) -> ScreenResult:
    """
    Run engineered k-mer screening (Kraken2) and optional codon adaptation (DIAMOND + CSDB).

    Engineered detection uses Kraken2 minimizer blocks (taxid 32630) in a sliding window.
    Codon CAI runs on reads labeled **Natural** only, using DIAMOND blastx for ORF coordinates
    and host taxids, then pre-built codon usage tables for CAI.

    When ``filter_hits=True``, **engineered** reads are aligned with minimap2 against a
    nucleotide reference (FASTA or ``.mmi``) and top tiled hits are attached as
    ``ScreenResult.database_hits`` / per-read ``ReadFlagDetail.database_hits``.

    Parameters
    ----------
    fasta_file
        Input FASTA/FASTQ path.
    kraken_db
        Kraken2 database directory (required when ``run_kraken=True``).
    engineered_report_path
        If set, writes the engineered k-mer scan TSV to this path.
    kraken_output_path
        Save or load raw Kraken2 classifications (required when ``run_kraken=False`` or
        when ``debug_write_kraken_output=True``).
    threads
        Threads for Kraken / DIAMOND / minimap2 unless overridden elsewhere.
    run_kraken
        Run Kraken2 in-process; if False, ``kraken_output_path`` must point to existing output.
    debug_write_kraken_report
        Debug option to write kraken report
    debug_write_kraken_output
        Debug to write the output of kraken
    codon_usage_dir
        Directory with ``codon_tables.json``. Required when ``run_codon_usage=True``.
    run_codon_usage
        If False, skip DIAMOND and CAI (engineered scan only).
    codon_usage_output_path
        Optional path for codon adaptation TSV (Natural reads only).
    codon_cai_engineered_threshold
        If set, reads with CAI below this value get ``engineered_by_codon_cai=True`` and
        may be marked ``engineered_overall=True`` on :class:`~plasmidScreen.lib.models.ReadFlagDetail`
        even when the k-mer scan labeled them Natural.
    engineered_kmer_threshold
        Min engineered (32630) k-mers in a window to label a read Synthetic by k-mer scan.
    window_size
        Sliding window size (bp) for the k-mer scan.
    diamond_db
        DIAMOND protein database (``.dmnd``) for codon CAI only.
    diamond_output_path
        Save or load DIAMOND outfmt 6 TSV for codon CAI.
    debug_write_diamond_output
        Persist DIAMOND alignments to ``diamond_output_path``.
    run_diamond
        Run DIAMOND blastx for codon CAI; if False, load precomputed TSV.
    filter_hits
        After screening, minimap2-align engineered reads and attach best hits.
    reference_db
        Nucleotide FASTA or minimap2 ``.mmi`` index for engineered-read lookup.
    reference_top_n
        Maximum significant hits retained per engineered query.
    reference_min_identity
        Minimum percent identity (0–100) for retained hits.
    reference_min_mapq
        Minimum minimap2 mapping quality.
    reference_min_bitscore
        Optional minimum match-score threshold (for minimap2: matching bases).
    reference_database_source
        Human-readable database name/version stored on each :class:`ReferenceHit`.
    reference_output_path
        Save/load minimap2 PAF for the reference-search step.
    run_reference_minimap
        Run minimap2; if False, load precomputed PAF from ``reference_output_path``.
    reference_data_type
        minimap2 ``-x`` preset (:class:`~plasmidScreen.lib.models.SupportDataTypes`).

    Returns
    -------
    ScreenResult
        Includes ``per_read`` with ``engineered_overall`` / ``overall_label`` per read,
        optional ``database_hits`` from minimap2 reference search, and stored thresholds.
    """
    if run_codon_usage:
        if run_diamond and diamond_db is None:
            raise ValueError(
                "diamond_db is required when run_codon_usage=True and run_diamond=True."
            )
        if not run_diamond and diamond_output_path is None:
            raise ValueError(
                "diamond_output_path is required when run_codon_usage=True and run_diamond=False."
            )
        if debug_write_diamond_output and diamond_output_path is None:
            raise ValueError(
                "diamond_output_path is required when debug_write_diamond_output=True."
            )
    if filter_hits:
        if run_reference_minimap and reference_db is None:
            raise ValueError(
                "reference_db (FASTA or .mmi) is required when "
                "filter_hits=True and run_reference_minimap=True."
            )
        if not run_reference_minimap and reference_output_path is None:
            raise ValueError(
                "reference_output_path is required when run_reference_minimap=False."
            )

    workflow = Workflow(
        str(fasta_file),
        str(engineered_report_path) if engineered_report_path else None,
        str(kraken_db),
        threads,
        str(kraken_output_path) if kraken_output_path else None,
        window_size=window_size,
        engineered_kmer_threshold=engineered_kmer_threshold,
        codon_usage_output_path=str(codon_usage_output_path) if codon_usage_output_path else None,
        codon_usage_dir=str(codon_usage_dir) if codon_usage_dir else None,
        run_kraken=run_kraken,
        debug_write_kraken_output=debug_write_kraken_output,
        debug_write_kraken_report=debug_write_kraken_report,
        run_codon_usage=run_codon_usage,
        codon_cai_engineered_threshold=codon_cai_engineered_threshold,
        diamond_db=str(diamond_db) if diamond_db else None,
        diamond_output_path=str(diamond_output_path) if diamond_output_path else None,
        debug_write_diamond_output=debug_write_diamond_output,
        run_diamond=run_diamond,
        run_reference_search=filter_hits,
        reference_db=str(reference_db) if reference_db else None,
        reference_top_n=reference_top_n,
        reference_min_identity=reference_min_identity,
        reference_min_mapq=reference_min_mapq,
        reference_min_bitscore=reference_min_bitscore,
        reference_database_source=reference_database_source,
        reference_output_path=(
            str(reference_output_path) if reference_output_path else None
        ),
        run_reference_minimap=run_reference_minimap,
        reference_data_type=reference_data_type,
    )
    screen_results = workflow.run()
    if visualize:
        run_visualization_routine(screen_results)
    return screen_results


def build_codon_database(
    *,
    output_dir: str | Path | None = None,
    taxids: Iterable[str | int] | None = None,
    taxids_file: str | Path | None = None,
    csdb_archive: str | Path | None = None,
    download_csdb: bool = True,
    gene_set: GeneSet = "nuclear",
) -> BuildCodonReferenceResult:
    """
    Serves as an API wrapper for build_codon_reference so users can easily create codon database as needed.

    Build the codon usage reference for CAI scoring. It attempts to use a
    taxids file if provided to grab the codon references. If none is provided, it imports
    **every** taxid in the CSDB archive for ``gene_set``. Downloads CSDB archive if download_csdb is provided as true.
    Writes ``codon_tables.json`` and optionally ``taxonomy_parents.json``
    under ``output_dir`` (default: PlasmidScreen user data ``codon_usage/``).
    """
    data_dir = Path(output_dir) if output_dir else default_codon_usage_dir()

    resolved: set[str] = set()
    if taxids:
        resolved.update(str(t) for t in taxids if str(t) not in ("0", ""))

    if taxids_file:
        path = Path(taxids_file)
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip().split("#")[0].strip()
            if line:
                resolved.add(line)

    taxid_list = sorted(resolved) if resolved else None

    return build_codon_reference(
        data_dir,
        taxid_list,
        csdb_archive=csdb_archive,
        download_csdb=download_csdb,
        gene_set=gene_set,
    )
