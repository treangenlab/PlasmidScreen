"""
PlasmidScreen public library API.

Screening and codon scoring are airgapped-safe once ``codon_tables.json`` has been
built offline with :func:`build_codon_reference` or :func:`build_codon_database`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

from plasmidScreen.lib.codon_usage_build import build_codon_reference
from plasmidScreen.lib.codon_usage_sources import all_csdb_taxids
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
from plasmidScreen.lib.types import GeneSet
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
    write_codon_adaptation_tsv,
)
from plasmidScreen.src.plasmidScreen import Workflow

__all__ = [
    "analyze_codon_adaptation",
    "build_codon_reference",
    "build_codon_database",
    "all_csdb_taxids",
    "default_codon_usage_dir",
    "run_screen",
    "write_codon_adaptation_tsv",
    "BuildCodonReferenceResult",
    "CodonAdaptationResult",
    "CodonUsageStore",
    "ScreenResult",
]


def run_screen(
    fasta_file: str | Path,
    kraken_db: str | Path,
    *,
    engineered_report_path: str | Path | None = None,
    kraken_output_path: str | Path | None = None,
    threads: int = 4,
    window_size: int = 200,
    engineered_kmer_threshold: int = 25,
    codon_usage_dir: str | Path | None = None,
    run_kraken: bool = True,
    debug_write_kraken_output: bool = False,
    debug_write_kraken_report: bool = False,
    run_codon_usage: bool = True,
    codon_usage_output_path: str | Path | None = None,
    codon_cai_engineered_threshold: float | None = None,
    diamond_db: str | Path | None = None,
    diamond_threads: int = 4,
    diamond_extra_args: list[str] | None = None,
    diamond_output_path: str | Path | None = None,
    debug_write_diamond_output: bool = False,
    run_diamond: bool = True,
) -> ScreenResult:
    """
    Run engineered k-mer screening (Kraken2) and optional codon adaptation (DIAMOND + CSDB).

    Engineered detection uses Kraken2 minimizer blocks (taxid 32630) in a sliding window.
    Codon CAI runs on reads labeled **Natural** only, using DIAMOND blastx for ORF coordinates
    and host taxids, then pre-built codon usage tables for CAI.

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
    run_kraken
        Run Kraken2 in-process; if False, ``kraken_output_path`` must point to existing output.
    codon_usage_dir
        Directory with ``codon_tables.json``. Required when ``run_codon_usage=True``.
    run_codon_usage
        If False, skip DIAMOND and CAI (engineered scan only).
    codon_usage_output_path
        Optional path for codon adaptation TSV (Natural reads only).
    codon_cai_engineered_threshold
        If set, reads with CAI below this value get ``engineered_by_codon_cai=True`` and
        may be marked ``engineered_overall=True`` on :class:`~plasmidScreen.lib.models.ReadFlagDetail`
        even when the k-mer scan labeled them Natural. Also adds columns to the codon TSV when written.
    engineered_kmer_threshold
        Min engineered (32630) k-mers in a window to label a read Synthetic by k-mer scan.
    window_size
        Sliding window size (bp) for the k-mer scan.
    diamond_db
        DIAMOND protein database (``.dmnd``). Required when ``run_codon_usage=True`` and
        ``run_diamond=True``.
    diamond_output_path
        Save or load DIAMOND outfmt 6 TSV (for ``debug_write_diamond_output`` or
        ``run_diamond=False``).
    debug_write_diamond_output
        Persist DIAMOND alignments to ``diamond_output_path``.
    run_diamond
        Run DIAMOND blastx; if False, load precomputed TSV from ``diamond_output_path``.

    Returns
    -------
    ScreenResult
        Includes ``per_read`` with ``engineered_overall`` / ``overall_label`` per read,
        run-level ``overall_synthetic_count`` / ``engineered_read_ids``, and stored
        threshold values used for the combined decision.
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
    workflow = Workflow(
        str(fasta_file),
        str(engineered_report_path) if engineered_report_path else None,
        str(kraken_db),
        threads,
        str(kraken_output_path) if kraken_output_path else None,
        write_engineered_report=engineered_report_path is not None,
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
        diamond_threads=diamond_threads,
        diamond_extra_args=diamond_extra_args,
        diamond_output_path=str(diamond_output_path) if diamond_output_path else None,
        debug_write_diamond_output=debug_write_diamond_output,
        run_diamond=run_diamond,
    )
    return workflow.run()


def build_codon_database(
    *,
    output_dir: str | Path | None = None,
    taxids: Iterable[str | int] | None = None,
    taxids_file: str | Path | None = None,
    kraken_output: str | Path | None = None,
    include_taxonomy: bool = True,
    taxdump_dir: str | Path | None = None,
    csdb_archive: str | Path | None = None,
    download_csdb: bool = True,
    gene_set: GeneSet = "nuclear",
) -> BuildCodonReferenceResult:
    """
    Build the codon usage reference for airgapped CAI scoring (network required).

    Mirrors ``python plasmidScreen.py build``: unions explicit taxids from ``taxids``,
    ``taxids_file``, and classified taxids from ``kraken_output``. When no taxids are
    given, imports **every** taxid in the CSDB archive for ``gene_set``.

    Writes ``codon_tables.json`` and optionally ``taxonomy_parents.json`` under
    ``output_dir`` (default: PlasmidScreen user data ``codon_usage/``).
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

    if kraken_output:
        resolved.update(taxids_from_kraken_output(Path(kraken_output)))

    taxid_list = sorted(resolved) if resolved else None

    return build_codon_reference(
        data_dir,
        taxid_list,
        include_taxonomy=include_taxonomy,
        taxdump_dir=taxdump_dir,
        csdb_archive=csdb_archive,
        download_csdb=download_csdb,
        gene_set=gene_set,
    )
