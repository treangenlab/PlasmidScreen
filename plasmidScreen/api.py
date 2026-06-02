"""
PlasmidScreen public library API.

All screening assumes codon reference data was built offline (airgapped-safe).
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

from plasmidScreen.lib.codon_usage_build import build_codon_reference, default_reference_taxids
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
    "taxids_from_kraken_output",
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
    Run engineered k-mer screening and optional codon adaptation analysis.

    Parameters
    ----------
    fasta_file:str
        path of fasta file for reading and engineering detection
    kraken_db: str

    engineered_report_path
        When set, writes the engineered k-mer scan TSV to this path.
    kraken_output_path
    threads
    window_size
    engineered_kmer_threshold
    codon_usage_dir
        Pre-built reference directory (codon_tables.json). Required when
        run_codon_usage is True.
    run_kraken
        If False, kraken_output_path must be provided (airgapped Kraken step done).
    debug_write_diamond_output
    debug_write_kraken_output
    debug_write_kraken_report
    run_codon_usage
    codon_usage_output_path
    codon_cai_engineered_threshold
    diamond_db
    diamond_threads
    diamond_extra_args
    run_diamond

    diamond_output_path
        Path to save (``debug_write_diamond_output``) or load (``run_diamond=False``)
        DIAMOND outfmt 6 TSV.
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
    Convenience wrapper to build the codon usage reference database. Requires internet connection
    to fetch resources from Codon Statistics Database Krishnamurthy et al. 2022, Molecular Biology and Evolution.

    This mirrors the build CLI behavior:
    - union of explicit taxids + taxids_file + kraken_output taxids
    - when nothing is provided, imports every taxid in the CSDB archive
    - writes codon_tables.json (+ optional taxonomy_parents.json) under output_dir
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
