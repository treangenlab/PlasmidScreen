#!/usr/bin/env python3
"""CLI for offline codon reference build (network required for CSDB/taxdump download)."""
from __future__ import annotations

import logging
from pathlib import Path

import typer
from rich.logging import RichHandler

from plasmidScreen.api import build_codon_reference, default_codon_usage_dir, taxids_from_kraken_output
from plasmidScreen.lib.codon_usage_build import build_codon_reference
from plasmidScreen.lib.codon_usage_sources import default_csdb_archive_path
from plasmidScreen.lib.types import GeneSet

FORMAT = "%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])

app = typer.Typer(help="Build codon usage reference tables (offline/airgap prep step).")


def _parse_taxids(taxids: str | None, taxids_file: Path | None) -> list[str]:
    found: set[str] = set()
    if taxids:
        found.update(t.strip() for t in taxids.split(",") if t.strip())
    if taxids_file:
        for line in taxids_file.read_text().splitlines():
            line = line.strip().split("#")[0].strip()
            if line:
                found.add(line)
    return sorted(found)


@app.command()
def build(
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory for codon_tables.json and taxonomy_parents.json",
    ),
    taxids: str | None = typer.Option(None, "--taxids", help="Comma-separated NCBI taxonomy IDs"),
    taxids_file: Path | None = typer.Option(None, "--taxids-file", help="One taxid per line"),
    kraken_output: Path | None = typer.Option(
        None, "--kraken-output", help="Kraken2 output; all classified taxids are included"
    ),
    skip_taxonomy: bool = typer.Option(
        False, "--skip-taxonomy", help="Skip NCBI taxdump (no lineage resolution)"
    ),
    taxdump_dir: Path | None = typer.Option(None, "--taxdump-dir", help="NCBI taxdump cache directory"),
    csdb_archive: Path | None = typer.Option(
        None,
        "--csdb-archive",
        help="Path to codonstatsdb_March2022.tar.gz (default: PlasmidScreen data dir)",
    ),
    no_download_csdb: bool = typer.Option(
        False,
        "--no-download-csdb",
        help="Do not download CSDB; require --csdb-archive to exist",
    ),
    gene_set: GeneSet = typer.Option(
        "nuclear",
        "--gene-set",
        help="CSDB gene set: nuclear, ribosomal, mitochondrial, or plastid",
    ),
) -> None:
    """Build codon usage tables from the Codon Statistics Database for airgapped CAI scoring."""
    data_dir = output_dir or default_codon_usage_dir()
    archive = csdb_archive or default_csdb_archive_path()

    taxid_list = _parse_taxids(taxids, taxids_file)
    if kraken_output:
        taxid_list = sorted(set(taxid_list) | taxids_from_kraken_output(kraken_output))

    if taxid_list:
        typer.echo(
            f"Building codon reference at {data_dir} for {len(taxid_list)} taxid(s) "
            f"from CSDB archive {archive} ..."
        )
    else:
        typer.echo(
            f"Building codon reference at {data_dir} for all taxids in CSDB archive {archive} ..."
        )

    result = build_codon_reference(
        data_dir,
        taxid_list or None,
        include_taxonomy=not skip_taxonomy,
        taxdump_dir=taxdump_dir,
        csdb_archive=archive,
        download_csdb=not no_download_csdb,
        gene_set=gene_set,
    )
    typer.echo(
        f"Done. added={len(result.taxids_added)} skipped={len(result.taxids_skipped)} "
        f"failed={len(result.taxids_failed)} -> {result.data_dir / 'codon_tables.json'}"
    )


if __name__ == "__main__":
    app()
