#!/usr/bin/env python3
"""CLI for offline codon reference build (network required for CSDB/taxdump download)."""
import logging
from pathlib import Path
from typing import Optional

import typer
from rich.logging import RichHandler

from plasmidScreen.api import build_codon_reference, default_codon_usage_dir, taxids_from_kraken_output
from plasmidScreen.lib.codon_usage_build import default_reference_taxids
from plasmidScreen.lib.codon_usage_sources import default_csdb_archive_path

FORMAT = "%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])

app = typer.Typer(help="Build codon usage reference tables (offline/airgap prep step).")


def _parse_taxids(taxids: Optional[str], taxids_file: Optional[Path]) -> list[str]:
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
    output_dir: Path = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory for codon_tables.json and taxonomy_parents.json",
    ),
    taxids: Optional[str] = typer.Option(None, "--taxids", help="Comma-separated NCBI taxonomy IDs"),
    taxids_file: Optional[Path] = typer.Option(None, "--taxids-file", help="One taxid per line"),
    kraken_output: Optional[Path] = typer.Option(
        None, "--kraken-output", help="Kraken2 output; all classified taxids are included"
    ),
    skip_taxonomy: bool = typer.Option(
        False, "--skip-taxonomy", help="Skip NCBI taxdump (no lineage resolution)"
    ),
    taxdump_dir: Optional[Path] = typer.Option(None, "--taxdump-dir", help="NCBI taxdump cache directory"),
    csdb_archive: Optional[Path] = typer.Option(
        None,
        "--csdb-archive",
        help="Path to codonstatsdb_March2022.tar.gz (default: PlasmidScreen data dir)",
    ),
    no_download_csdb: bool = typer.Option(
        False,
        "--no-download-csdb",
        help="Do not download CSDB; require --csdb-archive to exist",
    ),
    gene_set: str = typer.Option(
        "nuclear",
        "--gene-set",
        help="CSDB gene set: nuclear, ribosomal, mitochondrial, or plastid",
    ),
):
    """Build codon usage tables from the Codon Statistics Database for airgapped CAI scoring."""
    data_dir = output_dir or default_codon_usage_dir()
    archive = csdb_archive or default_csdb_archive_path()

    taxid_list = _parse_taxids(taxids, taxids_file)
    if kraken_output:
        taxid_list = sorted(set(taxid_list) | taxids_from_kraken_output(kraken_output))

    if not taxid_list:
        taxid_list = default_reference_taxids()
        typer.echo(
            f"No taxids specified; using {len(taxid_list)} default reference taxid(s). "
            "Override with --taxids, --taxids-file, or --kraken-output."
        )

    typer.echo(
        f"Building codon reference at {data_dir} for {len(taxid_list)} taxid(s) "
        f"from CSDB archive {archive} ..."
    )
    result = build_codon_reference(
        data_dir,
        taxid_list,
        use_default_taxids=False,
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
