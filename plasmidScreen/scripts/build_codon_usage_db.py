#!/usr/bin/env python3
"""CLI for offline codon reference build (network required)."""
import logging
from pathlib import Path
from typing import Optional

import typer
from rich.logging import RichHandler

from plasmidScreen.api import build_codon_reference, default_codon_usage_dir, taxids_from_kraken_output

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
    timeout: int = typer.Option(30, "--timeout", help="HTTP timeout per Kazusa request (seconds)"),
):
    """Download codon usage tables and write JSON for airgapped screening."""
    data_dir = output_dir or default_codon_usage_dir()

    taxid_list = _parse_taxids(taxids, taxids_file)
    if kraken_output:
        taxid_list = sorted(set(taxid_list) | taxids_from_kraken_output(kraken_output))

    if not taxid_list:
        typer.echo("Provide --taxids, --taxids-file, or --kraken-output.", err=True)
        raise typer.Exit(code=1)

    typer.echo(f"Building codon reference at {data_dir} for {len(taxid_list)} taxid(s)...")
    result = build_codon_reference(
        data_dir,
        taxid_list,
        include_taxonomy=not skip_taxonomy,
        taxdump_dir=taxdump_dir,
        fetch_timeout=timeout,
    )
    typer.echo(
        f"Done. added={len(result.taxids_added)} skipped={len(result.taxids_skipped)} "
        f"failed={len(result.taxids_failed)} -> {result.data_dir / 'codon_tables.json'}"
    )


if __name__ == "__main__":
    app()
