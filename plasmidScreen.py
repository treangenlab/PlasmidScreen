import typer
import logging

from rich.logging import RichHandler
from typing import Annotated

from plasmidScreen.lib.funcs import get_default_db_path
from plasmidScreen.lib.codon_usage_db import default_codon_usage_dir
from plasmidScreen.api import run_screen
from plasmidScreen.scripts.build_codon_usage_db import app as build_codon_db_app

FORMAT = "%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
app = typer.Typer(pretty_exceptions_show_locals=False,
                  help="Welcome to PlasmidScreen. A k-mer based tool to detect"
                       " engineered DNA.")

DEFAULT_DB_PATH = get_default_db_path("PlasmidScreen")
DEFAULT_CODON_USAGE_DIR = str(default_codon_usage_dir())


@app.command()
def screen(ctx: typer.Context, fasta_file: Annotated[str, typer.Argument(help="Fasta file ")],
           output_report_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")],
           kraken_raw_output: Annotated[str, typer.Argument(help="Raw kraken output path")],
           window_size: Annotated[int, typer.Option("--window_size",
                                                    help="Window size to scan for engineered k-mers")] = 200,
           engineered_k_mer_threshold: Annotated[int, typer.Option("--threshold", help="Threshold of required "
                                                                                       "identified engineered DNA.")] = 25,
           codon_usage_output: Annotated[str | None, typer.Option("--codon-usage-output",
                                                                  help="Optional codon usage TSV output path. "
                                                                       "Only produced for reads labeled Natural.")] = None,
           codon_usage_dir: Annotated[str, typer.Option("--codon-usage-dir",
                                                        help="Directory with codon_tables.json reference data.")] = DEFAULT_CODON_USAGE_DIR,
           kraken_db_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")] = DEFAULT_DB_PATH,
           threads: Annotated[int, typer.Option("--threads", help="Available threads to use.")] = 4,
           verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable verbose logging.")] = False
           ):
    run_screen(
        fasta_file,
        output_report_path,
        kraken_raw_output,
        kraken_db_path,
        threads=threads,
        window_size=window_size,
        engineered_kmer_threshold=engineered_k_mer_threshold,
        codon_usage_output_path=codon_usage_output,
        codon_usage_dir=codon_usage_dir,
    )


app.add_typer(build_codon_db_app, name="build-codon-db")


@app.callback()
def main():
    pass


if __name__ == "__main__":
    app()
