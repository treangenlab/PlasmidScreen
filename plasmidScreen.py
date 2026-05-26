import typer
import logging

from rich.logging import RichHandler
from typing import Annotated

from plasmidScreen.lib.funcs import get_default_db_path
from plasmidScreen.src.plasmidScreen import Workflow

FORMAT = "%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
app = typer.Typer(pretty_exceptions_show_locals=False,
                  help="Welcome to PlasmidScreen. A k-mer based tool to detect"
                       " engineered DNA.")

DEFAULT_DB_PATH = get_default_db_path("PlasmidScreen")


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
                                                                       "Only produced when no engineered blocks are detected.")] = None,
           kraken_db_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")] = DEFAULT_DB_PATH,
           threads: Annotated[int, typer.Option("--threads", help="Available threads to use.")] = 4,
           verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable verbose logging.")] = False
           ):
    Workflow(fasta_file, output_report_path, kraken_db_path, threads, kraken_raw_output, window_size,
             engineered_k_mer_threshold, codon_usage_output_path=codon_usage_output).run()


@app.callback()
def main():
    pass


if __name__ == "__main__":
    app()
