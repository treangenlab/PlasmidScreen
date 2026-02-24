import typer
import logging

from pathlib import Path
from rich.logging import RichHandler
from typing import Annotated

from mofongo.lib.const import KrakenConfig, MetaMDBGConfig, MegahitConfig
from mofongo.lib.funcs import get_default_db_path
from mofongo.src.kraken.build_db import BuildDB
from mofongo.src.workflow.workflow import LRWorkflow, SRWorkflow

FORMAT = "%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
app = typer.Typer(pretty_exceptions_show_locals=False,
                  help="Welcome to KrakenEng. A tool to detect"
                       " engineered reads. We currently "
                       "support long reads and assemblies.")

DEFAULT_DB_PATH = get_default_db_path("KrakenEng")


@app.callback()
def main(ctx: typer.Context,
         verbose: Annotated[bool, typer.Option("--verbose", "-v",
                                               help="Enable verbose logging.")] = False):
    config = KrakenConfig()
    ctx.obj = config


@app.command()
def short_read(ctx: typer.Context, fasta_file: Annotated[str, typer.Argument(help="Fasta file ")],
               paired_end: Annotated[str, typer.Option(default="", help="Fasta file of paired end")],
               assembly: Annotated[bool, typer.Argument(help="Run on raw reads or assemble contigs")],
               output_report_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")],
               kraken_db_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")] = DEFAULT_DB_PATH):
    config: KrakenConfig = ctx.obj
    meta_config = MegahitConfig()
    SRWorkflow(fasta_file=fasta_file, paired_end=paired_end, output_report_path=output_report_path,
               kraken_db=kraken_db_path, kraken_config=config, assembly=assembly,
               assembly_config=meta_config).run_kraken()


@app.command()
def long_read(ctx: typer.Context, fasta_file: Annotated[str, typer.Argument(help="Fasta file ")],
              assembly: Annotated[bool, typer.Argument(help="Run on raw reads or assemble contigs")],
              output_report_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")],
              kraken_raw_output: Annotated[str, typer.Argument(help="Raw kraken output path")] ,
              kraken_db_path: Annotated[str, typer.Argument(help="Engineered KrakenDB")] = DEFAULT_DB_PATH,):
    kraken_config: KrakenConfig = ctx.obj
    meta_config = MetaMDBGConfig()
    LRWorkflow(fasta_file=fasta_file, output_report_path=output_report_path, kraken_db=kraken_db_path,
               kraken_raw_output=kraken_raw_output,
               kraken_config=kraken_config, assembly=assembly,
               assembly_config=meta_config).run_kraken()


@app.command(help="Build custom kraken database to build k-mer database. Recommend to add E. Coli")
def build_kraken_db(ctx: typer.Context,
                    natural_fastas: Annotated[str, typer.Argument(help="Fasta file for natural sequences")],
                    engineered_fastas: Annotated[str, typer.Argument(help="Fasta file for natural sequences")],
                    additional_fastas: Annotated[str, typer.Option(help="Fasta file of additional sequences for "
                                                                        "higher engineered sequence sensitivity. "
                                                                        "Recommended adding E. Coli genomes "
                                                                        "but will add runtime.")] = None,
                    db_name: Annotated[str, typer.Argument(help="Name of Custom kraken DB")] = "engineered_kraken_db"
                    ):
    config: KrakenConfig = ctx.obj
    build_db_obj = BuildDB(Path(DEFAULT_DB_PATH), db_name, natural_fastas, engineered_fastas,
                           additional_fastas,
                           config)
    build_db_obj.build_kraken_db()


if __name__ == "__main__":
    app()
