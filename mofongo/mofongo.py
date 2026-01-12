import typer
import logging
from rich.logging import RichHandler
from typing import Optional, Annotated

from mofongo.lib.const import KrakenConfig
from mofongo.lib.funcs import get_default_db_path
from mofongo.src.kraken.build_db import BuildDB
from mofongo.src.workflow.workflow import LRWorkflow

FORMAT = "%(message)s"
logging.basicConfig(level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
app = typer.Typer(no_args_is_help=True, pretty_exceptions_show_locals=False,
                  help="Welcome to KrakenEng. A tool to detect"
                       " engineered reads. We currently "
                       "support long reads and assemblies.")

app.callback()
DEFAULT_DB_PATH = get_default_db_path("KrakenEng")


def main(ctx: typer.Context,
         verbose: Annotated[bool, typer.Option("--verbose",
                                               "-v", help="Enable verbose logging.")] = False, ):
    config = KrakenConfig()
    ctx.obj = config


@app.command()
def long_read(ctx: typer.Context,
              fasta_file: Annotated[str, typer.Argument(help="Fasta file ")],
              kraken_db_path: Annotated[str, typer.Option(default=DEFAULT_DB_PATH, help="Engineered KrakenDB")]
              ):
    config: KrakenConfig = ctx.obj
    LRWorkflow(fasta_file=fasta_file, kraken_db=kraken_db_path, kraken_config=config)


@app.command(help="Build custom kraken database to build k-mer database. Recommend to add E. Coli")
def build_kraken_db(
        natural_fastas: Annotated[str, typer.Argument(help="Fasta file for natural sequences")],
        engineered_fastas: Annotated[str, typer.Argument(help="Fasta file for natural sequences")],
        additional_fastas: Annotated[str, typer.Argument(help="Fasta file of additional sequences for higher "
                                                              "engineered sequence sensitivity. Recommended adding "
                                                              "E. Coli genomes, but will add runtime.")],
        kraken_db_name: Annotated[str, typer.Option("engineered_kraken_db", help="Name of Custom kraken DB")]
):
    build_db_obj = BuildDB(DEFAULT_DB_PATH.parent, kraken_db_name, natural_fastas, engineered_fastas, additional_fastas)
    build_db_obj.build_kraken_db()
