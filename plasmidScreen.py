import typer
import logging

from rich.logging import RichHandler
from typing import Annotated

from plasmidScreen.lib.funcs import get_default_db_path
from plasmidScreen.src.codon_usage.codon_usage_db import default_codon_usage_dir
from plasmidScreen.src.codon_usage.codon_usage_sources import default_csdb_archive_path
from plasmidScreen.api import run_screen
from pathlib import Path
from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference
from plasmidScreen.lib.types import GeneSet


FORMAT = "%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
app = typer.Typer(
    pretty_exceptions_show_locals=False,
    help=(
        "PlasmidScreen: Kraken2 engineered k-mer screening and optional "
        "DIAMOND + codon adaptation (CAI) analysis."
    ),
)

DEFAULT_DB_PATH = get_default_db_path("PlasmidScreen")
DEFAULT_CODON_USAGE_DIR = str(default_codon_usage_dir())


@app.command()
def screen(ctx: typer.Context, fasta_file: Annotated[str, typer.Argument(help="Fasta file ")],
           output_report_path: Annotated[str, typer.Argument(help="Engineered k-mer scan report TSV path")],
           kraken_output_path: Annotated[str | None, typer.Option("--kraken-output-path",
                                                                  help="Optional path to read/write raw Kraken2 "
                                                                       "classifications. Required when "
                                                                       "--run-kraken=false or when "
                                                                       "--debug-write-kraken-out is set.")] = None,
           window_size: Annotated[int, typer.Option("--window_size",
                                                    help="Window size to scan for engineered k-mers")] = 200,
           engineered_k_mer_threshold: Annotated[int, typer.Option("--threshold",
                                                                   help="Threshold of required "
                                                                        "identified engineered DNA.")] = 10,
           codon_usage_output: Annotated[str | None, typer.Option("--codon-usage-output",
                                                                  help="Optional codon usage TSV output path. Only "
                                                                       "produced for reads labeled Natural.")] = None,
           codon_usage_dir: Annotated[str, typer.Option("--codon-usage-dir",
                                                        help="Directory with codon_tables.json "
                                                             "reference data.")] = DEFAULT_CODON_USAGE_DIR,
           codon_cai_engineered_threshold: Annotated[float | None, typer.Option("--codon-cai-threshold",
                                                                                help="If set, Natural reads with "
                                                                                     "CAI below this value are "
                                                                                     "engineered_overall in "
                                                                                     "ScreenResult (library).")] = None,
           diamond_db: Annotated[str | None, typer.Option("--diamond-db",
                                                          help="DIAMOND database (.dmnd) used to infer "
                                                               "ORFs/host taxid for codon CAI. Required when "
                                                               "codon usage is enabled.")] = None,
           diamond_threads: Annotated[int, typer.Option("--diamond-threads",
                                                        help="Threads for DIAMOND (defaults to "
                                                             "--threads if not set).")] = 0,
           diamond_extra_args: Annotated[str | None, typer.Option("--diamond-args",
                                                                  help="Extra DIAMOND args passed verbatim "
                                                                       "(comma-separated).")] = None,
           diamond_output_path: Annotated[str | None, typer.Option("--diamond-output-path",
                                                                   help="Save/load DIAMOND outfmt 6 TSV. "
                                                                        "Required with --debug-write-diamond-out "
                                                                        "or --no-run-diamond.")] = None,
           debug_write_diamond_out: Annotated[bool, typer.Option("--debug-write-diamond-out",
                                                                 help="Write DIAMOND TSV to --diamond-output-path "
                                                                      "(debug/reuse).")] = False,
           run_diamond: Annotated[bool, typer.Option("--run-diamond/--no-run-diamond",
                                                     help="Run DIAMOND blastx (default). If disabled, "
                                                          "--diamond-output-path must exist.")] = True,
           skip_codon_usage: Annotated[bool, typer.Option("--skip-codon-usage",
                                                          help="Skip codon adaptation/CAI analysis "
                                                               "(engineered k-mer scan only).")] = False,
           debug_write_kraken_out: Annotated[bool, typer.Option("--debug-write-kraken-out",
                                                                help="Write raw Kraken2 classifications to "
                                                                     "--kraken-output-path (debug).")] = False,
           debug_write_kraken_report: Annotated[bool, typer.Option("--debug-write-kraken-report",
                                                                   help="Write Kraken2 --report file (debug).")] = False,
           run_kraken: Annotated[bool, typer.Option("--run-kraken/--no-run-kraken",
                                                    help="Run Kraken2 as part of the pipeline (default). If disabled, "
                                                         "--kraken-output-path must point to an "
                                                         "existing classifications file.")] = True,
           kraken_db_path: Annotated[str, typer.Argument(help="Kraken2 database path")] = DEFAULT_DB_PATH,
           threads: Annotated[int, typer.Option("--threads", help="Available threads to use.")] = 4,
           ) -> None:
    diamond_thread_count = threads if diamond_threads <= 0 else diamond_threads
    diamond_args_list = None
    if diamond_extra_args:
        diamond_args_list = [a.strip() for a in diamond_extra_args.split(",") if a.strip()]
    if not run_kraken and not kraken_output_path:
        raise typer.BadParameter("--kraken-output-path is required when --no-run-kraken is set")
    if debug_write_kraken_out and not kraken_output_path:
        raise typer.BadParameter("--kraken-output-path is required when --debug-write-kraken-out is set")
    if not skip_codon_usage:
        if run_diamond and not diamond_db:
            raise typer.BadParameter(
                "--diamond-db is required when codon usage is enabled (omit --skip-codon-usage to disable)."
            )
        if not run_diamond and not diamond_output_path:
            raise typer.BadParameter(
                "--diamond-output-path is required when --no-run-diamond is set"
            )
        if debug_write_diamond_out and not diamond_output_path:
            raise typer.BadParameter(
                "--diamond-output-path is required when --debug-write-diamond-out is set"
            )
    run_screen(
        fasta_file,
        kraken_db_path,
        engineered_report_path=output_report_path,
        kraken_output_path=kraken_output_path,
        threads=threads,
        window_size=window_size,
        engineered_kmer_threshold=engineered_k_mer_threshold,
        codon_usage_output_path=codon_usage_output,
        codon_usage_dir=codon_usage_dir,
        codon_cai_engineered_threshold=codon_cai_engineered_threshold,
        run_codon_usage=not skip_codon_usage,
        debug_write_kraken_output=debug_write_kraken_out,
        debug_write_kraken_report=debug_write_kraken_report,
        run_kraken=run_kraken,
        diamond_db=diamond_db,
        diamond_threads=diamond_thread_count,
        diamond_extra_args=diamond_args_list,
        diamond_output_path=diamond_output_path,
        debug_write_diamond_output=debug_write_diamond_out,
        run_diamond=run_diamond,
    )


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
        csdb_archive=archive,
        download_csdb=not no_download_csdb,
        gene_set=gene_set,
    )
    typer.echo(
        f"Done. added={len(result.taxids_added)} skipped={len(result.taxids_skipped)} "
        f"failed={len(result.taxids_failed)} -> {result.data_dir / 'codon_tables.json'}"
    )


@app.callback()
def main() -> None:
    pass


if __name__ == "__main__":
    app()
