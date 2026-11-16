"""PlasmidScreen command-line interface."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

from plasmidScreen.api import (
    default_genbank_cache_dir,
    map_and_visualize_engineered,
    map_engineered_reads,
    prefetch_genbank_from_fasta,
    run_screen,
    screen_result_from_engineered_report,
)
from plasmidScreen.lib.exceptions import PlasmidScreenError
from plasmidScreen.lib.funcs import get_default_db_path
from plasmidScreen.lib.logging_util import configure_logging
from plasmidScreen.lib.models import BACKGROUND_RATE, PVALUE_ALPHA
from plasmidScreen.lib.types import GeneSet
from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference
from plasmidScreen.src.codon_usage.codon_usage_db import default_codon_usage_dir
from plasmidScreen.src.codon_usage.codon_usage_sources import default_csdb_archive_path

app = typer.Typer(
    pretty_exceptions_show_locals=False,
    help=(
        "PlasmidScreen: Kraken2 engineered k-mer screening, optional "
        "DIAMOND + codon adaptation (CAI), and optional blastn map/visualize."
    ),
)

DEFAULT_DB_PATH = get_default_db_path("PlasmidScreen")
DEFAULT_CODON_USAGE_DIR = str(default_codon_usage_dir())
logger = logging.getLogger(__name__)


def _exit_on_error(exc: BaseException) -> None:
    """Log a user-facing error and exit with code 1."""
    logger.error("%s", exc)
    raise typer.Exit(code=1) from exc


@app.command()
def screen(
        ctx: typer.Context,
        fasta_file: Annotated[
            str,
            typer.Argument(help="Input FASTA/FASTQ (plain or .gz)"),
        ],
        output_report_path: Annotated[
            str, typer.Argument(help="Engineered k-mer scan report TSV path")
        ],
        kraken_output_path: Annotated[
            str | None,
            typer.Option(
                "--kraken-output-path",
                help=(
                    "Optional path to read/write raw Kraken2 classifications. "
                    "Required when --run-kraken=false or when "
                    "--debug-write-kraken-out is set."
                ),
            ),
        ] = None,
        window_size: Annotated[
            int,
            typer.Option(
                "--window_size", help="Window size to scan for engineered k-mers"
            ),
        ] = 200,
        engineered_k_mer_threshold: Annotated[
            int,
            typer.Option(
                "--threshold",
                help="Threshold of required identified engineered DNA.",
            ),
        ] = 10,
        kmer_size: Annotated[
            int,
            typer.Option(
                "--kmer-size",
                help=(
                    "k-mer length for the engineered window scan (must match "
                    "the Kraken2 database; default 21)."
                ),
            ),
        ] = 21,
        minimizer_size: Annotated[
            int,
            typer.Option(
                "--minimizer-size",
                help=(
                    "Minimizer length for the engineered window scan (must "
                    "match the Kraken2 database; default 17)."
                ),
            ),
        ] = 17,
        codon_usage_output: Annotated[
            str | None,
            typer.Option(
                "--codon-usage-output",
                help=(
                    "Optional codon usage TSV output path. Only produced for "
                    "reads labeled Natural."
                ),
            ),
        ] = None,
        codon_usage_dir: Annotated[
            str,
            typer.Option(
                "--codon-usage-dir",
                help="Directory with codon_tables.json reference data.",
            ),
        ] = DEFAULT_CODON_USAGE_DIR,
        codon_cai_engineered_threshold: Annotated[
            float | None,
            typer.Option(
                "--codon-cai-threshold",
                help=(
                    "If set, Natural reads with CAI above this value are "
                    "engineered_overall in ScreenResult (library)."
                ),
            ),
        ] = None,
        diamond_db: Annotated[
            str | None,
            typer.Option(
                "--diamond-db",
                help=(
                    "DIAMOND database (.dmnd) used to infer ORFs/host taxid "
                    "for codon CAI. Required when codon usage is enabled."
                ),
            ),
        ] = None,
        diamond_output_path: Annotated[
            str | None,
            typer.Option(
                "--diamond-output-path",
                help=(
                    "Save/load DIAMOND outfmt 6 TSV. Required with "
                    "--debug-write-diamond-out or --no-run-diamond."
                ),
            ),
        ] = None,
        seq_type: Annotated[
            int,
            typer.Option(
                "--seq-type",
                help=(
                    "Specify if reads are long read(0) or short read(1). Default behavior is 1 set as short read."
                )
            )
        ] = 1,
        debug_write_diamond_out: Annotated[
            bool,
            typer.Option(
                "--debug-write-diamond-out",
                help="Write DIAMOND TSV to --diamond-output-path (debug/reuse).",
            ),
        ] = False,
        run_diamond: Annotated[
            bool,
            typer.Option(
                "--run-diamond/--no-run-diamond",
                help=(
                    "Run DIAMOND blastx (default). If disabled, "
                    "--diamond-output-path must exist."
                ),
            ),
        ] = True,
        skip_codon_usage: Annotated[
            bool,
            typer.Option(
                "--skip-codon-usage",
                help="Skip codon adaptation/CAI analysis (engineered k-mer scan only).",
            ),
        ] = False,
        debug_write_kraken_out: Annotated[
            bool,
            typer.Option(
                "--debug-write-kraken-out",
                help="Write raw Kraken2 classifications to --kraken-output-path (debug).",
            ),
        ] = False,
        debug_write_kraken_report: Annotated[
            bool,
            typer.Option(
                "--debug-write-kraken-report",
                help="Write Kraken2 --report file (debug).",
            ),
        ] = False,
        run_kraken: Annotated[
            bool,
            typer.Option(
                "--run-kraken/--no-run-kraken",
                help=(
                    "Run Kraken2 as part of the pipeline (default). If disabled, "
                    "--kraken-output-path must point to an existing classifications file."
                ),
            ),
        ] = True,
        # Perseus disabled for now
        # run_perseus: Annotated[
        #     bool,
        #     typer.Option(
        #         "--run-perseus/--no-run-perseus",
        #         help=(
        #             "After the Kraken2 k-mer scan, run Perseus on Natural reads. "
        #             "Reads with perseus_taxid 32630 are marked engineered. "
        #             "Requires perseus on PATH and --perseus-db."
        #         ),
        #     ),
        # ] = False,
        # perseus_db: Annotated[
        #     str | None,
        #     typer.Option(
        #         "--perseus-db",
        #         help=(
        #             "Perseus ETE3 taxonomy directory (from `perseus setup`). "
        #             "Required with --run-perseus."
        #         ),
        #     ),
        # ] = None,
        # perseus_output_path: Annotated[
        #     str | None,
        #     typer.Option(
        #         "--perseus-output-path",
        #         help=(
        #             "Save/load Perseus filter TSV. Required with "
        #             "--debug-write-perseus-out, or with --no-run-perseus when "
        #             "loading a precomputed file (omit both to skip Perseus)."
        #         ),
        #     ),
        # ] = None,
        # debug_write_perseus_out: Annotated[
        #     bool,
        #     typer.Option(
        #         "--debug-write-perseus-out",
        #         help="Write Perseus filter TSV to --perseus-output-path (debug).",
        #     ),
        # ] = False,
        kraken_db_path: Annotated[
            str, typer.Argument(help="Kraken2 database path")
        ] = DEFAULT_DB_PATH,
        threads: Annotated[
            int, typer.Option("--threads", help="Available threads to use.")
        ] = 4,
        background_rate: Annotated[
            float,
            typer.Option(
                "--bg",
                "--background-rate",
                help=(
                    "Null background rate for the engineered-kmer binomial p-value "
                    f"(default: {BACKGROUND_RATE})."
                ),
            ),
        ] = BACKGROUND_RATE,
        pvalue_alpha: Annotated[
            float,
            typer.Option(
                "--pvalue-alpha",
                help=(
                    "Keep window-threshold Engineered hits as engineered only if "
                    f"P-Value ≤ this cutoff (default: {PVALUE_ALPHA})."
                ),
            ),
        ] = PVALUE_ALPHA,
        mem_mode: Annotated[
            int,
            typer.Option(
                "--mem-mode",
                help=(
                    "Memory configuration bins: low(0)=32GB, medium(1)=64GB, "
                    "high(2)=128GB RAM"
                ),
            ),
        ] = 0,
        quiet: Annotated[
            bool,
            typer.Option(
                "--quiet",
                "-q",
                help="Suppress progress bars and INFO logs (warnings/errors only).",
            ),
        ] = False,
        map_hits: Annotated[
            bool,
            typer.Option(
                "--map-hits/--no-map-hits",
                help=(
                    "After screening, BLAST engineered reads against "
                    "--reference-db and attach top reference hits."
                ),
            ),
        ] = False,
        reference_db: Annotated[
            str | None,
            typer.Option(
                "--reference-db",
                help=(
                    "Nucleotide FASTA (blastn subject) for engineered-read "
                    "mapping. Required with --map-hits or --visualize."
                ),
            ),
        ] = None,
        visualize: Annotated[
            bool,
            typer.Option(
                "--visualize/--no-visualize",
                help=(
                    "After screening (and mapping), write pileup figures under "
                    "--viz-dir. Implies --map-hits when mapping has not run yet "
                    "(map already fills BLAST/feature report columns)."
                ),
            ),
        ] = False,
        viz_dir: Annotated[
            str | None,
            typer.Option(
                "--viz-dir",
                help=(
                    "Output directory for pileup PNGs. Required with "
                    "--visualize (default: <report>_viz beside the report)."
                ),
            ),
        ] = None,
        viz_top_n: Annotated[
            int,
            typer.Option(
                "--viz-top-n",
                help="Keep the sample's top-N reference plasmids by read coverage.",
            ),
        ] = 5,
        viz_min_bitscore: Annotated[
            float,
            typer.Option(
                "--viz-min-bitscore",
                help="Minimum blastn bitscore retained for engineered-read hits.",
            ),
        ] = 50.0,
        blastn_binary: Annotated[
            str,
            typer.Option(
                "--blastn",
                help="blastn executable on PATH (or absolute path).",
            ),
        ] = "blastn",
        blast_output_path: Annotated[
            str | None,
            typer.Option(
                "--blast-output-path",
                help="Optional path to save raw blastn TSV from --map-hits.",
            ),
        ] = None,
        genbank_path: Annotated[
            str | None,
            typer.Option(
                "--genbank",
                help=(
                    "Optional local GenBank file for map-hits overlapping "
                    "features and pileup tracks."
                ),
            ),
        ] = None,
        genbank_cache_dir: Annotated[
            str | None,
            typer.Option(
                "--genbank-cache-dir",
                help=(
                    "Local GenBank cache for map-hits overlapping features and "
                    f"pileup tracks (default: {default_genbank_cache_dir()}). "
                    "Prefill with `plasmidscreen prefetch-genbank`."
                ),
            ),
        ] = None,
        fetch_genbank: Annotated[
            bool,
            typer.Option(
                "--fetch-genbank/--no-fetch-genbank",
                help=(
                    "With --map-hits / --visualize, download missing GenBank "
                    "records from NCBI into the cache (default: off; prefer "
                    "prefetch-genbank)."
                ),
            ),
        ] = False,
        ncbi_api_key: Annotated[
            str | None,
            typer.Option(
                "--ncbi-api-key",
                help="Optional NCBI API key used when --fetch-genbank is set.",
            ),
        ] = None,
) -> None:
    """Screening of engineered sequences."""
    configure_logging(quiet=quiet, force=True)
    if not run_kraken and not kraken_output_path:
        raise typer.BadParameter(
            "--kraken-output-path is required when --no-run-kraken is set"
        )
    if debug_write_kraken_out and not kraken_output_path:
        raise typer.BadParameter(
            "--kraken-output-path is required when --debug-write-kraken-out is set"
        )
    if not skip_codon_usage:
        if run_diamond and not diamond_db:
            raise typer.BadParameter(
                "--diamond-db is required when codon usage is enabled "
                "(use --skip-codon-usage to disable)."
            )
        if not run_diamond and not diamond_output_path:
            raise typer.BadParameter(
                "--diamond-output-path is required when --no-run-diamond is set"
            )
        if debug_write_diamond_out and not diamond_output_path:
            raise typer.BadParameter(
                "--diamond-output-path is required when --debug-write-diamond-out is set"
            )
    # Perseus disabled for now
    # if run_perseus and not perseus_db:
    #     raise typer.BadParameter(
    #         "--perseus-db is required when --run-perseus is set"
    #     )
    # if debug_write_perseus_out and not perseus_output_path:
    #     raise typer.BadParameter(
    #         "--perseus-output-path is required when --debug-write-perseus-out is set"
    #     )
    if visualize:
        map_hits = True
    if map_hits and not reference_db:
        raise typer.BadParameter(
            "--reference-db is required when --map-hits or --visualize is set"
        )
    if visualize and viz_dir is None:
        viz_dir = str(Path(output_report_path).with_suffix("")) + "_viz"
    try:
        result = run_screen(
            fasta_file=fasta_file,
            seq_type=seq_type,
            kraken_db=kraken_db_path,
            engineered_report_path=output_report_path,
            kraken_output_path=kraken_output_path,
            threads=threads,
            window_size=window_size,
            engineered_kmer_threshold=engineered_k_mer_threshold,
            kmer_size=kmer_size,
            minimizer_size=minimizer_size,
            codon_usage_output_path=codon_usage_output,
            codon_usage_dir=codon_usage_dir,
            codon_cai_engineered_threshold=codon_cai_engineered_threshold,
            run_codon_usage=not skip_codon_usage,
            debug_write_kraken_output=debug_write_kraken_out,
            debug_write_kraken_report=debug_write_kraken_report,
            run_kraken=run_kraken,
            diamond_db=diamond_db,
            diamond_output_path=diamond_output_path,
            debug_write_diamond_output=debug_write_diamond_out,
            run_diamond=run_diamond,
            # Perseus disabled for now
            # run_perseus=run_perseus,
            # perseus_db=perseus_db,
            # perseus_output_path=perseus_output_path,
            # debug_write_perseus_output=debug_write_perseus_out,
            background_rate=background_rate,
            pvalue_alpha=pvalue_alpha,
            memory_mode=mem_mode,
            quiet_mode=quiet,
        )
    except (PlasmidScreenError, FileNotFoundError, OSError, ValueError) as exc:
        _exit_on_error(exc)
        return  # pragma: no cover — typer.Exit above

    if not quiet:
        n_eng = result.overall_synthetic_count
        n_total = len(result.per_read)
        logger.info(
            "Screen complete: %d/%d read(s) engineered overall → %s",
            n_eng,
            n_total,
            output_report_path,
        )

    if map_hits or visualize:
        try:
            if visualize:
                assert reference_db is not None and viz_dir is not None
                result, viz_report = map_and_visualize_engineered(
                    result,
                    query_fasta=fasta_file,
                    reference_db=reference_db,
                    viz_dir=viz_dir,
                    threads=threads,
                    top_n=viz_top_n,
                    min_bitscore=viz_min_bitscore,
                    blast_output_path=blast_output_path,
                    blastn_binary=blastn_binary,
                    genbank_path=genbank_path,
                    genbank_cache_dir=genbank_cache_dir,
                    fetch_genbank=fetch_genbank,
                    ncbi_api_key=ncbi_api_key,
                )
                if not quiet:
                    report_note = (
                        f"; report {result.engineered_report_path}"
                        if result.engineered_report_path
                        else ""
                    )
                    logger.info(
                        "Visualization complete: %d pileup(s) → %s%s",
                        len(viz_report.pileup_paths),
                        viz_report.output_dir,
                        report_note,
                    )
            else:
                assert reference_db is not None
                result = map_engineered_reads(
                    result,
                    query_fasta=fasta_file,
                    reference_db=reference_db,
                    threads=threads,
                    top_n=viz_top_n,
                    min_bitscore=viz_min_bitscore,
                    output_path=blast_output_path,
                    blastn_binary=blastn_binary,
                    genbank_path=genbank_path,
                    genbank_cache_dir=genbank_cache_dir,
                    fetch_genbank=fetch_genbank,
                    ncbi_api_key=ncbi_api_key,
                )
                if not quiet:
                    logger.info(
                        "Mapped engineered reads: %d top reference(s)",
                        len(result.top_engineered_references),
                    )
        except (PlasmidScreenError, FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
            _exit_on_error(exc)


@app.command()
def visualize(
        fasta_file: Annotated[
            str,
            typer.Argument(help="Input FASTA/FASTQ used for screening (plain or .gz)"),
        ],
        engineered_report: Annotated[
            str,
            typer.Argument(help="Engineered report TSV from a prior `screen` run"),
        ],
        reference_db: Annotated[
            str,
            typer.Argument(help="Nucleotide FASTA (blastn subject) for engineered-read mapping"),
        ],
        viz_dir: Annotated[
            str,
            typer.Argument(help="Output directory for pileup PNGs"),
        ],
        threads: Annotated[
            int, typer.Option("--threads", help="Threads for blastn")
        ] = 4,
        top_n: Annotated[
            int,
            typer.Option(
                "--top-n",
                help="Keep the sample's top-N reference plasmids by read coverage.",
            ),
        ] = 5,
        min_bitscore: Annotated[
            float,
            typer.Option(
                "--min-bitscore",
                help="Minimum blastn bitscore retained for engineered-read hits.",
            ),
        ] = 50.0,
        blastn_binary: Annotated[
            str,
            typer.Option(
                "--blastn",
                help="blastn executable on PATH (or absolute path).",
            ),
        ] = "blastn",
        blast_output_path: Annotated[
            str | None,
            typer.Option(
                "--blast-output-path",
                help="Optional path to save raw blastn TSV.",
            ),
        ] = None,
        genbank_path: Annotated[
            str | None,
            typer.Option(
                "--genbank",
                help=(
                    "Optional local GenBank file for overlapping features and "
                    "pileup tracks."
                ),
            ),
        ] = None,
        genbank_cache_dir: Annotated[
            str | None,
            typer.Option(
                "--genbank-cache-dir",
                help=(
                    "Local GenBank cache for overlapping features and pileup "
                    f"tracks (default: {default_genbank_cache_dir()}). "
                    "Prefill with `plasmidscreen prefetch-genbank`."
                ),
            ),
        ] = None,
        fetch_genbank: Annotated[
            bool,
            typer.Option(
                "--fetch-genbank/--no-fetch-genbank",
                help=(
                    "Download missing GenBank records from NCBI into the cache "
                    "(default: off; prefer prefetch-genbank for bulk installs)."
                ),
            ),
        ] = False,
        ncbi_api_key: Annotated[
            str | None,
            typer.Option(
                "--ncbi-api-key",
                help="Optional NCBI API key used when --fetch-genbank is set.",
            ),
        ] = None,
        threshold: Annotated[
            int,
            typer.Option(
                "--threshold",
                help="K-mer threshold recorded on reconstructed ScreenResult metadata.",
            ),
        ] = 10,
        window_size: Annotated[
            int,
            typer.Option(
                "--window_size",
                help="Window size recorded on reconstructed ScreenResult metadata.",
            ),
        ] = 200,
        quiet: Annotated[
            bool,
            typer.Option(
                "--quiet",
                "-q",
                help="Suppress INFO logs (warnings/errors only).",
            ),
        ] = False,
) -> None:
    """Map engineered reads from a prior screen report and write pileup visualizations."""
    configure_logging(quiet=quiet, force=True)
    try:
        result = screen_result_from_engineered_report(
            engineered_report,
            engineered_kmer_threshold=threshold,
            engineered_kmer_window_size=window_size,
        )
        _enriched, viz_report = map_and_visualize_engineered(
            result,
            query_fasta=fasta_file,
            reference_db=reference_db,
            viz_dir=viz_dir,
            threads=threads,
            top_n=top_n,
            min_bitscore=min_bitscore,
            blast_output_path=blast_output_path,
            blastn_binary=blastn_binary,
            genbank_path=genbank_path,
            genbank_cache_dir=genbank_cache_dir,
            fetch_genbank=fetch_genbank,
            ncbi_api_key=ncbi_api_key,
        )
    except (PlasmidScreenError, FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        _exit_on_error(exc)
        return  # pragma: no cover

    if not quiet:
        report_note = (
            f"; report {result.engineered_report_path}"
            if result.engineered_report_path
            else ""
        )
        logger.info(
            "Visualization complete: %d/%d engineered read(s), %d pileup(s) → %s%s",
            result.overall_synthetic_count,
            len(result.per_read),
            len(viz_report.pileup_paths),
            viz_report.output_dir,
            report_note,
        )


@app.command("prefetch-genbank")
def prefetch_genbank(
        reference_fasta: Annotated[
            str,
            typer.Argument(
                help="Nucleotide reference FASTA (same DB used for --reference-db / visualize)",
            ),
        ],
        cache_dir: Annotated[
            str | None,
            typer.Option(
                "--cache-dir",
                help=(
                    "Directory for cached .gb files "
                    f"(default: {default_genbank_cache_dir()})"
                ),
            ),
        ] = None,
        ncbi_api_key: Annotated[
            str | None,
            typer.Option(
                "--ncbi-api-key",
                help="Optional NCBI API key (higher eutils rate limit).",
            ),
        ] = None,
        quiet: Annotated[
            bool,
            typer.Option(
                "--quiet",
                "-q",
                help="Suppress INFO logs (warnings/errors only).",
            ),
        ] = False,
) -> None:
    """Download GenBank records for NCBI accessions in a reference FASTA into the local cache.

    Run this once (networked) before airgapped visualization. ``screen --visualize`` /
    ``visualize`` then overlay features from the cache without contacting NCBI.
    """
    configure_logging(quiet=quiet, force=True)
    try:
        result = prefetch_genbank_from_fasta(
            reference_fasta,
            cache_dir=cache_dir,
            api_key=ncbi_api_key,
        )
    except (PlasmidScreenError, FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        _exit_on_error(exc)
        return  # pragma: no cover

    if not quiet:
        logger.info(
            "GenBank prefetch done: requested=%d cached=%d missing=%d fetched_this_run=%d → %s",
            result.ids_requested,
            result.ids_cached,
            result.ids_missing,
            result.ids_fetched_this_run,
            result.cache_dir,
        )


def _parse_taxids(taxids: str | None, taxids_file: Path | None) -> list[str]:
    found: set[str] = set()
    if taxids:
        found.update(t.strip() for t in taxids.split(",") if t.strip())
    if taxids_file:
        try:
            text = taxids_file.read_text(encoding="utf-8")
        except OSError as exc:
            raise typer.BadParameter(
                f"Could not read --taxids-file {taxids_file}: {exc}"
            ) from exc
        for line in text.splitlines():
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
        taxids: str | None = typer.Option(
            None, "--taxids", help="Comma-separated NCBI taxonomy IDs"
        ),
        taxids_file: Path | None = typer.Option(
            None, "--taxids-file", help="One taxid per line"
        ),
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
            GeneSet.NUCLEAR,
            "--gene-set",
            help="CSDB gene set: nuclear, ribosomal, mitochondrial, or plastid",
        ),
        quiet: bool = typer.Option(
            False,
            "--quiet",
            "-q",
            help="Suppress INFO logs (warnings/errors only).",
        ),
) -> None:
    """Build codon usage tables from the Codon Statistics Database for CAI scoring."""
    configure_logging(quiet=quiet, force=True)
    data_dir = output_dir or default_codon_usage_dir()
    archive = csdb_archive or default_csdb_archive_path()

    taxid_list = _parse_taxids(taxids, taxids_file)

    if not quiet:
        if taxid_list:
            logger.info(
                "Building codon reference at %s for %d taxid(s) from CSDB archive %s",
                data_dir,
                len(taxid_list),
                archive,
            )
        else:
            logger.info(
                "Building codon reference at %s for all taxids in CSDB archive %s",
                data_dir,
                archive,
            )

    try:
        result = build_codon_reference(
            data_dir,
            taxid_list or None,
            csdb_archive=archive,
            download_csdb=not no_download_csdb,
            gene_set=gene_set,
        )
    except (PlasmidScreenError, FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        _exit_on_error(exc)

    if not quiet:
        logger.info(
            "Done. added=%d skipped=%d failed=%d → %s",
            len(result.taxids_added),
            len(result.taxids_skipped),
            len(result.taxids_failed),
            result.data_dir / "codon_tables.json",
        )


@app.callback()
def main() -> None:
    pass


if __name__ == "__main__":
    app()
