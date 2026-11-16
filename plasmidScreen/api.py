"""
PlasmidScreen public library API.

Screening and codon scoring are airgapped-safe once ``codon_tables.json`` has been
built offline with :func:`build_codon_reference` or :func:`build_codon_database`.

Optional post-screen steps:
- :func:`map_engineered_reads` — blastn enrichment (top references on ``ScreenResult``)
- :func:`visualize_engineered_mappings` — pileups + mapping JSON (viz mode only)
- :func:`map_and_visualize_engineered` — convenience chain of the two above
- :func:`screen_result_from_engineered_report` — rebuild a minimal ``ScreenResult``
  from a written engineered-report TSV (for standalone viz CLI / reruns)
- :func:`prefetch_genbank_from_fasta` — preload GenBank features into the local cache
  before visualization (viz reads the cache; it does not fetch NCBI by default)
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable, Literal

from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference
from plasmidScreen.src.codon_usage.codon_usage_db import (
    default_codon_usage_dir,
)
from plasmidScreen.lib.genbank_fetch import (
    PrefetchResult,
    default_genbank_cache_dir,
    prefetch_genbank_from_fasta,
)
from plasmidScreen.lib.models import (
    BACKGROUND_RATE,
    BuildCodonReferenceResult,
    CodonAdaptationResult,
    EngineeredScanResult,
    MEMORY_CONFIG,
    ReadEngineeringLabel,
    ReadFlagDetail,
    ScreenResult,
    VisualizationReport,
)
from plasmidScreen.lib.types import GeneSet
from plasmidScreen.src.plasmidScreen import Workflow

__all__ = [
    "build_codon_database",
    "run_screen",
    "map_engineered_reads",
    "visualize_engineered_mappings",
    "map_and_visualize_engineered",
    "screen_result_from_engineered_report",
    "prefetch_genbank_from_fasta",
    "default_genbank_cache_dir",
    "PrefetchResult",
    "ScreenResult",
    "VisualizationReport",
    "BuildCodonReferenceResult",
]


def run_screen(
        fasta_file: str | Path,
        seq_type: int,
        kraken_db: str | Path,
        engineered_report_path: str | Path | None = None,
        kraken_output_path: str | Path | None = None,
        threads: int = 4,
        window_size: int = 200,
        engineered_kmer_threshold: int = 10,
        codon_usage_dir: str | Path | None = None,
        run_kraken: bool = True,
        debug_write_kraken_output: bool = False,
        debug_write_kraken_report: bool = False,
        run_codon_usage: bool = True,
        codon_usage_output_path: str | Path | None = None,
        codon_cai_engineered_threshold: float = 0.6,
        diamond_db: str | Path | None = None,
        diamond_output_path: str | Path | None = None,
        debug_write_diamond_output: bool = False,
        run_diamond: bool = True,
        # Perseus disabled for now
        # run_perseus: bool = False,
        # perseus_db: str | Path | None = None,
        # perseus_output_path: str | Path | None = None,
        # debug_write_perseus_output: bool = False,
        quiet_mode: bool = True,
        memory_mode: int = 0,
        background_rate: float = BACKGROUND_RATE
) -> ScreenResult:
    """
    Run engineered k-mer screening (Kraken2) and optional codon adaptation (DIAMOND + CSDB).

    Engineered detection uses Kraken2 minimizer blocks (taxid 32630) in a sliding window.
    Codon CAI runs on reads still labeled **Natural** after the k-mer scan.
    (Perseus Natural-read reclassification is currently disabled.)

    Parameters
    ----------
    fasta_file
        Input FASTA/FASTQ path (plain or gzip-compressed ``.gz``).
    seq_type
        Sequence technology type, if short read specify as 1 or set 0 if your input is long read.
    kraken_db
        Kraken2 database directory (required when ``run_kraken=True``).
    engineered_report_path
        If set, writes the engineered k-mer scan TSV to this path.
    kraken_output_path
        Save or load raw Kraken2 classifications (required when ``run_kraken=False`` or
        when ``debug_write_kraken_output=True``).
    threads
        Number of threads to be used for the tool for engineered k-mer scanning to codon optimization,
        unless specific thread counts are given to diamond
    run_kraken
        Run Kraken2 in-process; if False, ``kraken_output_path`` must point to existing output.
    debug_write_kraken_report
        Debug option to write kraken report
    debug_write_kraken_output
        Debug to write the output of kraken
    codon_usage_dir
        Directory with ``codon_tables.json``. Required when ``run_codon_usage=True``.
    run_codon_usage
        If False, skip DIAMOND and CAI (engineered scan only).
    codon_usage_output_path
        Optional path for codon adaptation TSV (Natural reads only).
    codon_cai_engineered_threshold
        If set, reads with CAI above this value get ``engineered_by_codon_cai=True`` and
        may be marked ``engineered_overall=True`` on :class:`~plasmidScreen.lib.models.ReadFlagDetail`
        even when the k-mer scan labeled them Natural. Also adds columns to the codon TSV when written.
    engineered_kmer_threshold
        Min engineered (32630) k-mers in a window to label a read Synthetic by k-mer scan.
    window_size
        Sliding window size (bp) for the k-mer scan.
    diamond_db
        DIAMOND protein database (``.dmnd``). Required when ``run_codon_usage=True`` and
        ``run_diamond=True``.
    diamond_output_path
        Save or load DIAMOND outfmt 6 TSV (for ``debug_write_diamond_output`` or
        ``run_diamond=False``).
    debug_write_diamond_output
        Persist DIAMOND alignments to ``diamond_output_path``.
    run_diamond
        Run DIAMOND blastx; if False, load precomputed TSV from ``diamond_output_path``.
    quiet_mode
        If True (library default), suppress progress bars. Pair with logging level
        WARNING+ for a quiet run; the CLI uses ``--quiet`` for both.
    memory_mode
        specify the memory configuration for deploying Plasmidscreen with a few discrete bins:
            low(0)    - 32GB  RAM
            medium(1) - 64GB  RAM
            high(2)   - 128GB RAM
    background_rate
        probability rate of how likely any given kmer is classified as engineered.

    Returns
    -------
    ScreenResult
        Includes ``per_read`` with ``engineered_overall`` / ``overall_label`` per read,
        run-level ``overall_synthetic_count`` / ``engineered_read_ids``, and stored
        threshold values used for the combined decision.
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
    # Perseus disabled for now
    # if run_perseus and perseus_db is None:
    #     raise ValueError("perseus_db is required when run_perseus=True.")
    # if debug_write_perseus_output and perseus_output_path is None:
    #     raise ValueError(
    #         "perseus_output_path is required when debug_write_perseus_output=True."
    #     )
    if not memory_mode in MEMORY_CONFIG.values():
        raise ValueError("memory_mode was expecting a value between 0-2 where:"
                         "\nlow(0) <= 32GB  RAM "
                         "\nmedium(1) <= 64GB  RAM "
                         "\nhigh(2) <= 128GB RAM")
    workflow = Workflow(
        fasta_file=str(fasta_file),
        seq_type=seq_type,
        output_report_path= str(engineered_report_path) if engineered_report_path else None,
        kraken_db=str(kraken_db),
        threads=threads,
        kraken_raw_output=str(kraken_output_path) if kraken_output_path else None,
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
        diamond_output_path=str(diamond_output_path) if diamond_output_path else None,
        debug_write_diamond_output=debug_write_diamond_output,
        run_diamond=run_diamond,
        # Perseus disabled for now
        # run_perseus=run_perseus,
        # perseus_db=str(perseus_db) if perseus_db else None,
        # perseus_output_path=str(perseus_output_path) if perseus_output_path else None,
        # debug_write_perseus_output=debug_write_perseus_output,
        quiet_mode=quiet_mode,
        mem_mode=memory_mode,
        background_rate=background_rate
    )
    return workflow.run()


def map_engineered_reads(
        screen_result: ScreenResult,
        *,
        query_fasta: str | Path,
        reference_db: str | Path,
        threads: int = 4,
        top_n: int = 5,
        min_bitscore: float = 50.0,
        max_evalue: float | None = None,
        database_source: str | None = None,
        output_path: str | Path | None = None,
        run_blast: bool = True,
        blastn_binary: str = "blastn",
        blast_task: str = "blastn",
) -> ScreenResult:
    """
    BLAST engineered reads against a nucleotide reference and attach top hits.

    Updates ``ScreenResult`` with per-read ``database_match_*`` fields and
    sample-wide ``top_engineered_references``. Does **not** generate pileups or
    ``engineered_mapping_hits`` — call :func:`visualize_engineered_mappings` for that.
    """
    from plasmidScreen.src.post_hit_analysis import post_analysis_report

    return post_analysis_report(
        screen_result,
        query_fasta=query_fasta,
        reference_db=reference_db,
        threads=threads,
        top_n=top_n,
        min_bitscore=min_bitscore,
        max_evalue=max_evalue,
        database_source=database_source,
        output_path=output_path,
        run_blast=run_blast,
        blastn_binary=blastn_binary,
        blast_task=blast_task,
    )


def visualize_engineered_mappings(
        screen_result: ScreenResult,
        output_dir: str | Path,
        *,
        genbank_path: str | Path | None = None,
        fetch_genbank: bool = False,
        genbank_cache_dir: str | Path | None = None,
        ncbi_api_key: str | None = None,
        dpi: int = 180,
        write_json: bool = True,
) -> VisualizationReport:
    """
    Optional viz mode: pileup figures + mapping summaries for top engineered refs.

    Uses ``screen_result.top_engineered_references`` (from
    :func:`map_engineered_reads`) when present. Writes PNGs under ``output_dir``,
    attaches ``engineered_mapping_hits`` on ``screen_result``, and optionally
    writes a sidecar JSON mapping summary.

    GenBank feature tracks come from the local cache by default (see
    :func:`prefetch_genbank_from_fasta`). On-demand NCBI fetch is off unless
    ``fetch_genbank=True``.
    """
    from plasmidScreen.src.visualization_code import create_pile_ups

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    pileup_paths = create_pile_ups(
        screen_result,
        out,
        genbank_path=genbank_path,
        fetch_genbank=fetch_genbank,
        genbank_cache_dir=genbank_cache_dir,
        ncbi_api_key=ncbi_api_key,
        dpi=dpi,
        engineered_only=True,
    )

    report_json_path: Path | None = None
    if write_json:
        report_json_path = out / "visualization_report.json"
        payload = {
            "top_engineered_references": [
                {
                    "accession_id": r.accession_id,
                    "n_reads": r.n_reads,
                    "sum_bitscore": r.sum_bitscore,
                }
                for r in screen_result.top_engineered_references
            ],
            "engineered_mapping_hits": [
                {
                    "accession_id": m.accession_id,
                    "target_length": m.target_length,
                    "n_queries": m.n_queries,
                    "n_alignments": m.n_alignments,
                    "query_ids": list(m.query_ids),
                    "covered_bp": m.covered_bp,
                    "mean_coverage": m.mean_coverage,
                    "max_coverage": m.max_coverage,
                    "alignments": [
                        {
                            "query_id": a.query_id,
                            "target_start": a.target_start,
                            "target_end": a.target_end,
                            "sequence_identity": a.sequence_identity,
                            "match_score": a.match_score,
                            "overlapping_features": list(a.overlapping_features),
                        }
                        for a in m.alignments
                    ],
                    "features": [
                        {
                            "label": f.label,
                            "feature_type": f.feature_type,
                            "start": f.start,
                            "end": f.end,
                            "strand": f.strand,
                        }
                        for f in m.features
                    ],
                }
                for m in screen_result.engineered_mapping_hits
            ],
            "pileups": [str(p) for p in pileup_paths],
            "reference_database_source": screen_result.reference_database_source,
        }
        report_json_path.write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )

    return VisualizationReport(
        output_dir=out,
        pileup_paths=tuple(pileup_paths),
        engineered_mapping_hits=screen_result.engineered_mapping_hits,
        top_engineered_references=screen_result.top_engineered_references,
        report_json_path=report_json_path,
    )


def map_and_visualize_engineered(
        screen_result: ScreenResult,
        *,
        query_fasta: str | Path,
        reference_db: str | Path,
        viz_dir: str | Path,
        threads: int = 4,
        top_n: int = 5,
        min_bitscore: float = 50.0,
        max_evalue: float | None = None,
        database_source: str | None = None,
        blast_output_path: str | Path | None = None,
        run_blast: bool = True,
        blastn_binary: str = "blastn",
        blast_task: str = "blastn",
        genbank_path: str | Path | None = None,
        fetch_genbank: bool = False,
        genbank_cache_dir: str | Path | None = None,
        ncbi_api_key: str | None = None,
        dpi: int = 180,
        write_json: bool = True,
) -> tuple[ScreenResult, VisualizationReport]:
    """
    BLAST-enrich engineered reads, then write pileup figures / mapping summaries.

    Convenience wrapper around :func:`map_engineered_reads` followed by
    :func:`visualize_engineered_mappings`. Returns the enriched ``ScreenResult``
    and the ``VisualizationReport``. Prefetch GenBank features beforehand so
    pileup annotation tracks resolve offline from the local cache.
    """
    enriched = map_engineered_reads(
        screen_result,
        query_fasta=query_fasta,
        reference_db=reference_db,
        threads=threads,
        top_n=top_n,
        min_bitscore=min_bitscore,
        max_evalue=max_evalue,
        database_source=database_source,
        output_path=blast_output_path,
        run_blast=run_blast,
        blastn_binary=blastn_binary,
        blast_task=blast_task,
    )
    viz = visualize_engineered_mappings(
        enriched,
        viz_dir,
        genbank_path=genbank_path,
        fetch_genbank=fetch_genbank,
        genbank_cache_dir=genbank_cache_dir,
        ncbi_api_key=ncbi_api_key,
        dpi=dpi,
        write_json=write_json,
    )
    return enriched, viz


def screen_result_from_engineered_report(
        report_path: str | Path,
        *,
        engineered_kmer_threshold: int = 10,
        engineered_kmer_window_size: int = 200,
) -> ScreenResult:
    """
    Rebuild a minimal ``ScreenResult`` from a written engineered-report TSV.

    Used by the ``visualize`` CLI when re-running map/viz without a fresh
    ``run_screen`` call. Only overall engineered / natural labels (and optional
    Methods / P-Value / coverage columns) are restored; codon and BLAST fields
    are empty until :func:`map_engineered_reads` runs.
    """
    path = Path(report_path)
    if not path.is_file():
        raise FileNotFoundError(f"Engineered report not found: {path}")

    labels: list[ReadEngineeringLabel] = []
    per_read: list[ReadFlagDetail] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Read_ID" not in reader.fieldnames:
            raise ValueError(
                f"Engineered report {path} is missing a Read_ID column "
                f"(found: {reader.fieldnames})"
            )
        for row in reader:
            read_id = (row.get("Read_ID") or "").strip()
            if not read_id:
                continue
            raw_label = (row.get("Label") or "").strip()
            engineered = raw_label.lower() in {"engineered", "synthetic"}
            overall_label: Literal["Natural", "Synthetic"] = (
                "Synthetic" if engineered else "Natural"
            )
            methods_raw = (row.get("Methods") or "").strip()
            methods = {
                m.strip()
                for m in methods_raw.split(",")
                if m.strip() and m.strip().upper() != "NA"
            }
            by_kmer = engineered and (
                not methods or "engineered_kmer_scan" in methods
            )
            # Prefer Methods when present; otherwise treat Engineered rows as k-mer.
            if engineered and not methods:
                by_kmer = True
            kmer_label: Literal["Natural", "Synthetic"] = (
                "Synthetic" if by_kmer else "Natural"
            )
            coverage = 0.0
            cov_raw = (row.get("Eng-Coverage Percentage") or "").strip()
            if cov_raw.endswith("%"):
                try:
                    coverage = float(cov_raw[:-1]) / 100.0
                except ValueError:
                    coverage = 0.0
            elif cov_raw:
                try:
                    coverage = float(cov_raw)
                except ValueError:
                    coverage = 0.0
            p_value: float | None = None
            p_raw = (row.get("P-Value") or "").strip()
            if p_raw and p_raw.upper() != "NA":
                try:
                    p_value = float(p_raw)
                except ValueError:
                    p_value = None

            labels.append(
                ReadEngineeringLabel(
                    read_id=read_id,
                    label=kmer_label,
                    eng_kmer_coverage=coverage,
                )
            )
            per_read.append(
                ReadFlagDetail(
                    read_id=read_id,
                    kmer_label=kmer_label,
                    engineered_by_kmer_scan=by_kmer,
                    engineered_overall=engineered,
                    overall_label=overall_label,
                    engineered_by_codon_cai=(
                        True if "codon_cai" in methods else None
                    ),
                    # engineered_by_perseus=(
                    #     True if "perseus" in methods else None
                    # ),
                    eng_kmer_coverage=coverage,
                    p_value=p_value,
                    engineered_kmer_threshold=engineered_kmer_threshold,
                    engineered_kmer_window_size=engineered_kmer_window_size,
                )
            )

    synthetic = sum(1 for lbl in labels if lbl.label == "Synthetic")
    return ScreenResult(
        engineered_scan=EngineeredScanResult(
            labels=labels,
            synthetic_count=synthetic,
            natural_count=len(labels) - synthetic,
        ),
        codon_adaptation=CodonAdaptationResult(),
        engineered_kmer_threshold=engineered_kmer_threshold,
        engineered_kmer_window_size=engineered_kmer_window_size,
        codon_cai_engineered_threshold=None,
        per_read=per_read,
        engineered_report_path=path,
    )


def build_codon_database(
        *,
        output_dir: str | Path | None = None,
        taxids: Iterable[str | int] | None = None,
        taxids_file: str | Path | None = None,
        csdb_archive: str | Path | None = None,
        download_csdb: bool = True,
        gene_set: GeneSet = GeneSet.NUCLEAR,
) -> BuildCodonReferenceResult:
    """
    Serves as an API wrapper for build_codon_reference so users can easily create codon database as needed.

    Build the codon usage reference for CAI scoring. It attempts to use a
    taxids file if provided to grab the codon references. If none is provided, it imports
    **every** taxid in the CSDB archive for ``gene_set``. Downloads CSDB archive if download_csdb is provided as true.
    Writes ``codon_tables.json`` and optionally ``taxonomy_parents.json``
    under ``output_dir`` (default: PlasmidScreen user data ``codon_usage/``).
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

    taxid_list = sorted(resolved) if resolved else None

    return build_codon_reference(
        data_dir,
        taxid_list,
        csdb_archive=csdb_archive,
        download_csdb=download_csdb,
        gene_set=gene_set,
    )
