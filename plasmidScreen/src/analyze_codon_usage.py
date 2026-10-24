"""Codon adaptation analysis via DIAMOND blastx and pre-built CSDB reference tables."""
from __future__ import annotations

import logging
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable, Optional, Union

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from plasmidScreen.src.codon_usage.codon_usage_db import CodonUsageStore, default_codon_usage_dir
from plasmidScreen.lib.models import CodonAdaptationResult, CodonAdaptationRead
from plasmidScreen.lib.seqio_util import (
    is_fastq_path,
    open_text_maybe_gzip,
    sequence_format,
)
from plasmidScreen.lib.types import CdsOrf, KrakenReadInfo
from plasmidScreen.lib.diamond_host_taxonomy import (
    infer_orfs_and_host_taxids,
    resolve_diamond_lines,
)


def parse_kraken_lines(lines: Iterable[str]) -> dict[str, KrakenReadInfo]:
    """Parse Kraken2 classification lines (same layout as :func:`parse_kraken_file`)."""
    kraken_data: dict[str, KrakenReadInfo] = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        status = parts[0]
        read_id = parts[1]
        taxid = parts[2]
        try:
            r_len = int(parts[3])
        except ValueError:
            continue
        k_info = parts[-1] if len(parts) >= 5 else ""
        kraken_data[read_id] = (status, taxid, r_len, k_info)

    return kraken_data


def parse_kraken_file(path: Union[str, Path]) -> dict[str, KrakenReadInfo]:
    """Parse a Kraken2 classifications file into read_id → fields."""
    text = Path(path).read_text(encoding="utf-8")
    return parse_kraken_lines(text.splitlines())


def compute_cai(cds_seq: str, weights: dict[str, float]) -> float:
    """Sharp & Li CAI: geometric mean of per-codon adaptiveness weights (in-frame triplets)."""
    log_w_sum = 0.0
    codon_count = 0
    for i in range(0, len(cds_seq) - 2, 3):
        codon = cds_seq[i: i + 3]
        if codon in weights and weights[codon] > 0:
            log_w_sum += math.log(weights[codon])
            codon_count += 1
    return math.exp(log_w_sum / codon_count) if codon_count > 0 else 0.0


def _iter_codon_reads(
        results: CodonAdaptationResult | Iterable[CodonAdaptationRead],
) -> Iterable[CodonAdaptationRead]:
    if isinstance(results, CodonAdaptationResult):
        return results.labels
    return results


def codon_adaptation_to_tsv_lines(
        results: CodonAdaptationResult | Iterable[CodonAdaptationRead],
) -> list[str]:
    """Format codon results as TSV lines (header + one row per result)."""
    header = (
        "Read_ID\tCDS_Strand\tCDS_Range\tHost_TaxID\t"
        "Reference_TaxID\tCDS_Len_bp\tCAI_vs_Host"
    )
    lines = [header]
    for r in _iter_codon_reads(results):
        cai = f"{r.cai_vs_host:.4f}" if r.cai_vs_host is not None else "NA"
        ref = r.reference_taxid or "NA"
        lines.append(
            f"{r.read_id}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
            f"{r.host_taxid}\t{ref}\t{r.cds_len_bp}\t{cai}"
        )
    return lines


def write_codon_adaptation_results_tsv(
        output_path: Union[str, Path],
        results: CodonAdaptationResult | Iterable[CodonAdaptationRead],
        *,
        cai_engineered_threshold: Optional[float] = None,
) -> str:
    """
    Write codon adaptation results to a TSV file.

    When ``cai_engineered_threshold`` is set, adds columns for threshold-based
    engineered flagging (used by the full ``screen`` workflow).
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        if cai_engineered_threshold is None:
            for line in codon_adaptation_to_tsv_lines(results):
                f.write(line + "\n")
        else:
            header = (
                "Read_ID\tCDS_Strand\tCDS_Range\tHost_TaxID\t"
                "Reference_TaxID\tCDS_Len_bp\tCAI_vs_Host\t"
                "Codon_CAI_Threshold\tEngineered_By_Codon_CAI"
            )
            f.write(header + "\n")
            for r in _iter_codon_reads(results):
                cai = r.cai_vs_host
                cai_s = f"{cai:.4f}" if cai is not None else "NA"
                engineered = (
                    "NA" if cai is None else str(cai > cai_engineered_threshold)
                )
                ref = r.reference_taxid or "NA"
                f.write(
                    f"{r.read_id}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
                    f"{r.host_taxid}\t{ref}\t{r.cds_len_bp}\t"
                    f"{cai_s}\t{cai_engineered_threshold}\t{engineered}\n"
                )
    return str(out)


def write_codon_adaptation_tsv(
        output_path: Union[str, Path],
        reads_path: Union[str, Path],
        *,
        diamond_db: Union[str, Path] | None = None,
        diamond_threads: int = 4,
        mem_mode: int = 0,
        run_diamond: bool = True,
        diamond_output_path: Union[str, Path] | None = None,
        debug_write_diamond_output: bool = False,
        codon_usage_store: Optional[CodonUsageStore] = None,
        codon_usage_dir: Optional[Union[str, Path]] = None,
        include_read_ids: Optional[set[str]] = None,
        cai_engineered_threshold: Optional[float] = None,
) -> tuple[str, CodonAdaptationResult]:
    """Run codon adaptation analysis and write results to a TSV file."""
    results, _diamond_path = analyze_codon_adaptation(
        reads_path,
        diamond_db=diamond_db,
        diamond_threads=diamond_threads,
        mem_mode=mem_mode,
        run_diamond=run_diamond,
        diamond_output_path=diamond_output_path,
        debug_write_diamond_output=debug_write_diamond_output,
        codon_usage_store=codon_usage_store,
        codon_usage_dir=codon_usage_dir,
        include_read_ids=include_read_ids,
    )
    path = write_codon_adaptation_results_tsv(
        output_path,
        results,
        cai_engineered_threshold=cai_engineered_threshold,
    )
    return path, results


def write_subset_reads(
        reads_path: Path,
        include_read_ids: set[str],
        dest: Path,
) -> int:
    """
    Write only ``include_read_ids`` from ``reads_path`` into ``dest``.

    Returns the number of records written. Preserves FASTA vs FASTQ format.
    Supports gzip-compressed inputs (``.gz``); ``dest`` is always plain text.
    """
    written = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    if is_fastq_path(reads_path):
        with open_text_maybe_gzip(reads_path) as in_h, open(
            dest, "w", encoding="utf-8"
        ) as out_h:
            for title, seq, qual in FastqGeneralIterator(in_h):
                read_id = title.split(None, 1)[0]
                if read_id in include_read_ids:
                    out_h.write(f"@{title}\n{seq}\n+\n{qual}\n")
                    written += 1
    else:
        with open_text_maybe_gzip(reads_path) as in_h, open(
            dest, "w", encoding="utf-8"
        ) as out_h:
            for title, seq in SimpleFastaParser(in_h):
                read_id = title.split(None, 1)[0]
                if read_id in include_read_ids:
                    out_h.write(f">{title}\n{seq}\n")
                    written += 1
    return written


def analyze_codon_adaptation(
        reads_path: Union[str, Path],
        *,
        diamond_db: Union[str, Path] | None = None,
        diamond_threads: int = 4,
        mem_mode: int = 0,
        run_diamond: bool = True,
        diamond_output_path: Union[str, Path] | None = None,
        debug_write_diamond_output: bool = False,
        codon_usage_store: Optional[CodonUsageStore] = None,
        codon_usage_dir: Optional[Union[str, Path]] = None,
        include_read_ids: Optional[set[str]] = None,
        seq_type: int = 1,
) -> tuple[CodonAdaptationResult, Path | None]:
    """
    Score codon adaptation (CAI) for reads using DIAMOND and a pre-built CSDB reference.

    Pipeline per read with DIAMOND hits:

    1. Run DIAMOND blastx (``--min-orf``) or load a saved outfmt 6 TSV.
       When ``include_read_ids`` is set and DIAMOND is run, only those reads are
       written to a temporary query file (avoids blasting Synthetic reads).
    2. Merge hit ``qstart``/``qend`` into ORF intervals; pick the longest interval.
    3. Infer host taxid from hit ``staxids`` (majority over alignments).
    4. Resolve host taxid to a reference codon table (NCBI lineage fallback).
    5. Compute CAI on the ORF nucleotide slice vs host weights.

    Parameters
    ----------
    reads_path
        FASTA or FASTQ path (plain or ``.gz``).
    diamond_db
        DIAMOND protein database (``.dmnd``). Required when ``run_diamond=True``.
    run_diamond
        If False, load alignments from ``diamond_output_path`` instead of running DIAMOND.
    diamond_output_path
        Path to save (with ``debug_write_diamond_output=True``) or load DIAMOND TSV.
    debug_write_diamond_output
        Write DIAMOND outfmt 6 to ``diamond_output_path`` after blastx.
    codon_usage_dir / codon_usage_store
        Pre-built ``codon_tables.json`` (and optional ``taxonomy_parents.json``).
    include_read_ids
        If set, only these read IDs are scored (e.g. Natural reads from engineered scan).
    seq_type
        ``0`` long-read/ONT DIAMOND flags; ``1`` short-read flags (default).

    Returns
    -------
    results
        One :class:`~plasmidScreen.lib.models.CodonAdaptationResult` per scored read.
    diamond_path
        Path to saved DIAMOND TSV when written or loaded from disk; else ``None``.
    """
    reads_path = Path(reads_path)

    if codon_usage_store is not None:
        store = codon_usage_store
    else:
        data_dir = Path(codon_usage_dir) if codon_usage_dir else default_codon_usage_dir()
        store = CodonUsageStore.load(data_dir)

    file_format = sequence_format(reads_path)
    tmp_subset: Path | None = None
    diamond_query = reads_path
    seq_source = reads_path

    try:
        # Only blast the reads we will score (Natural subset in full screen).
        if run_diamond and include_read_ids is not None:
            suffix = ".fastq" if file_format == "fastq" else ".fasta"
            tmp_fd, tmp_name = tempfile.mkstemp(
                prefix="plasmidScreen_natural_", suffix=suffix
            )
            # Close fd; we reopen by path for Bio parsers / DIAMOND.
            os.close(tmp_fd)
            tmp_subset = Path(tmp_name)
            n_written = write_subset_reads(reads_path, include_read_ids, tmp_subset)
            logging.info(
                "DIAMOND query subset: %d/%d requested read(s) written to temp file",
                n_written,
                len(include_read_ids),
            )
            if n_written == 0:
                return CodonAdaptationResult(labels=[]), None
            diamond_query = tmp_subset
            seq_source = tmp_subset

        diamond_lines, diamond_path = resolve_diamond_lines(
            diamond_query,
            diamond_db,
            mem_mode=mem_mode,
            run_diamond=run_diamond,
            output_path=diamond_output_path,
            debug_write_output=debug_write_diamond_output,
            threads=diamond_threads,
            sequencing_technology=seq_type,
        )
        orfs_by_read, host_by_read = infer_orfs_and_host_taxids(diamond_lines)
        pending: list[tuple[str, CdsOrf, str]] = []
        n_seen = 0
        n_no_host = 0
        n_no_orf = 0

        with open_text_maybe_gzip(seq_source) as seq_handle:
            for record in SeqIO.parse(seq_handle, file_format):
                read_id = record.id
                if include_read_ids is not None and read_id not in include_read_ids:
                    continue
                n_seen += 1
                host_taxid = host_by_read.get(read_id)
                if not host_taxid:
                    n_no_host += 1
                    continue

                orfs = orfs_by_read.get(read_id) or []
                if not orfs:
                    n_no_orf += 1
                    continue

                best = max(orfs, key=lambda o: (o.length_bp, o.end))
                seq_str = str(record.seq).upper()
                cds_seq = seq_str[best.start: best.end]
                if len(cds_seq) < 3:
                    n_no_orf += 1
                    continue

                cds: CdsOrf = {
                    "strand": best.strand,
                    "start": best.start,
                    "end": best.end,
                    "length": len(cds_seq) // 3,
                    "seq": cds_seq,
                }
                pending.append((read_id, cds, host_taxid))

        if n_no_host or n_no_orf:
            logging.info(
                "Codon CAI candidates: %d considered, %d pending; "
                "skipped no-host=%d no-orf=%d",
                n_seen,
                len(pending),
                n_no_host,
                n_no_orf,
            )

        host_taxids = {p[2] for p in pending if p[2] not in ("0", "")}
        store.require_host_taxids(host_taxids)

        results: CodonAdaptationResult = CodonAdaptationResult(labels=[])
        for read_id, cds, host_taxid in pending:
            ref_taxid: Optional[str] = None
            cai_score: Optional[float] = None
            if host_taxid not in ("0", "") and len(cds["seq"]) >= 3:
                ref_taxid = store.resolve_reference_taxid(host_taxid)
                weights = store.get_cai_weights_for_host(host_taxid)
                if weights:
                    cai_score = compute_cai(cds["seq"], weights)
                else:
                    logging.warning(
                        "No CAI weights for host taxid %s (read %s); CAI left unset",
                        host_taxid,
                        read_id,
                    )

            results.labels.append(
                CodonAdaptationRead(
                    read_id=read_id,
                    cds_strand=cds["strand"],
                    cds_start=cds["start"],
                    cds_end=cds["end"],
                    host_taxid=host_taxid,
                    reference_taxid=ref_taxid,
                    cds_len_bp=len(cds["seq"]),
                    cai_vs_host=cai_score,
                    host_taxid_method="diamond",
                )
            )

        logging.info("Codon CAI scored %d read(s)", len(results.labels))
        return results, diamond_path
    finally:
        if tmp_subset is not None and tmp_subset.exists():
            try:
                tmp_subset.unlink(missing_ok=True)
            except OSError as exc:
                logging.warning("Could not remove temp DIAMOND query %s: %s", tmp_subset, exc)

