"""Codon adaptation analysis via DIAMOND blastx and pre-built CSDB reference tables."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Optional, Union

from Bio import SeqIO

from plasmidScreen.src.codon_usage.codon_usage_db import CodonUsageStore, default_codon_usage_dir
from plasmidScreen.lib.models import CodonAdaptationResult
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


def codon_adaptation_to_tsv_lines(results: Iterable[CodonAdaptationResult]) -> list[str]:
    """Format codon results as TSV lines (header + one row per result)."""
    header = (
        "Read_ID\tCDS_Strand\tCDS_Range\tHost_TaxID\t"
        "Reference_TaxID\tCDS_Len_bp\tCAI_vs_Host"
    )
    lines = [header]
    for r in results:
        cai = f"{r.cai_vs_host:.4f}" if r.cai_vs_host is not None else "NA"
        ref = r.reference_taxid or "NA"
        lines.append(
            f"{r.read_id}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
            f"{r.host_taxid}\t{ref}\t{r.cds_len_bp}\t{cai}"
        )
    return lines


def write_codon_adaptation_results_tsv(
        output_path: Union[str, Path],
        results: Iterable[CodonAdaptationResult],
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
            for r in results:
                cai = r.cai_vs_host
                cai_s = f"{cai:.4f}" if cai is not None else "NA"
                engineered = (
                    "NA" if cai is None else str(cai < cai_engineered_threshold)
                )
                ref = r.reference_taxid or "NA"
                f.write(
                    f"{r.read_id}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
                    f"{r.host_taxid}\t{ref}\t{r.cds_len_bp}\t"
                    f"{cai_s}\t{cai_engineered_threshold}\t{engineered}\n"
                )
    return str(out)


def analyze_codon_adaptation(
        reads_path: Union[str, Path],
        *,
        diamond_db: Union[str, Path] | None = None,
        diamond_threads: int = 4,
        run_diamond: bool = True,
        diamond_output_path: Union[str, Path] | None = None,
        debug_write_diamond_output: bool = False,
        codon_usage_store: Optional[CodonUsageStore] = None,
        codon_usage_dir: Optional[Union[str, Path]] = None,
        include_read_ids: Optional[set[str]] = None,
) -> tuple[list[CodonAdaptationResult], Path | None]:
    """
    Score codon adaptation (CAI) for reads using DIAMOND and a pre-built CSDB reference.

    Pipeline per read with DIAMOND hits:

    1. Run DIAMOND blastx (``--min-orf``) or load a saved outfmt 6 TSV.
    2. Merge hit ``qstart``/``qend`` into ORF intervals; pick the longest interval.
    3. Infer host taxid from hit ``staxids`` (majority over alignments).
    4. Resolve host taxid to a reference codon table (NCBI lineage fallback).
    5. Compute CAI on the ORF nucleotide slice vs host weights.

    Parameters
    ----------
    reads_path
        FASTA or FASTQ path.
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

    diamond_lines, diamond_path = resolve_diamond_lines(
        reads_path,
        diamond_db,
        run_diamond=run_diamond,
        output_path=diamond_output_path,
        debug_write_output=debug_write_diamond_output,
        threads=diamond_threads,
    )
    orfs_by_read, host_by_read = infer_orfs_and_host_taxids(diamond_lines)
    file_format = "fastq" if reads_path.suffix.lower() in (".fastq", ".fq") else "fasta"
    pending: list[tuple[str, CdsOrf, str]] = []

    for record in SeqIO.parse(reads_path, file_format):
        read_id = record.id
        if include_read_ids is not None and read_id not in include_read_ids:
            continue
        host_taxid = host_by_read.get(read_id)
        if not host_taxid:
            continue

        orfs = orfs_by_read.get(read_id) or []
        if not orfs:
            continue

        best = max(orfs, key=lambda o: (o.length_bp, o.end))
        seq_str = str(record.seq).upper()
        cds_seq = seq_str[best.start: best.end]
        if len(cds_seq) < 3:
            continue

        cds: CdsOrf = {
            "strand": best.strand,
            "start": best.start,
            "end": best.end,
            "length": len(cds_seq) // 3,
            "seq": cds_seq,
        }
        pending.append((read_id, cds, host_taxid))

    host_taxids = {p[2] for p in pending if p[2] not in ("0", "")}
    store.require_host_taxids(host_taxids)

    results: list[CodonAdaptationResult] = []
    for read_id, cds, host_taxid in pending:
        ref_taxid: Optional[str] = None
        cai_score: Optional[float] = None
        if host_taxid not in ("0", "") and len(cds["seq"]) >= 3:
            ref_taxid = store.resolve_reference_taxid(host_taxid)
            weights = store.get_cai_weights_for_host(host_taxid)
            if weights:
                cai_score = compute_cai(cds["seq"], weights)

        results.append(
            CodonAdaptationResult(
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

    return results, diamond_path

