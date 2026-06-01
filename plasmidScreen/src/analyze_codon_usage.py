"""Codon adaptation analysis using a pre-built reference store (airgapped-safe)."""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional, Union

from Bio import SeqIO
from Bio.Seq import Seq

from plasmidScreen.lib.codon_usage_db import CodonUsageStore, default_codon_usage_dir
from plasmidScreen.lib.models import CodonAdaptationResult
from plasmidScreen.lib.types import CdsOrf, KrakenReadInfo, PendingCodonRead
from plasmidScreen.lib.diamond_host_taxonomy import (
    infer_orfs_and_host_taxids,
    resolve_diamond_lines,
)


def parse_kraken_file(kraken_path: str | Path) -> dict[str, KrakenReadInfo]:
    """Parse Kraken2 output into read_id -> (status, taxid, length, kmer_info)."""
    kraken_data: dict[str, KrakenReadInfo] = {}
    with open(kraken_path, "r", buffering=1024 * 1024) as f:
        for line in f:
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


def parse_kraken_lines(lines: Iterable[str]) -> dict[str, KrakenReadInfo]:
    """Parse Kraken2 output lines into read_id -> (status, taxid, length, kmer_info)."""
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


def expand_kmer_taxids(kmer_info: str, read_length: int, k: int = 35) -> list[str]:
    kmers: list[str] = []
    for block in kmer_info.split(" "):
        if not block:
            continue
        taxid, count = block.split(":")
        kmers.extend([taxid] * int(count))

    expected_kmers = read_length - k + 1
    if len(kmers) < expected_kmers:
        kmers.extend(["0"] * (expected_kmers - len(kmers)))
    return kmers[:expected_kmers]


def map_nucleotides_to_taxids_fast(kmers: list[str], read_length: int, k: int = 35) -> list[str]:
    nuc_taxids: list[str] = []
    num_kmers = len(kmers)
    current_counts: dict[str, int] = {}

    for i in range(read_length):
        if i < num_kmers:
            tax_in = kmers[i]
            current_counts[tax_in] = current_counts.get(tax_in, 0) + 1

        leaving_idx = i - k
        if leaving_idx >= 0 and leaving_idx < num_kmers:
            tax_out = kmers[leaving_idx]
            current_counts[tax_out] -= 1
            if current_counts[tax_out] == 0:
                del current_counts[tax_out]

        if current_counts:
            if len(current_counts) > 1 and "0" in current_counts:
                best_tax = "0"
                max_val = -1
                for tax, count in current_counts.items():
                    if tax != "0" and count > max_val:
                        max_val = count
                        best_tax = tax
                nuc_taxids.append(best_tax)
            else:
                nuc_taxids.append(max(current_counts, key=current_counts.get))
        else:
            nuc_taxids.append("0")

    return nuc_taxids


def find_longest_cds_optimized(seq_str: str) -> CdsOrf:
    best_orf: CdsOrf = {"strand": "+", "start": 0, "end": 0, "length": 0, "seq": ""}
    L = len(seq_str)

    for frame in range(3):
        rem = (L - frame) % 3
        sub_seq = seq_str[frame:-rem] if rem > 0 else seq_str[frame:]
        if not sub_seq:
            continue

        aa_seq = str(Seq(sub_seq).translate())
        orfs = aa_seq.split("*")
        current_pos = 0
        for orf in orfs:
            orf_len = len(orf)
            if orf_len > best_orf["length"]:
                start_nuc = frame + (current_pos * 3)
                best_orf = {
                    "strand": "+",
                    "start": start_nuc,
                    "end": start_nuc + (orf_len * 3),
                    "length": orf_len,
                    "seq": sub_seq[current_pos * 3 : current_pos * 3 + orf_len * 3],
                }
            current_pos += orf_len + 1

    rev_seq_str = str(Seq(seq_str).reverse_complement())
    for frame in range(3):
        rem = (L - frame) % 3
        sub_seq = rev_seq_str[frame:-rem] if rem > 0 else rev_seq_str[frame:]
        if not sub_seq:
            continue

        aa_seq = str(Seq(sub_seq).translate())
        orfs = aa_seq.split("*")
        current_pos = 0
        for orf in orfs:
            orf_len = len(orf)
            if orf_len > best_orf["length"]:
                start_nuc_rev = frame + (current_pos * 3)
                end_nuc_rev = start_nuc_rev + (orf_len * 3)
                best_orf = {
                    "strand": "-",
                    "start": L - end_nuc_rev,
                    "end": L - start_nuc_rev,
                    "length": orf_len,
                    "seq": sub_seq[current_pos * 3 : current_pos * 3 + orf_len * 3],
                }
            current_pos += orf_len + 1

    return best_orf


def compute_cai(cds_seq: str, weights: dict[str, float]) -> float:
    log_w_sum = 0.0
    codon_count = 0
    for i in range(0, len(cds_seq) - 2, 3):
        codon = cds_seq[i : i + 3]
        if codon in weights and weights[codon] > 0:
            log_w_sum += math.log(weights[codon])
            codon_count += 1
    return math.exp(log_w_sum / codon_count) if codon_count > 0 else 0.0


def analyze_codon_adaptation(
    reads_path: Union[str, Path],
    kraken_path: Union[str, Path],
    *,
    kraken_data: Optional[dict[str, KrakenReadInfo]] = None,
    codon_usage_store: Optional[CodonUsageStore] = None,
    codon_usage_dir: Optional[Union[str, Path]] = None,
    include_read_ids: Optional[set[str]] = None,
    kmer_len: int = 35,
) -> list[CodonAdaptationResult]:
    """
    Score codon adaptation (CAI) vs pre-built host reference tables.

    Requires codon_tables.json from build_codon_reference(); no network access.
    """
    reads_path = Path(reads_path)
    kraken_path = Path(kraken_path)
    if kraken_data is None:
        kraken_data = parse_kraken_file(kraken_path)

    if codon_usage_store is not None:
        store = codon_usage_store
    else:
        data_dir = Path(codon_usage_dir) if codon_usage_dir else default_codon_usage_dir()
        store = CodonUsageStore.load(data_dir)

    pending: list[PendingCodonRead] = []
    file_format = "fastq" if reads_path.suffix.lower() in (".fastq", ".fq") else "fasta"

    for record in SeqIO.parse(reads_path, file_format):
        read_id = record.id
        if include_read_ids is not None and read_id not in include_read_ids:
            continue
        if read_id not in kraken_data:
            continue

        status, taxid, r_len, k_info = kraken_data[read_id]
        if status == "U":
            continue

        seq_str = str(record.seq).upper()
        if not k_info:
            continue

        kmers = expand_kmer_taxids(k_info, r_len, k=kmer_len)
        nuc_taxids = map_nucleotides_to_taxids_fast(kmers, r_len, k=kmer_len)

        cds = find_longest_cds_optimized(seq_str)
        if cds["length"] == 0:
            continue

        host_tax_votes = Counter(nuc_taxids[: cds["start"]] + nuc_taxids[cds["end"] :])
        host_taxid = host_tax_votes.most_common(1)[0][0] if host_tax_votes else "0"

        pending.append((read_id, cds, host_taxid))

    host_taxids = {p[2] for p in pending if p[2] != "0"}
    store.require_host_taxids(host_taxids)

    results: list[CodonAdaptationResult] = []
    for read_id, cds, host_taxid in pending:
        ref_taxid: Optional[str] = None
        cai_score: Optional[float] = None

        if host_taxid != "0" and len(cds["seq"]) >= 3:
            weights, ref_taxid = store.get_cai_weights_for_host(host_taxid)
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
            )
        )

    return results


def codon_adaptation_to_tsv_lines(results: Iterable[CodonAdaptationResult]) -> list[str]:
    header = (
        "Read_ID\tCDS_Strand\tCDS_Range\tHost_TaxID\t"
        "Reference_TaxID\tCDS_Len_bp\tCAI_vs_Host"
    )
    lines = [header]
    for r in results:
        cai = f"{r.cai_vs_host:.4f}" if r.cai_vs_host is not None else "NA"
        lines.append(
            f"{r.read_id}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
            f"{r.host_taxid}\t{r.reference_taxid or 'NA'}\t{r.cds_len_bp}\t{cai}"
        )
    return lines


def write_codon_adaptation_results_tsv(
    output_path: Union[str, Path],
    results: Iterable[CodonAdaptationResult],
    *,
    cai_engineered_threshold: Optional[float] = None,
) -> str:
    """Write precomputed codon adaptation results to TSV; returns output path."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        if cai_engineered_threshold is None:
            for line in codon_adaptation_to_tsv_lines(results):
                f.write(line + "\n")
        else:
            header = (
                "Read_ID\tCDS_Strand\tCDS_Range\tHost_TaxID\t"
                "Reference_TaxID\tCDS_Len_bp\tCAI_vs_Host\tCodon_CAI_Threshold\tEngineered_By_Codon_CAI"
            )
            f.write(header + "\n")
            for r in results:
                cai = r.cai_vs_host
                cai_s = f"{cai:.4f}" if cai is not None else "NA"
                engineered = (
                    "NA" if cai is None else str(cai < cai_engineered_threshold)
                )
                f.write(
                    f"{r.read_id}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
                    f"{r.host_taxid}\t{r.reference_taxid or 'NA'}\t{r.cds_len_bp}\t"
                    f"{cai_s}\t{cai_engineered_threshold}\t{engineered}\n"
                )
    return str(out)


def write_codon_adaptation_tsv(
    reads_path: Union[str, Path],
    kraken_path: Union[str, Path],
    output_path: Union[str, Path],
    *,
    kmer_len: int = 35,
    include_read_ids: Optional[set[str]] = None,
    codon_usage_dir: Optional[Union[str, Path]] = None,
    codon_usage_store: Optional[CodonUsageStore] = None,
    cai_engineered_threshold: Optional[float] = None,
) -> tuple[str, list[CodonAdaptationResult]]:
    """Write Kraken-ORF codon adaptation TSV (legacy); returns (output path, results)."""
    results = analyze_codon_adaptation(
        reads_path,
        kraken_path,
        codon_usage_dir=codon_usage_dir,
        codon_usage_store=codon_usage_store,
        include_read_ids=include_read_ids,
        kmer_len=kmer_len,
    )
    out_path = write_codon_adaptation_results_tsv(
        output_path,
        results,
        cai_engineered_threshold=cai_engineered_threshold,
    )
    return out_path, results


def analyze_codon_adaptation_with_diamond(
    reads_path: Union[str, Path],
    *,
    diamond_db: Union[str, Path] | None = None,
    diamond_threads: int = 4,
    diamond_extra_args: Optional[list[str]] = None,
    run_diamond: bool = True,
    diamond_output_path: Union[str, Path] | None = None,
    debug_write_diamond_output: bool = False,
    taxonomy_parents: Optional[dict[str, str]] = None,
    codon_usage_store: Optional[CodonUsageStore] = None,
    codon_usage_dir: Optional[Union[str, Path]] = None,
    include_read_ids: Optional[set[str]] = None,
) -> tuple[list[CodonAdaptationResult], Path | None]:
    """
    DIAMOND-based ORF + host-taxid inference for CAI (SeqScreen-Nano style).

    - ORF intervals and CDS coordinates come from DIAMOND qstart/qend (--min-orf).
    - Host taxid is inferred from DIAMOND hit staxids only (majority/LCA over all hits).
    - CAI is computed vs the inferred host taxid (CSDB reference store).
    - Set ``debug_write_diamond_output`` + ``diamond_output_path`` to save TSV for reuse.
    - Set ``run_diamond=False`` + ``diamond_output_path`` to load a saved TSV.
    """
    reads_path = Path(reads_path)

    if codon_usage_store is not None:
        store = codon_usage_store
    else:
        data_dir = Path(codon_usage_dir) if codon_usage_dir else default_codon_usage_dir()
        store = CodonUsageStore.load(data_dir)

    parents = taxonomy_parents if taxonomy_parents is not None else store.taxonomy_parents()

    diamond_lines, diamond_path = resolve_diamond_lines(
        reads_path,
        diamond_db,
        run_diamond=run_diamond,
        output_path=diamond_output_path,
        debug_write_output=debug_write_diamond_output,
        threads=diamond_threads,
        extra_args=diamond_extra_args,
    )
    orfs_by_read, host_by_read = infer_orfs_and_host_taxids(
        diamond_lines, taxonomy_parents=parents
    )

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
        cds_seq = seq_str[best.start : best.end]
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
            weights, ref_taxid = store.get_cai_weights_for_host(host_taxid)
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


def main(reads_path: str, kraken_path: str, kmer_len: int = 35) -> list[CodonAdaptationResult]:
    results = analyze_codon_adaptation(reads_path, kraken_path, kmer_len=kmer_len)
    for line in codon_adaptation_to_tsv_lines(results):
        print(line)
    return results


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Usage: python -m plasmidScreen.src.analyze_codon_usage "
            "<reads.fasta/.fastq> <kraken2_output.txt> [kmer_length=35]"
        )
        sys.exit(1)

    k_len = int(sys.argv[3]) if len(sys.argv) > 3 else 35
    main(sys.argv[1], sys.argv[2], k_len)
