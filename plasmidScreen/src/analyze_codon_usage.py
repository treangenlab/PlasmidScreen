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


def parse_kraken_file(kraken_path: str) -> dict[str, tuple[str, str, int, str]]:
    """Parse Kraken2 output into read_id -> (status, taxid, length, kmer_info)."""
    kraken_data: dict[str, tuple[str, str, int, str]] = {}
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


def find_longest_cds_optimized(seq_str: str) -> dict:
    best_orf = {"strand": "+", "start": 0, "end": 0, "length": 0, "seq": ""}
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
    kraken_data = parse_kraken_file(str(kraken_path))

    if codon_usage_store is not None:
        store = codon_usage_store
    else:
        data_dir = Path(codon_usage_dir) if codon_usage_dir else default_codon_usage_dir()
        store = CodonUsageStore.load(data_dir)

    pending: list[tuple] = []
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

        cds_tax_votes = Counter(nuc_taxids[cds["start"] : cds["end"]])
        host_tax_votes = Counter(nuc_taxids[: cds["start"]] + nuc_taxids[cds["end"] :])

        cds_taxid = cds_tax_votes.most_common(1)[0][0] if cds_tax_votes else taxid
        host_taxid = host_tax_votes.most_common(1)[0][0] if host_tax_votes else "0"

        pending.append(
            (read_id, taxid, cds, cds_taxid, host_taxid)
        )

    host_taxids = {p[4] for p in pending if p[4] != "0"}
    store.require_host_taxids(host_taxids)

    results: list[CodonAdaptationResult] = []
    for read_id, taxid, cds, cds_taxid, host_taxid in pending:
        ref_taxid: Optional[str] = None
        cai_score: Optional[float] = None

        if host_taxid != "0" and len(cds["seq"]) >= 3:
            weights, ref_taxid = store.get_cai_weights_for_host(host_taxid)
            if weights:
                cai_score = compute_cai(cds["seq"], weights)

        results.append(
            CodonAdaptationResult(
                read_id=read_id,
                overall_taxid=taxid,
                cds_strand=cds["strand"],
                cds_start=cds["start"],
                cds_end=cds["end"],
                cds_taxid=cds_taxid,
                host_taxid=host_taxid,
                reference_taxid=ref_taxid,
                cds_len_bp=len(cds["seq"]),
                cai_vs_host=cai_score,
            )
        )

    return results


def codon_adaptation_to_tsv_lines(results: Iterable[CodonAdaptationResult]) -> list[str]:
    header = (
        "Read_ID\tOverall_TaxID\tCDS_Strand\tCDS_Range\tCDS_TaxID\tHost_TaxID\t"
        "Reference_TaxID\tCDS_Len_bp\tCAI_vs_Host"
    )
    lines = [header]
    for r in results:
        cai = f"{r.cai_vs_host:.4f}" if r.cai_vs_host is not None else "NA"
        lines.append(
            f"{r.read_id}\t{r.overall_taxid}\t{r.cds_strand}\t{r.cds_start}-{r.cds_end}\t"
            f"{r.cds_taxid}\t{r.host_taxid}\t{r.reference_taxid or 'NA'}\t{r.cds_len_bp}\t{cai}"
        )
    return lines


def write_codon_adaptation_tsv(
    reads_path: Union[str, Path],
    kraken_path: Union[str, Path],
    output_path: Union[str, Path],
    *,
    kmer_len: int = 35,
    include_read_ids: Optional[set[str]] = None,
    codon_usage_dir: Optional[Union[str, Path]] = None,
    codon_usage_store: Optional[CodonUsageStore] = None,
) -> tuple[str, list[CodonAdaptationResult]]:
    """Write codon adaptation TSV; returns (output path, structured results)."""
    results = analyze_codon_adaptation(
        reads_path,
        kraken_path,
        codon_usage_dir=codon_usage_dir,
        codon_usage_store=codon_usage_store,
        include_read_ids=include_read_ids,
        kmer_len=kmer_len,
    )
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for line in codon_adaptation_to_tsv_lines(results):
            f.write(line + "\n")
    return str(out), results


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
