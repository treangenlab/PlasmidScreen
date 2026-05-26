import sys
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from typing import Iterable, Optional

# 61 Sense Codons (Excluding typical translation termination stops: TAA, TAG, TGA)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'CAT': 'H', 'CAC': 'H', 
    'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 
    'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 
    'GGG': 'G'
}

AA_TO_CODONS = defaultdict(list)
for codon, aa in GENETIC_CODE.items():
    AA_TO_CODONS[aa].append(codon)

def parse_kraken_file(kraken_path: str):
    """Parse Kraken2 output.

    Expected format is tab-delimited per record. With minimizer data enabled, the
    k-mer runlength tokens are in the last column. This function stores only the
    minimal fields required for downstream codon usage analysis.
    """
    kraken_data = {}
    with open(kraken_path, 'r', buffering=1024*1024) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            # Kraken2 standard: status, read_id, taxid, length, [optional minimizer/kmer info...]
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

def expand_kmer_taxids(kmer_info, read_length, k=35):
    """Decompresses Kraken2 runlength tokens rapidly into an array."""
    kmers = []
    for block in kmer_info.split(' '):
        if not block:
            continue
        taxid, count = block.split(':')
        kmers.extend([taxid] * int(count))
    
    expected_kmers = read_length - k + 1
    if len(kmers) < expected_kmers:
        kmers.extend(['0'] * (expected_kmers - len(kmers)))
    return kmers[:expected_kmers]

def map_nucleotides_to_taxids_fast(kmers, read_length, k=35):
    """Computes exact consensus taxonomy per nucleotide in strict O(N) linear time."""
    nuc_taxids = []
    num_kmers = len(kmers)
    current_counts = {}
    
    for i in range(read_length):
        # 1. Add entering k-mer to the window
        if i < num_kmers:
            tax_in = kmers[i]
            current_counts[tax_in] = current_counts.get(tax_in, 0) + 1
            
        # 2. Remove leaving k-mer from the window
        leaving_idx = i - k
        if leaving_idx >= 0 and leaving_idx < num_kmers:
            tax_out = kmers[leaving_idx]
            current_counts[tax_out] -= 1
            if current_counts[tax_out] == 0:
                del current_counts[tax_out]
        
        # 3. Compute the dominant tax ID with strict priority matching
        if current_counts:
            if len(current_counts) > 1 and '0' in current_counts:
                # Find the highest vote that isn't '0' (unclassified)
                best_tax = '0'
                max_val = -1
                for tax, count in current_counts.items():
                    if tax != '0' and count > max_val:
                        max_val = count
                        best_tax = tax
                nuc_taxids.append(best_tax)
            else:
                # Optimized standard max retrieval
                nuc_taxids.append(max(current_counts, key=current_counts.get))
        else:
            nuc_taxids.append('0')
            
    return nuc_taxids

def find_longest_cds_optimized(seq_str):
    """Finds the longest unbroken 6-frame ORF using low-overhead string manipulations."""
    best_orf = {'strand': '+', 'start': 0, 'end': 0, 'length': 0, 'seq': ''}
    L = len(seq_str)
    
    # Forward frames
    for frame in range(3):
        rem = (L - frame) % 3
        sub_seq = seq_str[frame:-rem] if rem > 0 else seq_str[frame:]
        if not sub_seq: continue
        
        aa_seq = str(Seq(sub_seq).translate())
        orfs = aa_seq.split('*')
        
        current_pos = 0
        for orf in orfs:
            orf_len = len(orf)
            if orf_len > best_orf['length']:
                start_nuc = frame + (current_pos * 3)
                best_orf = {
                    'strand': '+', 'start': start_nuc, 'end': start_nuc + (orf_len * 3),
                    'length': orf_len, 'seq': sub_seq[current_pos*3 : current_pos*3 + orf_len*3]
                }
            current_pos += orf_len + 1

    # Reverse complement frames
    rev_seq_str = str(Seq(seq_str).reverse_complement())
    for frame in range(3):
        rem = (L - frame) % 3
        sub_seq = rev_seq_str[frame:-rem] if rem > 0 else rev_seq_str[frame:]
        if not sub_seq: continue
        
        aa_seq = str(Seq(sub_seq).translate())
        orfs = aa_seq.split('*')
        
        current_pos = 0
        for orf in orfs:
            orf_len = len(orf)
            if orf_len > best_orf['length']:
                start_nuc_rev = frame + (current_pos * 3)
                end_nuc_rev = start_nuc_rev + (orf_len * 3)
                best_orf = {
                    'strand': '-', 'start': L - end_nuc_rev, 'end': L - start_nuc_rev,
                    'length': orf_len, 'seq': sub_seq[current_pos*3 : current_pos*3 + orf_len*3]
                }
            current_pos += orf_len + 1
            
    return best_orf


def codon_adaptation_rows(
    reads_path: str,
    kraken_path: str,
    kmer_len: int = 35,
    include_read_ids: Optional[set[str]] = None,
) -> Iterable[str]:
    """Generate TSV rows for codon adaptation vs host from reads + Kraken output."""
    kraken_data = parse_kraken_file(kraken_path)

    host_profiles = defaultdict(Counter)
    processed_reads = []

    file_format = "fastq" if reads_path.endswith((".fastq", ".fq")) else "fasta"

    for record in SeqIO.parse(reads_path, file_format):
        read_id = record.id
        if include_read_ids is not None and read_id not in include_read_ids:
            continue
        if read_id not in kraken_data:
            continue

        status, taxid, r_len, k_info = kraken_data[read_id]
        if status == 'U':
            continue

        seq_str = str(record.seq).upper()
        if not k_info:
            continue

        kmers = expand_kmer_taxids(k_info, r_len, k=kmer_len)
        nuc_taxids = map_nucleotides_to_taxids_fast(kmers, r_len, k=kmer_len)

        cds = find_longest_cds_optimized(seq_str)
        if cds['length'] == 0:
            continue

        cds_tax_votes = Counter(nuc_taxids[cds['start']:cds['end']])
        host_tax_votes = Counter(nuc_taxids[:cds['start']] + nuc_taxids[cds['end']:])

        cds_taxid = cds_tax_votes.most_common(1)[0][0] if cds_tax_votes else taxid
        host_taxid = host_tax_votes.most_common(1)[0][0] if host_tax_votes else '0'

        if host_taxid != '0':
            h_prof = host_profiles[host_taxid]
            for i in range(0, cds['start'] - 2, 3):
                h_prof[seq_str[i:i+3]] += 1
            for i in range(cds['end'], r_len - 2, 3):
                h_prof[seq_str[i:i+3]] += 1

        processed_reads.append(
            (read_id, taxid, cds['strand'], cds['start'], cds['end'], cds_taxid, host_taxid, cds['seq'])
        )

    cai_weights = {}
    for taxid, counts in host_profiles.items():
        smoothed = {codon: counts.get(codon, 0) + 1 for codon in GENETIC_CODE}
        cai_weights[taxid] = {}
        for aa, codons in AA_TO_CODONS.items():
            max_c = max(smoothed[c] for c in codons)
            for c in codons:
                cai_weights[taxid][c] = smoothed[c] / max_c

    import math

    header = "Read_ID\tOverall_TaxID\tCDS_Strand\tCDS_Range\tCDS_TaxID\tHost_TaxID\tCDS_Len_bp\tAdaptive_CAI_vs_Host"
    yield header

    for r_id, ov_tax, strand, start, end, c_tax, h_tax, cds_seq in processed_reads:
        if h_tax in cai_weights and len(cds_seq) >= 3:
            w = cai_weights[h_tax]
            log_w_sum = 0.0
            codon_count = 0
            for i in range(0, len(cds_seq) - 2, 3):
                codon = cds_seq[i:i+3]
                if codon in w:
                    log_w_sum += math.log(w[codon])
                    codon_count += 1
            cai_score = math.exp(log_w_sum / codon_count) if codon_count > 0 else 0.0
        else:
            cai_score = 1.0

        yield f"{r_id}\t{ov_tax}\t{strand}\t{start}-{end}\t{c_tax}\t{h_tax}\t{len(cds_seq)}\t{cai_score:.4f}"


def write_codon_adaptation_tsv(
    reads_path: str,
    kraken_path: str,
    output_path: str,
    kmer_len: int = 35,
    include_read_ids: Optional[set[str]] = None,
) -> str:
    """Write codon adaptation TSV and return output path."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for row in codon_adaptation_rows(
            reads_path, kraken_path, kmer_len=kmer_len, include_read_ids=include_read_ids
        ):
            f.write(row + "\n")
    return str(out)

def main(reads_path, kraken_path, kmer_len=35):
    for row in codon_adaptation_rows(reads_path, kraken_path, kmer_len=kmer_len):
        print(row)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analyze_codon_adaptation_fast.py <reads.fasta/.fastq> <kraken2_output.txt> [kmer_length=35]")
        sys.exit(1)
    
    k_len = int(sys.argv[3]) if len(sys.argv) > 3 else 35
    main(sys.argv[1], sys.argv[2], k_len)
