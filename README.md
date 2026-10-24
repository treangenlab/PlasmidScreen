# PlasmidScreen: An Accurate and Fast Engineered DNA Detection

### Description

PlasmidScreen detects engineered DNA in sequencing reads using a Kraken2 minimizer/kmer scan 
(taxid 32630 in sliding windows) with a custom database buiilt for engineering detection (preferably plasmids).
Optionally, PlasmidScreen additionally searches for signal of codon optimization of coding regions provided through 
codon adaptation index (CAI) by computing the nearest neighbor in UniRef powered by*DIAMOND blastx. The identified taxonomy is 
compared against the codon usage tables from the [Codon Statistics Database](http://codonstatsdb.unr.edu/). Computed CAI above a 
threshold is identified as engineered.

**Design specification (inputs/outputs):** [docs/DESIGN.md](docs/DESIGN.md)

### Pipeline overview

| Step | Tool | Purpose                                                                           |
|------|------|-----------------------------------------------------------------------------------|
| 1 | CSDB build | `codon_tables.json` + optional `taxonomy_parents.json`                            |
| 2 | Kraken2 | Engineered k-mer labels (Natural / Synthetic)                                     |
| 3 | DIAMOND blastx | ORF intervals (`qstart`/`qend`) and host taxid (`staxids`)                        |
| 4 | CAI | Codon adaptation vs host reference (K-mer scan natural reads only in full screen) |
| 5 | Combined call | `ScreenResult.per_read[].engineered_overall` / `overall_label`                    |

### Overall engineered decision (`run_screen` / `ScreenResult`)

Each read in `screen_result.per_read` gets an **`overall_label`** (`Natural` or `Synthetic`) from the thresholds you 
pass to `run_screen` which can be a result from the two components K-mer scan, and Codon CAI which each are optionally 
ran. If either are found to contain engineering, `ScreenResult` object will report engineering.

| Signal | Parameter | Counts toward overall engineered when |
|--------|-----------|----------------------------------------|
| K-mer scan | `--threshold`, `--window_size` | Engineered minimizers (taxid 32630) in a window ≥ threshold → k-mer **Synthetic** |
| Codon CAI | `--codon-cai-threshold` (optional) | CAI &gt; threshold on a scored Natural read (requires codon step) |

### Library usage

```python
from plasmidScreen import (
    build_codon_database,
    run_screen,
)

# Step 1 — build reference JSON from Codon Statistics Database.
# If output_dir is omitted, data is written under ~/.local/share/PlasmidScreen/
build_codon_database(output_dir="/path/to/codon_usage_db")

# Full pipeline: Kraken engineered scan + codon CAI on Natural reads
# seq_type: 0 = long read, 1 = short read (default)
screen_result = run_screen(
    "reads.fa",
    1,
    "/path/to/kraken/db",
    diamond_db="/path/to/protein.dmnd",
    codon_usage_dir="/path/to/codon_usage_db",
    engineered_report_path="engineered_report.txt",  # omit for in-memory only
    codon_cai_engineered_threshold=0.6,
)

print(screen_result.engineered_read_ids)                   # overall (k-mer ∪ CAI)
print(screen_result.engineered_scan.engineered_read_ids)   # k-mer Synthetic
print(screen_result.codon_cai_engineered_read_ids)         # CAI above threshold
print(screen_result.natural_read_ids_overall)              # overall Natural

for r in screen_result.per_read:
    print(
        r.read_id,
        r.overall_label,            # "Natural" or "Synthetic"
        r.engineered_overall,       # True if overall engineered
        r.engineered_methods,       # e.g. ["engineered_kmer_scan"] or ["codon_cai"]
        r.engineered_by_kmer_scan,
        r.engineered_by_codon_cai,  # True/False/None vs CAI threshold
        r.cai_vs_host,              # float CAI score, or None
    )
```

### CLI — build reference (network required)

```bash
# Default: import all taxids in the CSDB archive (~5.2 GB download on first run)
plasmidscreen build
# (or: python -m plasmidScreen build)

# Use an existing CSDB archive (no download)
plasmidscreen build \
  --csdb-archive /data/codonstatsdb_March2022.tar.gz \
  --no-download-csdb

```
The archive is cached under `~/.local/share/PlasmidScreen/` unless `--csdb-archive` is set.
Output: `codon_usage/codon_tables.json` and `taxonomy_parents.json`.


### Installation

```bash
# Installs package + console script `plasmidscreen`
pip install .
```

Requires Python ≥3.10. Also need `kraken2` and (for codon CAI) `diamond` on `PATH`.
Kraken2 installation process can be found on their wiki: https://github.com/DerrickWood/kraken2/wiki/Manual.
Similarly installing diamond can be found in their wiki: https://github.com/bbuchfink/diamond/wiki/2.-Installation,

### Recommended Databases

We recommend custom databases for both Kraken2 and DIAMOND. The suggested Kraken2 database contains a curated set of
sequences tailored for engineered-DNA signal. For codon optimization, host taxonomy of the ORF matters, so a
comprehensive DIAMOND database is recommended (e.g. UniRef100 with taxids).

- K-mer scan: Custom Kraken2 database (~20GB)
- Codon optimization: UniRef100 DIAMOND database (~235GB)

### CLI — screen

Input reads may be plain or gzip-compressed FASTA/FASTQ
(`.fa`, `.fasta`, `.fq`, `.fastq`, and the same with `.gz`).
By default the CLI logs progress at INFO; pass `--quiet` / `-q` to hide
progress bars and INFO messages (warnings and errors still show).

Use `--seq-type 0` for long reads (Nanopore/PacBio-style) or `--seq-type 1`
for short reads (default). Use `--mem-mode 0|1|2` to match available RAM for DIAMOND
(32 / 64 / 128 GB bins).

```bash
plasmidscreen screen reads.fa report.txt /path/to/kraken/db \
  --diamond-db /path/to/protein.dmnd \
  --codon-usage-dir ~/.local/share/PlasmidScreen/codon_usage \
  --codon-usage-output codon_usage.tsv \
  --seq-type 1

# Reuse Kraken classifications (no Kraken subprocess)
plasmidscreen screen reads.fa report.txt /path/to/kraken/db \
  --no-run-kraken --kraken-output-path kraken.out \
  --diamond-db /path/to/protein.dmnd

# Save DIAMOND TSV for debugging / codon-only reruns
plasmidscreen screen reads.fa report.txt /path/to/kraken/db \
  --diamond-db /path/to/protein.dmnd \
  --debug-write-diamond-out --diamond-output-path diamond.tsv

# Re-run codon step from saved DIAMOND output
plasmidscreen screen reads.fa report.txt /path/to/kraken/db \
  --no-run-diamond --diamond-output-path diamond.tsv \
  --no-run-kraken --kraken-output-path kraken.out

# CAI above threshold can flag overall engineered (on Natural reads with CAI scores)
plasmidscreen screen reads.fa report.txt /path/to/kraken/db \
  --diamond-db /path/to/protein.dmnd \
  --codon-cai-threshold 0.6 \
  --threshold 10

# Engineered k-mer scan only (no codon / DIAMOND; overall == k-mer labels)
plasmidscreen screen reads.fa report.txt /path/to/kraken/db --skip-codon-usage
```


#### P-Value Documentation
The output contains two columns representing the confidence that a given read is of engineered origin. 
The first column reports the percentage of k-mers classified as engineered out of all k-mers in the read. 
The second column provides a p-value estimating the probability of observing this fraction of engineered k-mers by chance.
The p-value is calculated using a binomial distribution based on the total number of classified k-mers,
the number of classified engineered k-mers, and a background rate of engineered k-mers in non-engineered sequences. 
This background rate represents the probability that a random, naturally occurring k-mer is misclassified as engineered.
We recommend a conservative default background rate of 0.1%; lower values result in more aggressive filtering. 
This rate was established on a curated dataset of non-engineered sequences, 
comprising bacterial chromosomes (>22,000) and naturally occurring incomplete plasmids from IMG-PR (>561,000), 
each chopped into 1000 bp fragments.

Reported p-values are adjusted for multiple testing across all classified reads using the Benjamini-Hochberg False 
Discovery Rate (FDR) correction. 
