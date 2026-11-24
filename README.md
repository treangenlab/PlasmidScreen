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
| 4 | CAI | Codon adaptation vs host reference (remaining Natural reads only in full screen)  |
| 5 | Combined call | `ScreenResult.per_read[].engineered_overall` / `overall_label`                    |
| 6 | blastn map (optional) | Top engineered-read reference hits on a nucleotide DB (`--map-hits`)              |
| 7 | Visualize (optional) | Pileup PNGs + mapping JSON (`--visualize` / `plasmidscreen visualize`)       |

### Overall engineered decision (`run_screen` / `ScreenResult`)

Each read in `screen_result.per_read` gets an **`overall_label`** (`Natural` or `Synthetic`) from the thresholds you 
pass to `run_screen` which can be a result from the two components K-mer scan, and Codon CAI which each are optionally 
ran. If either are found to contain engineering, `ScreenResult` object will report engineering.

| Signal | Parameter | Counts toward overall engineered when |
|--------|-----------|----------------------------------------|
| K-mer scan | `--threshold`, `--window_size`, `--kmer-size`, `--minimizer-size`, `--bg` / `--background-rate`, `--pvalue-alpha` | Window ≥ threshold **and** binomial p-value ≤ `--pvalue-alpha` (default 0.05) → k-mer engineered |
| Codon CAI | `--codon-cai-threshold` (optional) | CAI &gt; threshold on a scored Natural read (requires codon step) |

### Library usage

```python
from plasmidScreen import (
    build_codon_database,
    map_and_visualize_engineered,
    map_engineered_reads,
    run_screen,
    visualize_engineered_mappings,
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

# Optional: BLAST engineered reads against a nucleotide reference, then pileups
screen_result, viz = map_and_visualize_engineered(
    screen_result,
    query_fasta="reads.fa",
    reference_db="/path/to/plasmids.fa",
    viz_dir="viz_out",
    top_n=5,
)
# Or stepwise: map_engineered_reads(...) then visualize_engineered_mappings(...)
print(screen_result.top_engineered_references)
print(viz.pileup_paths)
assert screen_result.visualization is viz
for r in screen_result.per_read_engineered:
    print(
        r.read_id,
        r.matched_plasmid,
        r.match_target_start,
        r.match_target_end,
        r.match_sequence_identity,
        r.match_score,
        r.overlapping_features,
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
Optional map/visualize needs `blastn` on `PATH` (NCBI BLAST+).
Kraken2 installation process can be found on their wiki: https://github.com/DerrickWood/kraken2/wiki/Manual.
Similarly installing diamond can be found in their wiki: https://github.com/bbuchfink/diamond/wiki/2.-Installation,

### Recommended Databases

We recommend custom databases for both Kraken2 and DIAMOND. The suggested Kraken2 database contains a curated set of
sequences tailored for engineered-DNA signal. For codon optimization, host taxonomy of the ORF matters, so a
comprehensive DIAMOND database is recommended (e.g. UniRef100 with taxids). A BLAST database is also needed for the engineered signal alone for --map-hits option.

- K-mer scan: Custom Kraken2 database (~20GB)
- Codon optimization: UniRef100 DIAMOND database (~235GB)
- BLAST step: EngineeredSeqsBlastDB
### CLI — screen

Input reads may be plain or gzip-compressed FASTA/FASTQ
(`.fa`, `.fasta`, `.fq`, `.fastq`, and the same with `.gz`).
By default the CLI logs short progress milestones at INFO. Pass `--verbose` /
`-v` for DEBUG detail (commands, timings, per-file writes), or `--quiet` /
`-q` to hide progress bars and INFO messages (warnings and errors still show).

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
  --threshold 15

# Engineered k-mer scan only (no codon / DIAMOND; overall == k-mer labels)
plasmidscreen screen reads.fa report.txt /path/to/kraken/db --skip-codon-usage

# Optional: BLAST map engineered reads + write pileup visualizations
plasmidscreen screen reads.fa report.txt /path/to/kraken/db --skip-codon-usage \
  --visualize --reference-db /path/to/plasmids.fa --viz-dir viz_out \
  --viz-top-n 5

# Map only (no figures): attach top_engineered_references on ScreenResult
plasmidscreen screen reads.fa report.txt /path/to/kraken/db --skip-codon-usage \
  --map-hits --reference-db /path/to/plasmids.fa
```

#### Engineered report TSV (`report.txt`)

Written at the end of `screen`. Columns:

| Column | Meaning |
|--------|---------|
| `Label` | Overall call (`Engineered` / `Natural`) |
| `Read_ID` | Read identifier |
| `Methods` | Contributing methods (`engineered_kmer_scan`, `codon_cai`, …) or `NA` |
| `P-Value` | Binomial p-value for engineered k-mer count vs `--bg` |
| `Eng-Coverage Percentage` | Fraction of mapped k-mers that are engineered |
| `Plasmid` | Primary BLAST accession (best bitscore), or `NA` until map/viz |
| `Target_Start` / `Target_End` | Alignment coords on that plasmid, or `NA` |
| `Sequence_Identity` | Percent identity of the primary hit, or `NA` |
| `Match_Score` | Primary hit score (blastn bitscore), or `NA` |
| `Overlapping_Features` | GenBank feature labels overlapped by the primary hit (after map; needs GenBank cache), or `NA` |

Screen fills the first five columns immediately. `--map-hits` (or `--visualize`,
which maps first) rewrites the same path to populate BLAST columns and
overlapping features. Visualization only adds pileup PNGs under `--viz-dir`.

Example after `--visualize`:

```text
Label	Read_ID	Methods	P-Value	Eng-Coverage Percentage	Plasmid	Target_Start	Target_End	Sequence_Identity	Match_Score	Overlapping_Features
Natural	read_natural_1	NA	1.0	0.00%	NA	NA	NA	NA	NA	NA
Engineered	read_synthetic_1	engineered_kmer_scan	1.2e-10	14.75%	MK178573.1	5846	8105	99.115	3995	AmpR,ColE1 origin
```

### CLI — prefetch GenBank (networked, once)

Pileup feature tracks read a **local GenBank cache** only. Prefill it from the
same nucleotide reference FASTA you pass to `--reference-db` / `visualize`:

```bash
plasmidscreen prefetch-genbank /path/to/plasmids.fa
# optional: --cache-dir /path/to/cache --ncbi-api-key YOUR_KEY
```

Default cache: `~/.local/share/PlasmidScreen/genbank_cache`. Visualization does
not contact NCBI by default; missing accessions simply omit feature overlays.
For a quick install without a bulk prefetch, pass `--fetch-genbank` on
`screen --visualize` / `visualize` to download only the accessions needed for
that run (optional `--ncbi-api-key`).

### CLI — visualize (from a prior screen report)

Re-run blastn enrichment and pileups without re-screening. Rebuilds a minimal
`ScreenResult` from the engineered report TSV, then maps and visualizes.

```bash
plasmidscreen visualize reads.fa report.txt /path/to/plasmids.fa viz_out \
  --top-n 5 --threads 8
```

Outputs under `--viz-dir` / `viz_out`:
- `*.png` pileup figures for top references

Also remaps and rewrites the engineered report TSV (same path as the prior
`screen` run) with BLAST / feature columns, then writes pileups. Those fields
are available on `ScreenResult.per_read` (`matched_plasmid`, `match_target_*`,
`match_score`, `overlapping_features`).

Library equivalents: `prefetch_genbank_from_fasta`, `map_engineered_reads`,
`visualize_engineered_mappings`, `map_and_visualize_engineered`, and
`screen_result_from_engineered_report`.


#### P-Value Documentation
The engineered report includes engineered k-mer coverage and a p-value for each read.
The p-value estimates the probability of observing at least this many engineered
k-mers by chance under a binomial null with background rate `--bg` (default 0.1%).
A window-threshold Synthetic call is kept as engineered only when
`p ≤ 0.05`; otherwise the k-mer component is reassigned Natural (CAI may still
flag the read). Lower `--bg` makes significance easier to reach.
This rate was established on a curated dataset of non-engineered sequences,
comprising bacterial chromosomes (>22,000) and naturally occurring incomplete plasmids from IMG-PR (>561,000),
each chopped into 1000 bp fragments.

The engineered call gates on the raw p-value threshold above (BH FDR helpers
remain available in the library but are not used on the live screen path).

