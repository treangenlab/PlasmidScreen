# PlasmidScreen: An Accurate and Fast Engineered DNA Detection

### Description

PlasmidScreen detects engineered DNA in sequencing reads using a **Kraken2 minimizer/kmer scan** (taxid 32630 in sliding windows). Optional **codon adaptation index (CAI)** scoring uses **DIAMOND blastx** for ORF coordinates and host taxonomy, compared against host codon usage tables from the [Codon Statistics Database](http://codonstatsdb.unr.edu/).

**Design specification (inputs/outputs):** [docs/DESIGN.md](docs/DESIGN.md)

### Pipeline overview

| Step | Tool | Purpose |
|------|------|---------|
| 1 | CSDB build | `codon_tables.json` + optional `taxonomy_parents.json` |
| 2 | Kraken2 | Engineered k-mer labels (Natural / Synthetic) |
| 3 | DIAMOND blastx | ORF intervals (`qstart`/`qend`) and host taxid (`staxids`) |
| 4 | CAI | Codon adaptation vs host reference (Natural reads only in full screen) |
| 5 | Combined call | `ScreenResult.per_read[].engineered_overall` / `overall_label` |

### Overall engineered decision (`run_screen` / `ScreenResult`)

Each read in `screen_result.per_read` gets a single **`engineered_overall`** flag and **`overall_label`** (`Natural` or `Synthetic`) from the thresholds you pass to `run_screen` which can be a result from the two components K-mer scan, and Codon CAI which each are optionally ran. If either are found to contain engineering, `ScreenResult` object will report engineering.

| Signal | Parameter | Counts toward overall engineered when |
|--------|-----------|----------------------------------------|
| K-mer scan | `--threshold`, `--window_size` | Engineered minimizers (taxid 32630) in a window ≥ threshold → k-mer **Synthetic** |
| Codon CAI | `--codon-cai-threshold` (optional) | CAI &lt; threshold on a scored Natural read (requires codon step) |


- Without `--codon-cai-threshold`, overall matches the k-mer scan only (`engineered_scan.synthetic_count` == `overall_synthetic_count`).
- With `--codon-cai-threshold`, a read can be k-mer **Natural** but overall **Synthetic** if CAI is low.
- Use **`screen_result.overall_synthetic_count`**, **`engineered_read_ids`**, and **`natural_read_ids_overall`** for run-level summaries (not only `engineered_scan.synthetic_count`).

The written engineered report TSV (`report.txt`) still records **k-mer scan** labels only. The combined decision lives in the library `ScreenResult` / `ReadFlagDetail` objects.

### Library usage

```python
from plasmidScreen import (
    build_codon_database,
    run_screen,
)

# Step 1 — build reference JSON from Codon Statistics Database, 
# please specify path otherwise download will occur in ~/.local/share/PlasmidScreen/
build_codon_database(output_dir="/path/to/codon_usage_db")


# Full pipeline: Kraken engineered scan + codon optimization detection on Natural reads
screen_result = run_screen(
    "reads.fa",
    kraken_db="/path/to/kraken/db",
    diamond_db="/path/to/protein.dmnd",
    codon_usage_dir="/path/to/codon_usage_db",
    engineered_report_path="engineered_report.txt",  # omit for in-memory results only
    codon_cai_engineered_threshold=0.7,
)
# Here is the found engineered reads
print(screen_result.engineered_read_ids)
# Here is the engineered reads based on kmer scanning
print(screen_result.engineered_scan.engineered_read_ids)
# Here is the engineered reads from codon optimization detection
print(screen_result.codon_adaptation.engineered_read_ids)
# Here is the natural read ids found:
print(screen_result.natural_read_ids_overall)
# If you want to get more granular information of each read's call you can look at the list of reads
# This contains the ReadFlagDetail object that holds additonal info

for r in screen_result.per_read:
    print(
        r.read_id,
        r.overall_label, # Returns literal of natural or synthetic
        r.engineered_overall, # boolean True if engineered
        r.engineered_methods, # returns list of strings of type of engineering detection performed on the read
        r.engineered_by_kmer_scan, # boolean if engineered by kmer scan
        r.cai_vs_host,  # boolean if engineered by codon optimization threshold
    )
```

### CLI — build reference (network required)

```bash
# Default: import all taxids in the CSDB archive (~5.2 GB download on first run)
python plasmidScreen.py build

# Use an existing CSDB archive (no download)
python plasmidScreen.py build \
  --csdb-archive /data/codonstatsdb_March2022.tar.gz \
  --no-download-csdb

```

The archive is cached under `~/.local/share/PlasmidScreen/` unless `--csdb-archive` is set.
Output: `codon_usage/codon_tables.json` and optionally `taxonomy_parents.json`.

### CLI — screen

```bash
python plasmidScreen.py screen reads.fa report.txt /path/to/kraken/db \
  --diamond-db /path/to/protein.dmnd \
  --codon-usage-dir ~/.local/share/PlasmidScreen/codon_usage \
  --codon-usage-output codon_usage.tsv

# Reuse Kraken classifications (no Kraken subprocess)
python plasmidScreen.py screen reads.fa report.txt /path/to/kraken/db \
  --no-run-kraken --kraken-output-path kraken.out \
  --diamond-db /path/to/protein.dmnd

# Save DIAMOND TSV for debugging / codon-only reruns
python plasmidScreen.py screen reads.fa report.txt /path/to/kraken/db \
  --diamond-db /path/to/protein.dmnd \
  --debug-write-diamond-out --diamond-output-path diamond.tsv

# Re-run codon step from saved DIAMOND output
python plasmidScreen.py screen reads.fa report.txt /path/to/kraken/db \
  --no-run-diamond --diamond-output-path diamond.tsv \
  --no-run-kraken --kraken-output-path kraken.out

# Low CAI can flag overall engineered (on Natural reads with CAI scores)
python plasmidScreen.py screen reads.fa report.txt /path/to/kraken/db \
  --diamond-db /path/to/protein.dmnd \
  --codon-cai-threshold 0.7 \
  --threshold 25

# Engineered k-mer scan only (no codon / DIAMOND; overall == k-mer labels)
python plasmidScreen.py screen reads.fa report.txt /path/to/kraken/db --skip-codon-usage
```

If a host taxid has no resolvable codon table, screening raises `MissingCodonReferenceError` — build or extend the reference first (no online fetch at runtime).
