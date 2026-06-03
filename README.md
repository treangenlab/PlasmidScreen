# PlasmidScreen: An Accurate and Fast Engineered DNA Detection

### Description

PlasmidScreen detects engineered DNA in sequencing reads using a **Kraken2 minimizer scan** (taxid 32630 in sliding windows). Optional **codon adaptation index (CAI)** scoring uses **DIAMOND blastx** for ORF coordinates and host taxonomy, compared against **pre-built** host codon usage tables from the [Codon Statistics Database](http://codonstatsdb.unr.edu/) (airgapped-safe at runtime; no network during screening).

**Design specification (inputs/outputs):** [docs/DESIGN.md](docs/DESIGN.md)

### Pipeline overview

| Step | Tool | Purpose |
|------|------|---------|
| 1 (offline) | CSDB build | `codon_tables.json` + optional `taxonomy_parents.json` |
| 2 | Kraken2 | Engineered k-mer labels (Natural / Synthetic) |
| 3 | DIAMOND blastx | ORF intervals (`qstart`/`qend`) and host taxid (`staxids`) |
| 4 | CAI | Codon adaptation vs host reference (Natural reads only in full screen) |
| 5 | Combined call | `ScreenResult.per_read[].engineered_overall` / `overall_label` |

### Overall engineered decision (`run_screen` / `ScreenResult`)

Each read in `screen_result.per_read` gets a single **`engineered_overall`** flag and **`overall_label`** (`Natural` or `Synthetic`) from the thresholds you pass to `run_screen`:

| Signal | Parameter | Counts toward overall engineered when |
|--------|-----------|----------------------------------------|
| K-mer scan | `--threshold`, `--window_size` | Engineered minimizers (taxid 32630) in a window ≥ threshold → k-mer **Synthetic** |
| Codon CAI | `--codon-cai-threshold` (optional) | CAI &lt; threshold on a scored Natural read (requires codon step) |

**Rule:** `engineered_overall = engineered_by_kmer_scan OR engineered_by_codon_cai`.

- Without `--codon-cai-threshold`, overall matches the k-mer scan only (`engineered_scan.synthetic_count` == `overall_synthetic_count`).
- With `--codon-cai-threshold`, a read can be k-mer **Natural** but overall **Synthetic** if CAI is low.
- Use **`screen_result.overall_synthetic_count`**, **`engineered_read_ids`**, and **`natural_read_ids_overall`** for run-level summaries (not only `engineered_scan.synthetic_count`).

The written engineered report TSV (`report.txt`) still records **k-mer scan** labels only. The combined decision lives in the library `ScreenResult` / `ReadFlagDetail` objects.

### Library usage

```python
from plasmidScreen import (
    build_codon_reference,
    analyze_codon_adaptation,
    run_screen,
    taxids_from_kraken_output,
)

# Step 1 — networked machine: build reference JSON from Codon Statistics Database
build_codon_reference(
    "codon_usage/",
    csdb_archive="/path/to/codonstatsdb_March2022.tar.gz",  # optional; auto-download if omitted
)

# Step 2 — airgapped: codon-only analysis (DIAMOND + CSDB reference)
results, diamond_path = analyze_codon_adaptation(
    "reads.fa",
    diamond_db="/path/to/protein.dmnd",
    codon_usage_dir="codon_usage/",
    include_read_ids={"read1", "read2"},  # optional read filter
)
for row in results:
    print(row.read_id, row.host_taxid, row.reference_taxid, row.cai_vs_host)

# Full pipeline: Kraken engineered scan + codon CAI on Natural reads
screen_result = run_screen(
    "reads.fa",
    kraken_db="/path/to/kraken/db",
    diamond_db="/path/to/protein.dmnd",
    codon_usage_dir="codon_usage/",
    engineered_report_path="engineered_report.txt",  # omit for in-memory results only
    codon_cai_engineered_threshold=0.7,
)
# K-mer-only vs combined counts (differ when --codon-cai-threshold is set)
print("k-mer synthetic:", screen_result.engineered_scan.synthetic_count)
print("overall synthetic:", screen_result.overall_synthetic_count)
print(len(screen_result.codon_adaptation))
for r in screen_result.per_read:
    print(
        r.read_id,
        r.overall_label,
        r.engineered_overall,
        r.engineered_methods,
        r.cai_vs_host,
    )
print(screen_result.overall_synthetic_count, screen_result.engineered_read_ids)
```

### CLI — build reference (network required)

```bash
# Default: import all taxids in the CSDB archive (~5.2 GB download on first run)
python plasmidScreen.py build

# Use an existing CSDB archive (no download)
python plasmidScreen.py build \
  --csdb-archive /data/codonstatsdb_March2022.tar.gz \
  --no-download-csdb

# Subset: explicit taxids or taxids seen in a Kraken classifications file
python plasmidScreen.py build --taxids 9606,511145
python plasmidScreen.py build --taxids-file hosts.txt
```

The archive is cached under `~/.local/share/PlasmidScreen/` unless `--csdb-archive` is set.
Output: `codon_usage/codon_tables.json` and optionally `taxonomy_parents.json`.

### CLI — screen (airgapped after reference build)

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

### Tests

```bash
pip install -r requirements.txt
pytest
```

### Demo

```bash
PlasmidScreen example/full_length.fa
```
