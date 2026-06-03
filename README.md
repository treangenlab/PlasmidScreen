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
print(screen_result.engineered_scan.synthetic_count)
print(len(screen_result.codon_adaptation))
for r in screen_result.per_read:
    print(r.read_id, r.engineered_methods, r.engineered_kmer_max_in_window, r.cai_vs_host)
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

# Engineered k-mer scan only (no codon / DIAMOND)
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
