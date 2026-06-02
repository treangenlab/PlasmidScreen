# PlasmidScreen: An Accurate and Fast Engineered DNA Detection

### Description

PlasmidScreen is an engineered DNA detector capable of detecting inserts and small edits of reads and assemblies. Codon adaptation (CAI) scoring uses **pre-built** host codon usage tables (airgapped-safe; no runtime network).

**Design specification (inputs/outputs):** [docs/DESIGN.md](docs/DESIGN.md)

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

# Step 2 — airgapped: codon analysis (DIAMOND blastx + CSDB reference)
results, _diamond_path = analyze_codon_adaptation(
    "reads.fa",
    diamond_db="/path/to/protein.dmnd",
    codon_usage_dir="codon_usage/",
    include_read_ids={"read1", "read2"},  # optional read filter
)
for row in results:
    print(row.read_id, row.host_taxid, row.reference_taxid, row.cai_vs_host)

# Full pipeline (Kraken + engineered scan + codon usage; in-memory by default)
screen_result = run_screen(
    "reads.fa",
    kraken_db="/path/to/kraken/db",
    diamond_db="/path/to/protein.dmnd",
    codon_usage_dir="codon_usage/",
    engineered_report_path="engineered_report.txt",  # omit for library-only in-memory results
    codon_cai_engineered_threshold=0.7,
)
print(screen_result.engineered_scan.synthetic_count)
print(len(screen_result.codon_adaptation))
for r in screen_result.per_read:
    print(r.read_id, r.engineered_methods, r.engineered_kmer_max_in_window, r.cai_vs_host)
```

### CLI — build reference (network required)

```bash
# Default: all taxids in CSDB (downloads ~5.2 GB archive on first run)
python plasmidScreen.py build-codon-db build

# Use an existing CSDB archive (no download)
python plasmidScreen.py build-codon-db build --csdb-archive /data/codonstatsdb_March2022.tar.gz --no-download-csdb

# Or specify taxids / Kraken output explicitly
python plasmidScreen.py build-codon-db build --taxids 9606,511145
python plasmidScreen.py build-codon-db build --kraken-output kraken.out
```

Downloads the [Codon Statistics Database](http://codonstatsdb.unr.edu/) tar to
`~/.local/share/PlasmidScreen/codonstatsdb_March2022.tar.gz` unless `--csdb-archive` /
`--no-download-csdb` is set. Writes `codon_usage/codon_tables.json` and optionally
`taxonomy_parents.json`.

### CLI — screen (airgapped)

```bash
python plasmidScreen.py screen reads.fa report.txt \
  --diamond-db /path/to/protein.dmnd \
  --codon-usage-dir ~/.local/share/PlasmidScreen/codon_usage \
  --codon-usage-output codon_usage.tsv

# Save DIAMOND TSV for debugging / reuse (skip re-running DIAMOND later)
python plasmidScreen.py screen reads.fa report.txt \
  --diamond-db /path/to/protein.dmnd \
  --debug-write-diamond-out --diamond-output-path diamond.tsv

# Re-run codon step only from saved DIAMOND output
python plasmidScreen.py screen reads.fa report.txt \
  --no-run-diamond --diamond-output-path diamond.tsv
```

If a host taxid is missing from the reference, screening raises `MissingCodonReferenceError` (build the reference first; no online fetch at runtime).

### Tests

```bash
pip install -r requirements.txt
pytest
```

### Demo

```bash
PlasmidScreen example/full_length.fa
```
