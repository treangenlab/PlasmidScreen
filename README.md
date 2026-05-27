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
    taxids=["9606", "511145"],
    csdb_archive="/path/to/codonstatsdb_March2022.tar.gz",  # optional; auto-download if omitted
)

# Step 2 — airgapped: codon analysis only (structured results)
results = analyze_codon_adaptation(
    "reads.fa",
    "kraken.out",
    codon_usage_dir="codon_usage/",
    include_read_ids={"read1", "read2"},  # optional Natural-read filter
)
for row in results:
    print(row.read_id, row.host_taxid, row.reference_taxid, row.cai_vs_host)

# Full pipeline (Kraken + engineered scan + codon usage)
screen_result = run_screen(
    "reads.fa",
    "engineered_report.txt",
    "kraken.out",
    kraken_db="/path/to/kraken/db",
    codon_usage_dir="codon_usage/",
    # Optional: flag engineered by CAI threshold (CAI < threshold)
    codon_cai_engineered_threshold=0.7,
)
print(screen_result.engineered_scan.synthetic_count)
print(len(screen_result.codon_adaptation))
for r in screen_result.per_read:
    print(r.read_id, r.engineered_methods, r.engineered_kmer_max_in_window, r.cai_vs_host)
```

### CLI — build reference (network required)

```bash
# Default: ~150 curated taxids from CSDB (downloads ~5.2 GB archive on first run)
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
python plasmidScreen.py screen reads.fa report.txt kraken.out \
  --codon-usage-dir ~/.local/share/PlasmidScreen/codon_usage \
  --codon-usage-output codon_usage.tsv
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
