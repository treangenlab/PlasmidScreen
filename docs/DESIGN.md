# PlasmidScreen Design Specification

This document defines expected inputs, outputs, and data contracts for PlasmidScreen screening, codon adaptation analysis, and the offline codon reference build.

## Overview

```mermaid
flowchart LR
  subgraph prep [Offline prep - networked]
    A[Taxids / Kraken taxids] --> B[build_codon_reference]
    B --> C[codon_tables.json]
    B --> D[taxonomy_parents.json]
  end
  subgraph runtime [Runtime - airgapped]
    E[Reads FASTA/FASTQ] --> F[Kraken2]
    F --> G[kraken.out]
    G --> H[Engineered k-mer scan]
    H --> I[engineered_report.txt]
    G --> J[Codon adaptation]
    C --> J
    J --> K[codon_usage.tsv]
  end
```

| Phase | Network | Entry point |
|-------|---------|-------------|
| Reference build | Required | `build_codon_reference()` / `build-codon-db build` |
| Screening | Not required | `run_screen()` / `screen` |
| Codon-only | Not required | `analyze_codon_adaptation()` |

---

## 1. Inputs

### 1.1 Read sequences (FASTA / FASTQ)

| Property | Requirement |
|----------|-------------|
| Formats | `.fa`, `.fasta`, `.fq`, `.fastq` (suffix-based detection) |
| Encoding | Nucleotides A/C/G/T (case-insensitive; uppercased internally) |
| Read IDs | Must match Kraken2 output column 2 (first whitespace-delimited token in FASTA header) |

**Example FASTA**

```text
>read_natural_1
ATGAAATTTAAATAG
>read_synthetic_1
ATGAAATTTAAATAG
```

### 1.2 Kraken2 raw output

Produced with minimizer data enabled (required for engineered scan and codon analysis):

```bash
kraken2 --db <DB> --report-minimizer-data --output <kraken.out> \
  --use-names <reads.fa> --threads <N>
```

| Column (0-based) | Field | Description |
|------------------|-------|-------------|
| 0 | Status | `C` classified, `U` unclassified |
| 1 | Sequence ID | Read identifier (matches FASTA) |
| 2 | TaxID | NCBI taxonomy ID (`0` if unclassified) |
| 3 | Length | Sequence length (bp) |
| 4…n-2 | Optional | Additional Kraken fields |
| n-1 (last) | K-mer blocks | Run-length taxid assignments (see below) |

**K-mer block format (last column)**

Space-separated tokens: `TAXID:COUNT`

- Each token assigns `COUNT` consecutive k-mers (default k=35) the given taxonomy.
- `A:COUNT` marks ambiguous minimizers (handled as taxid `-1` internally).
- Engineered detection targets taxid **32630** in a sliding window.

**Example line**

```text
C	read_natural_1	562	100	0:35	562:30
```

### 1.3 Kraken2 database (engineered screen)

| Property | Description |
|----------|-------------|
| Path | Directory passed as `kraken_db` / positional argument |
| Usage | Local Kraken2 index; must include engineered taxon minimizers (taxid 32630) |

### 1.4 Codon usage reference (pre-built, airgapped)

Directory (default: `~/.local/share/PlasmidScreen/codon_usage/`):

| File | Required | Description |
|------|----------|-------------|
| `codon_tables.json` | Yes | Per-taxid codon relative frequencies |
| `taxonomy_parents.json` | No | NCBI `taxid → parent_taxid` for lineage fallback |

**`codon_tables.json` schema**

```json
{
  "<taxid>": {
    "source": "kazusa",
    "scientific_name": "optional string",
    "frequencies": {
      "TTT": 0.18,
      "TTC": 0.42
    }
  }
}
```

- Keys: DNA codons (T not U), sense codons only used for CAI.
- Frequencies: relative synonymous usage (0–1) per codon as in Kazusa CUTG.

**`taxonomy_parents.json` schema**

```json
{
  "562": "561",
  "561": "543"
}
```

### 1.5 Build-time inputs (`build_codon_reference`)

| Input | Description |
|-------|-------------|
| `taxids` | List of NCBI taxonomy IDs to fetch from Kazusa |
| `kraken_output` | Optional; union of all classified taxids from a Kraken file |
| `include_taxonomy` | If true, download NCBI taxdump and write `taxonomy_parents.json` |

---

## 2. Configuration parameters

### 2.1 Engineered k-mer scan

| Parameter | CLI flag | Default | Description |
|-----------|----------|---------|-------------|
| `window_size` | `--window_size` | 200 | Nucleotide window for k-mer voting |
| `engineered_kmer_threshold` | `--threshold` | 25 | Min engineered (32630) k-mers in window to label **Synthetic** |
| `threads` | `--threads` | 4 | Kraken2 threads (when `run_kraken=True`) |

Effective k-mers per window: `window_size - 21 + 1` (Kraken k=21 minimizers in window logic).

### 2.2 Codon adaptation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `kmer_len` | 35 | K-mer length for per-nucleotide taxid projection |
| `include_read_ids` | all Natural | If set, only these read IDs are scored |
| `codon_usage_dir` | user data dir | Pre-built reference directory |

---

## 3. Outputs

### 3.1 Engineered screening report (TSV-like text)

**Path:** user-provided `output_report_path` (e.g. `report.txt`)

| Column 1 | Column 2 | Meaning |
|----------|----------|---------|
| `Natural` or `Synthetic` | Read ID | Per-read engineering label |

**Example**

```text
Natural	read_natural_1
Synthetic	read_synthetic_1
```

### 3.2 Codon usage report (TSV)

**Path:** `--codon-usage-output` or `<report>.codon_usage.tsv`

Only **Natural** reads are included when run via full `screen` workflow.

| Column | Type | Description |
|--------|------|-------------|
| Read_ID | string | Read identifier |
| Overall_TaxID | string | Kraken assignment (column 2) |
| CDS_Strand | `+` / `-` | Strand of longest ORF |
| CDS_Range | string | `start-end` (0-based half-open on + strand) |
| CDS_TaxID | string | Dominant taxid over CDS region |
| Host_TaxID | string | Dominant taxid over flanking regions |
| Reference_TaxID | string | Taxid whose codon table was used (after lineage resolve) |
| CDS_Len_bp | int | CDS length in bp |
| CAI_vs_Host | float or `NA` | Codon Adaptation Index vs host reference (0–1) |

**Example**

```text
Read_ID	Overall_TaxID	CDS_Strand	CDS_Range	CDS_TaxID	Host_TaxID	Reference_TaxID	CDS_Len_bp	CAI_vs_Host
read_natural_1	562	+	0-12	562	562	562	12	0.7234
```

### 3.3 Library return types

#### `ScreenResult`

| Field | Type | Description |
|-------|------|-------------|
| `engineered_scan` | `EngineeredScanResult` | All read labels and counts |
| `codon_adaptation` | `list[CodonAdaptationResult]` | Natural reads only (if enabled) |
| `engineered_report_path` | `Path \| None` | Written engineered report |
| `codon_usage_report_path` | `Path \| None` | Written codon TSV |

#### `CodonAdaptationResult`

| Field | Type | Description |
|-------|------|-------------|
| `read_id` | str | Read ID |
| `overall_taxid` | str | Kraken taxid |
| `cds_strand` | str | `+` or `-` |
| `cds_start`, `cds_end` | int | CDS coordinates |
| `cds_taxid` | str | Taxid over CDS |
| `host_taxid` | str | Taxid over flanks |
| `reference_taxid` | str \| None | Resolved reference taxid |
| `cds_len_bp` | int | CDS length |
| `cai_vs_host` | float \| None | CAI score |

#### `ReadEngineeringLabel`

| Field | Type | Values |
|-------|------|--------|
| `read_id` | str | |
| `label` | str | `"Natural"` \| `"Synthetic"` |

#### `BuildCodonReferenceResult`

| Field | Type | Description |
|-------|------|-------------|
| `data_dir` | Path | Output directory |
| `taxids_requested` | list[str] | Input taxid set |
| `taxids_added` | list[str] | Newly fetched |
| `taxids_skipped` | list[str] | Already present |
| `taxids_failed` | list[str] | Kazusa fetch/parse failed |

---

## 4. Behavioural rules

### 4.1 Engineered vs Natural

- Each classified read with valid k-mer block data receives **Synthetic** if sliding-window count of taxid `32630` ≥ threshold; otherwise **Natural**.
- Unclassified reads (`U`) are still listed in the engineered report if present in Kraken output; codon analysis skips them.

### 4.2 Codon analysis gating

- Full pipeline: codon CAI runs **only** for reads labeled **Natural**.
- If zero Natural reads, codon output is omitted (empty list; no TSV unless explicitly requested elsewhere).

### 4.3 Airgapped / reference errors

| Condition | Exception |
|-----------|-----------|
| `codon_tables.json` missing | `CodonReferenceNotFoundError` |
| Host taxid not in store (no lineage match) | `MissingCodonReferenceError` |

No network access occurs during `analyze_codon_adaptation` or `run_screen` codon steps.

### 4.4 CAI computation

- Reference weights: Sharp & Li adaptiveness (per-codon weight / max within synonymous family).
- CAI = geometric mean of codon weights across in-frame codons in longest ORF (6-frame search).
- Longest ORF: longest segment between stop codons per frame; stops excluded from CDS sequence.

---

## 5. CLI command summary

| Command | Purpose |
|---------|---------|
| `plasmidScreen.py screen` | Kraken (optional) + engineered scan + codon usage |
| `plasmidScreen.py build-codon-db build` | Offline reference build |

### `screen` arguments

```text
screen FASTA_FILE OUTPUT_REPORT_PATH KRAKEN_RAW_OUTPUT [KRAKEN_DB_PATH]
  --window_size INT
  --threshold INT
  --codon-usage-output PATH
  --codon-usage-dir PATH
  --threads INT
```

### `build-codon-db build` arguments

```text
build-codon-db build
  --output-dir PATH
  --taxids "9606,511145"
  --taxids-file PATH
  --kraken-output PATH
  --skip-taxonomy
  --taxdump-dir PATH
  --timeout INT
```

---

## 6. Library API summary

```python
from plasmidScreen import (
    build_codon_reference,       # -> BuildCodonReferenceResult
    analyze_codon_adaptation,    # -> list[CodonAdaptationResult]
    run_screen,                  # -> ScreenResult
    write_codon_adaptation_tsv,  # -> tuple[str, list[CodonAdaptationResult]]
    CodonUsageStore,
    taxids_from_kraken_output,
)
```

Recommended airgapped workflow:

1. Networked: `build_codon_reference(dir, taxids=...)` or from Kraken taxid list.
2. Copy `dir` to airgapped environment.
3. Run Kraken2 locally → `analyze_codon_adaptation(reads, kraken.out, codon_usage_dir=dir)` or `run_screen(..., run_kraken=False)`.
