"""Import HIVE / CoCoPUTs species codon-count TSVs into CodonUsageStore.

Default cache: ``~/.local/share/PlasmidScreen/o537-genbank_species.tsv``
(DNAHIVE object 537). Build always imports ``Organelle=genomic`` rows with
``# CDS >= 1``, converting raw codon counts to synonymous fractions for
``codon_tables.json``.
"""
from __future__ import annotations

import csv
import logging
import urllib.request
from pathlib import Path
from typing import Iterable, Iterator, Optional

from Bio.Data import CodonTable

from plasmidScreen.lib.funcs import get_default_db_path
from plasmidScreen.src.codon_usage.codon_usage_db import CodonUsageStore

# FDA DNAHIVE CoCoPUTs bulk species table (object 537).
HIVE_SPECIES_FILENAME = "o537-genbank_species.tsv"
HIVE_SPECIES_URL = (
    "https://dnahive.fda.gov/dna.cgi?cmd=objFile&filename=genbank_species.tsv&ids=537"
)
HIVE_ORGANELLE = "genomic"
HIVE_MIN_CDS = 1

_CODON_TO_AA = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
for _stop in ("TAA", "TAG", "TGA"):
    _CODON_TO_AA[_stop] = "*"

SENSE_CODONS: tuple[str, ...] = tuple(
    sorted(c for c, aa in _CODON_TO_AA.items() if aa != "*")
)

REQUIRED_HIVE_COLUMNS: tuple[str, ...] = (
    "Taxid",
    "Species",
    "Organelle",
    "# CDS",
    "# Codons",
) + SENSE_CODONS


def default_hive_species_path(app_name: str = "PlasmidScreen") -> Path:
    """Default cache location for the CoCoPUTs/HIVE species codon TSV."""
    return Path(get_default_db_path(app_name)) / HIVE_SPECIES_FILENAME


def download_hive_species_tsv(
    dest_path: Path | None = None,
    *,
    url: str = HIVE_SPECIES_URL,
) -> Path:
    """Download the HIVE/CoCoPUTs genbank species TSV. Build-time only."""
    dest = Path(dest_path or default_hive_species_path())
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        logging.info("HIVE species TSV already present at %s", dest)
        return dest
    logging.info("Downloading HIVE/CoCoPUTs species TSV from %s", url)
    urllib.request.urlretrieve(url, dest)
    logging.info("HIVE species TSV saved to %s", dest)
    return dest


def validate_hive_species_header(fieldnames: Iterable[str] | None) -> list[str]:
    """Return header list or raise ValueError if required columns are missing."""
    if not fieldnames:
        raise ValueError("HIVE species TSV has no header row")
    header = list(fieldnames)
    missing = [c for c in REQUIRED_HIVE_COLUMNS if c not in header]
    if missing:
        preview = ", ".join(missing[:12])
        more = len(missing) - min(len(missing), 12)
        suffix = f" (+{more} more)" if more > 0 else ""
        raise ValueError(
            f"HIVE species TSV missing required column(s): {preview}{suffix}"
        )
    return header


def counts_to_synonymous_fractions(counts: dict[str, float]) -> dict[str, float]:
    """Convert raw codon counts to synonymous Fraction (0–1 per amino acid)."""
    by_aa: dict[str, dict[str, float]] = {}
    for codon, count in counts.items():
        aa = _CODON_TO_AA.get(codon)
        if aa is None or aa == "*":
            continue
        by_aa.setdefault(aa, {})[codon] = max(float(count), 0.0)

    frequencies: dict[str, float] = {}
    for codon_counts in by_aa.values():
        total = sum(codon_counts.values())
        if total <= 0:
            continue
        for codon, count in codon_counts.items():
            frequencies[codon] = count / total
    if not frequencies:
        raise ValueError("No sense-codon counts to convert to frequencies")
    return frequencies


def _parse_count(raw: str) -> float:
    text = (raw or "").strip()
    if not text:
        return 0.0
    return float(text)


def parse_hive_species_row(row: dict[str, str]) -> tuple[str, str, dict[str, float], int]:
    """
    Parse one HIVE species row into (taxid, scientific_name, frequencies, n_cds).

    Raises ValueError if the row cannot yield codon frequencies.
    """
    taxid = (row.get("Taxid") or "").strip()
    if not taxid or taxid in ("0",):
        raise ValueError("missing Taxid")
    name = (row.get("Species") or "").strip() or None
    try:
        n_cds = int(float((row.get("# CDS") or "0").strip() or "0"))
    except ValueError as exc:
        raise ValueError(f"bad # CDS for taxid {taxid}") from exc

    counts: dict[str, float] = {}
    for codon in SENSE_CODONS:
        counts[codon] = _parse_count(row.get(codon, ""))
    frequencies = counts_to_synonymous_fractions(counts)
    return taxid, (name or ""), frequencies, n_cds


def iter_hive_species_rows(
    path: Path | str,
    *,
    organelle: str = "genomic",
    min_cds: int = 1,
) -> Iterator[tuple[str, str, dict[str, float], int]]:
    """Yield (taxid, scientific_name, frequencies, n_cds) for matching rows."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"HIVE species TSV not found: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_hive_species_header(reader.fieldnames)
        for row in reader:
            if organelle and (row.get("Organelle") or "").strip() != organelle:
                continue
            try:
                taxid, name, frequencies, n_cds = parse_hive_species_row(row)
            except ValueError:
                continue
            if n_cds < min_cds:
                continue
            yield taxid, name, frequencies, n_cds


def index_hive_species_taxids(
    path: Path | str,
    *,
    organelle: str = "genomic",
    min_cds: int = 1,
) -> set[str]:
    """Taxids with a usable table in the HIVE species TSV."""
    return {taxid for taxid, _name, _freq, _n in iter_hive_species_rows(
        path, organelle=organelle, min_cds=min_cds
    )}


def resolve_hive_source_taxid(
    taxid: str,
    hive_taxids: set[str],
    parents: dict[str, str] | None,
) -> Optional[str]:
    """Nearest ancestor (including self) present in ``hive_taxids``."""
    if taxid in ("0", ""):
        return None
    if taxid in hive_taxids:
        return taxid
    if not parents:
        return None
    current = taxid
    seen: set[str] = set()
    while current and current not in seen:
        seen.add(current)
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        if parent in hive_taxids:
            return parent
        current = parent
    return None


def import_hive_species_tsv(
    store: CodonUsageStore,
    path: Path | str,
    *,
    taxids: Iterable[str | int] | None = None,
    parents: dict[str, str] | None = None,
    organelle: str = "genomic",
    min_cds: int = 1,
    only_missing: bool = True,
    save_every: int = 500,
) -> tuple[list[str], list[str], list[str]]:
    """
    Import HIVE/CoCoPUTs species codon tables into ``store``.

    Parameters
    ----------
    taxids
        If set, only these hosts are considered (lineage fallback via ``parents``).
        If None, stream every matching row in the file.
    only_missing
        Skip taxids that already have a ``codon_tables.json`` entry (CSDB wins).
    parents
        Optional ``taxonomy_parents.json`` map for lineage resolve when ``taxids``
        is set.

    Returns
    -------
    added, skipped, failed
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"HIVE species TSV not found: {path}")

    # Validate header early.
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_hive_species_header(reader.fieldnames)

    added: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []

    if taxids is None:
        processed = 0
        for taxid, name, frequencies, _n_cds in iter_hive_species_rows(
            path, organelle=organelle, min_cds=min_cds
        ):
            if only_missing and store.has_codon_table(taxid):
                skipped.append(taxid)
                continue
            store.set_codon_table(
                taxid,
                frequencies,
                scientific_name=name or None,
                source="hive",
            )
            added.append(taxid)
            processed += 1
            if save_every and processed % save_every == 0:
                store.save()
                logging.info("HIVE checkpoint: %d table(s) imported", processed)

        store.save()
        logging.info(
            "HIVE import complete: added=%d skipped=%d from %s",
            len(added),
            len(skipped),
            path,
        )
        return sorted(set(added)), sorted(set(skipped)), failed

    requested = sorted({str(t) for t in taxids if str(t) not in ("0", "")})
    to_resolve = []
    for taxid in requested:
        if only_missing and store.has_codon_table(taxid):
            skipped.append(taxid)
        else:
            to_resolve.append(taxid)

    if not to_resolve:
        return [], sorted(set(skipped)), []

    parents_map = parents if parents is not None else (
        store.taxonomy_parents() if store.has_taxonomy() else {}
    )

    # Candidate ids = requested hosts plus their ancestor chains (one TSV pass).
    watch: set[str] = set(to_resolve)
    if parents_map:
        for taxid in to_resolve:
            current = taxid
            seen: set[str] = set()
            while current and current not in seen:
                seen.add(current)
                parent = parents_map.get(current)
                if parent is None or parent == current:
                    break
                watch.add(parent)
                current = parent

    source_tables: dict[str, tuple[str, dict[str, float]]] = {}
    for taxid, name, frequencies, _n_cds in iter_hive_species_rows(
        path, organelle=organelle, min_cds=min_cds
    ):
        if taxid in watch and taxid not in source_tables:
            source_tables[taxid] = (name, frequencies)

    hive_taxids = set(source_tables)
    processed = 0
    for taxid in to_resolve:
        source_id = resolve_hive_source_taxid(taxid, hive_taxids, parents_map or None)
        if source_id is None:
            failed.append(taxid)
            continue
        name, frequencies = source_tables[source_id]
        store.set_codon_table(
            taxid,
            frequencies,
            scientific_name=name or None,
            source="hive",
        )
        added.append(taxid)
        if source_id != taxid:
            logging.info(
                "Imported HIVE taxid %s for requested taxid %s", source_id, taxid
            )
        processed += 1
        if save_every and processed % save_every == 0:
            store.save()

    store.save()
    logging.info(
        "HIVE selective import: added=%d skipped=%d failed=%d from %s",
        len(added),
        len(skipped),
        len(failed),
        path,
    )
    return sorted(set(added)), sorted(set(skipped)), sorted(set(failed))
