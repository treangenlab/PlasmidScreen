"""Codon Statistics Database (CSDB) import for offline reference builds."""
from __future__ import annotations

import csv
import io
import json
import logging
import tarfile
import urllib.request
from pathlib import Path
from typing import Iterable

from Bio.Data import CodonTable
from plasmidScreen.lib.codon_usage_db import CodonUsageStore
from plasmidScreen.lib.funcs import get_default_db_path

CSDB_ARCHIVE_URL = "http://codonstatsdb.unr.edu/codonstatsdb_March2022.tar.gz"
CSDB_ARCHIVE_FILENAME = "codonstatsdb_March2022.tar.gz"
CSDB_INDEX_FILENAME = "csdb_taxid_index.json"

GENE_SET_FILES = {
    "nuclear": "nuclear_codon_statistics.tsv",
    "ribosomal": "ribosomal_codon_statistics.tsv",
    "mitochondrial": "mitochondrial_codon_statistics.tsv",
    "plastid": "plastid_codon_statistics.tsv",
}

_CODON_TO_AA = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
for _stop in ("TAA", "TAG", "TGA"):
    _CODON_TO_AA[_stop] = "*"


def default_csdb_archive_path(app_name: str = "PlasmidScreen") -> Path:
    """Default cache location for the CSDB bulk tar archive."""
    return Path(get_default_db_path(app_name)) / CSDB_ARCHIVE_FILENAME


def default_csdb_index_path(archive_path: Path | None = None) -> Path:
    archive = archive_path or default_csdb_archive_path()
    return archive.parent / CSDB_INDEX_FILENAME


def download_csdb_archive(
    dest_path: Path | None = None,
    *,
    url: str = CSDB_ARCHIVE_URL,
) -> Path:
    """Download the CSDB flat-file release (~5.2 GB). Build-time only."""
    dest = Path(dest_path or default_csdb_archive_path())
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        logging.info("CSDB archive already present at %s", dest)
        return dest
    logging.info("Downloading CSDB archive from %s (this is ~5.2 GB)", url)
    urllib.request.urlretrieve(url, dest)
    logging.info("CSDB archive saved to %s", dest)
    return dest


def parse_csdb_codon_statistics(tsv_text: str) -> dict[str, float]:
    """
    Parse CSDB nuclear_codon_statistics.tsv into DNA codon -> Fraction.

    Uses the synonymous Fraction column (0–1 per amino acid), suitable for CAI weights.
    """
    lines = tsv_text.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("CODON\t") or line.startswith("Codon\t"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("CSDB codon statistics TSV missing CODON header row")

    reader = csv.DictReader(io.StringIO("\n".join(lines[header_idx:])), delimiter="\t")
    frequencies: dict[str, float] = {}
    for row in reader:
        codon = row.get("CODON", row.get("Codon", "")).strip().upper()
        if not codon or len(codon) != 3:
            continue
        aa = _CODON_TO_AA.get(codon)
        if aa is None or aa == "*":
            continue
        fraction = row.get("Fraction", row.get("fraction", "")).strip()
        if not fraction:
            continue
        frequencies[codon] = float(fraction)

    if not frequencies:
        raise ValueError("No sense codon frequencies parsed from CSDB TSV")
    return frequencies


def _member_path(taxid: str, gene_set: str) -> str:
    filename = GENE_SET_FILES[gene_set]
    return f"data/{taxid}/{filename}"


def index_csdb_archive(
    archive_path: Path,
    *,
    gene_set: str = "nuclear",
    index_path: Path | None = None,
    rebuild: bool = False,
) -> set[str]:
    """List taxids present in the CSDB archive (cached as JSON beside the tar)."""
    archive_path = Path(archive_path)
    cache = Path(index_path or default_csdb_index_path(archive_path))
    if cache.exists() and not rebuild:
        data = json.loads(cache.read_text(encoding="utf-8"))
        return set(data.get("taxids", []))

    suffix = GENE_SET_FILES[gene_set]
    taxids: set[str] = set()
    logging.info("Indexing CSDB archive at %s (one-time scan)", archive_path)
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar.getmembers():
            name = member.name
            if not name.endswith(suffix):
                continue
            parts = name.split("/")
            if len(parts) >= 3 and parts[0] == "data":
                taxids.add(parts[1])

    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(
        json.dumps(
            {
                "archive": str(archive_path.resolve()),
                "gene_set": gene_set,
                "taxids": sorted(taxids),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    logging.info("CSDB index: %d taxid(s) -> %s", len(taxids), cache)
    return taxids


def build_csdb_ancestor_map(
    csdb_taxids: Iterable[str | int],
    parents: dict[str, str],
) -> dict[str, str]:
    """
    Map any NCBI taxid to a CSDB table taxid by walking up the lineage.

    For each CSDB entry, every ancestor taxid maps to that entry (deepest CSDB child wins
    when multiple exist under the same ancestor).
    """
    ancestor_map: dict[str, str] = {}
    for csdb_id in csdb_taxids:
        current = str(csdb_id)
        seen: set[str] = set()
        while current and current not in seen:
            if current not in ancestor_map:
                ancestor_map[current] = str(csdb_id)
            seen.add(current)
            parent = parents.get(current)
            if parent is None or parent == current:
                break
            current = parent
    return ancestor_map


def resolve_csdb_source_taxid(
    taxid: str | int,
    csdb_taxids: set[str],
    ancestor_map: dict[str, str],
) -> str | None:
    """Resolve a requested taxid to a CSDB folder taxid."""
    key = str(taxid)
    if key in csdb_taxids:
        return key
    return ancestor_map.get(key)


def read_csdb_table_from_archive(
    archive_path: Path,
    source_taxid: str,
    *,
    gene_set: str = "nuclear",
) -> dict[str, float]:
    """Extract and parse one species codon table from the CSDB tar."""
    if gene_set not in GENE_SET_FILES:
        raise ValueError(f"Unknown gene_set={gene_set!r}; choose from {list(GENE_SET_FILES)}")

    member_name = _member_path(source_taxid, gene_set)
    with tarfile.open(archive_path, "r:gz") as tar:
        try:
            member = tar.getmember(member_name)
        except KeyError as exc:
            raise FileNotFoundError(
                f"No {member_name} in CSDB archive"
            ) from exc
        extracted = tar.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"Could not read {member_name} from CSDB archive")
        tsv_text = extracted.read().decode(errors="replace")

    return parse_csdb_codon_statistics(tsv_text)


def import_csdb_taxids(
    store: CodonUsageStore,
    archive_path: Path,
    taxids: Iterable[str | int],
    *,
    gene_set: str = "nuclear",
    parents: dict[str, str] | None = None,
    save_every: int = 50,
) -> tuple[list[str], list[str], list[str]]:
    """
    Import codon usage tables from a local CSDB tar archive.

    Returns (added, skipped, failed) for requested taxids.
    """
    archive_path = Path(archive_path)
    if not archive_path.is_file():
        raise FileNotFoundError(
            f"CSDB archive not found: {archive_path}. "
            f"Download from {CSDB_ARCHIVE_URL} or pass csdb_archive= to build_codon_reference()."
        )

    taxid_list = sorted({str(t) for t in taxids if str(t) not in ("0", "")})
    skipped = [t for t in taxid_list if store.has_codon_table(t)]
    to_import = [t for t in taxid_list if not store.has_codon_table(t)]

    csdb_taxids = index_csdb_archive(archive_path, gene_set=gene_set)
    ancestor_map: dict[str, str] = {}
    if parents:
        ancestor_map = build_csdb_ancestor_map(csdb_taxids, parents)

    added: list[str] = []
    failed: list[str] = []
    for i, taxid in enumerate(to_import, start=1):
        source_id = resolve_csdb_source_taxid(taxid, csdb_taxids, ancestor_map)
        if source_id is None:
            logging.warning(
                "No CSDB table for taxid %s (not in archive and no ancestor match)",
                taxid,
            )
            failed.append(taxid)
            continue
        try:
            frequencies = read_csdb_table_from_archive(
                archive_path, source_id, gene_set=gene_set
            )
            store.set_codon_table(taxid, frequencies, source="csdb")
            added.append(taxid)
            if source_id != taxid:
                logging.info(
                    "Imported CSDB taxid %s for requested taxid %s", source_id, taxid
                )
            else:
                logging.info("Imported CSDB codon table for taxid %s", taxid)
        except (FileNotFoundError, ValueError, OSError) as exc:
            logging.warning("CSDB import failed for taxid %s: %s", taxid, exc)
            failed.append(taxid)

        if save_every and i % save_every == 0:
            store.save()
            logging.info("Checkpoint: %d/%d CSDB imports processed", i, len(to_import))

    store.save()
    return added, skipped, failed
