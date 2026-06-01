"""Offline codon reference builder (network access only during build)."""
from __future__ import annotations

import logging
import tarfile
import urllib.request
from importlib import resources
from pathlib import Path
from typing import Iterable

from plasmidScreen.lib.codon_usage_db import (
    CODON_TABLES_FILE,
    CodonUsageStore,
    parse_taxonomy_nodes,
)
from plasmidScreen.lib.codon_usage_sources import (
    all_csdb_taxids,
    default_csdb_archive_path,
    download_csdb_archive,
    import_all_csdb_from_archive,
    import_csdb_taxids,
)
from plasmidScreen.lib.models import BuildCodonReferenceResult
from plasmidScreen.lib.types import GeneSet

NCBI_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def _load_taxids_from_package_file(filename: str) -> list[str]:
    try:
        text = resources.files("plasmidScreen.data").joinpath(filename).read_text(encoding="utf-8")
    except (FileNotFoundError, ModuleNotFoundError, TypeError, OSError):
        return []
    taxids: list[str] = []
    for line in text.splitlines():
        line = line.strip().split("#")[0].strip()
        if line:
            taxids.append(line)
    return taxids


def default_reference_taxids() -> list[str]:
    """
    Optional curated subset (~150 taxids) from package data files.

    Not used for default builds; see ``all_csdb_taxids()`` / full CSDB import instead.
    """
    combined = set(_load_taxids_from_package_file("common_codon_taxids.txt"))
    combined.update(_load_taxids_from_package_file("default_codon_taxids.txt"))
    if not combined:
        combined = {
            "9606", "10090", "511145", "4932", "7227", "6239",
            "287", "1282", "1313", "1288",
        }
    return sorted(combined)


def download_ncbi_taxdump(dest_dir: Path) -> Path:
    """Download and extract nodes.dmp from NCBI taxdump."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    archive = dest_dir / "taxdump.tar.gz"
    if not archive.exists():
        logging.info("Downloading NCBI taxdump from %s", NCBI_TAXDUMP_URL)
        urllib.request.urlretrieve(NCBI_TAXDUMP_URL, archive)

    nodes_path = dest_dir / "nodes.dmp"
    if not nodes_path.exists():
        logging.info("Extracting nodes.dmp from taxdump archive")
        with tarfile.open(archive, "r:gz") as tar:
            try:
                tar.extract("nodes.dmp", path=dest_dir, filter="data")
            except TypeError:
                tar.extract("nodes.dmp", path=dest_dir)

    return nodes_path


def build_codon_reference(
    data_dir: str | Path,
    taxids: Iterable[str | int] | None = None,
    *,
    include_taxonomy: bool = True,
    taxdump_dir: str | Path | None = None,
    csdb_archive: str | Path | None = None,
    download_csdb: bool = True,
    gene_set: GeneSet = "nuclear",
) -> BuildCodonReferenceResult:
    """
    Build codon_tables.json (and optional taxonomy_parents.json) for airgapped use.

    When ``taxids`` is omitted or empty, imports **every taxid present in the CSDB
    archive** for the selected ``gene_set``. Pass explicit taxids to import a subset.
    """
    data_dir = Path(data_dir)
    taxid_list = sorted({str(t) for t in taxids if str(t) not in ("0", "")}) if taxids else []

    archive_path = Path(csdb_archive) if csdb_archive else default_csdb_archive_path()
    if download_csdb and not archive_path.is_file():
        download_csdb_archive(archive_path)

    import_all = not taxid_list
    if import_all:
        taxid_list = all_csdb_taxids(archive_path, gene_set=gene_set)
        logging.info(
            "No taxids specified; importing all %d taxid(s) from CSDB archive",
            len(taxid_list),
        )

    added: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []

    store = CodonUsageStore.writable(data_dir)

    if include_taxonomy and not store.has_taxonomy():
        tdir = Path(taxdump_dir) if taxdump_dir else data_dir.parent / "taxdump"
        nodes = download_ncbi_taxdump(tdir)
        count = store.load_taxonomy_from_nodes(nodes)
        store.save()
        logging.info("Loaded %d taxonomy parent links", count)

    if import_all:
        added, skipped = import_all_csdb_from_archive(
            store, archive_path, gene_set=gene_set
        )
    else:
        parents = store.taxonomy_parents() if store.has_taxonomy() else {}
        csdb_added, csdb_skipped, csdb_failed = import_csdb_taxids(
            store,
            archive_path,
            taxid_list,
            gene_set=gene_set,
            parents=parents or None,
        )
        added.extend(csdb_added)
        skipped.extend(csdb_skipped)
        failed.extend(csdb_failed)

    store.save()

    if not (data_dir / CODON_TABLES_FILE).exists() and not added and not skipped:
        raise RuntimeError(
            f"No codon tables written to {data_dir}. "
            "Check CSDB archive path and requested taxids."
        )

    return BuildCodonReferenceResult(
        data_dir=data_dir,
        taxids_requested=taxid_list,
        taxids_added=added,
        taxids_skipped=skipped,
        taxids_failed=failed,
    )
