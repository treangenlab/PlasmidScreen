"""Offline codon reference builder (network access only during build)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

from plasmidScreen.lib.codon_usage_db import (
    CODON_TABLES_FILE,
    CodonUsageStore,
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


def build_codon_reference(
        data_dir: str | Path,
        taxids: Iterable[str | int] | None = None,
        *,
        csdb_archive: str | Path | None = None,
        download_csdb: bool = True,
        gene_set: GeneSet = "nuclear",
) -> BuildCodonReferenceResult:
    """
    Build codon_tables.json (and optional taxonomy_parents.json).



    :param csdb_archive: path of where the archieve should live, by default it's in the APPs directory,
    :param download_csdb: download the CSDB?
    :param gene_set: As there are several categories, by default plasmidscreen uses everything denoted as nuclear.
                     GeneSet specifies the rest.


    Build the codon usage reference for CAI scoring. It attempts to use a
    taxids file if provided to grab the codon references. If none is provided, it imports
    **every** taxid in the CSDB archive for ``gene_set``. Downloads CSDB archive if download_csdb is provided as true.
    Writes ``codon_tables.json`` and optionally ``taxonomy_parents.json``
    under ``output_dir`` (default: PlasmidScreen user data ``codon_usage/``).
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
