"""Offline codon reference builder (network access only during build)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

from plasmidScreen.src.codon_usage.codon_usage_db import (
    CODON_TABLES_FILE,
    CodonUsageStore,
)
from plasmidScreen.src.codon_usage.codon_usage_hive import (
    HIVE_MIN_CDS,
    HIVE_ORGANELLE,
    default_hive_species_path,
    download_hive_species_tsv,
    import_hive_species_tsv,
)
from plasmidScreen.src.codon_usage.codon_usage_sources import (
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
        gene_set: GeneSet = GeneSet.NUCLEAR,
) -> BuildCodonReferenceResult:
    """
    Build ``codon_tables.json`` (uses existing ``taxonomy_parents.json`` when present).

    Always imports CSDB first, then fills gaps from the CoCoPUTs/HIVE genomic
    species table (same data dir as CSDB: download on first build unless
    ``download_csdb=False``, which also skips HIVE download and requires the
    cached TSV to exist). HIVE rows use ``Organelle=genomic`` and ``# CDS >= 1``.
    """
    data_dir = Path(data_dir)
    requested = sorted({str(t) for t in taxids if str(t) not in ("0", "")}) if taxids else []
    # When the user did not pass taxids, CSDB imports everything and HIVE fills
    # any remaining gaps (taxids=None), rather than limiting HIVE to CSDB ids.
    selective = bool(requested)
    taxids_requested = list(requested)

    added: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []

    store = CodonUsageStore.writable(data_dir)
    parents = store.taxonomy_parents() if store.has_taxonomy() else {}

    archive_path = Path(csdb_archive) if csdb_archive else default_csdb_archive_path()
    if download_csdb and not archive_path.is_file():
        download_csdb_archive(archive_path)

    if not selective:
        csdb_taxid_list = all_csdb_taxids(archive_path, gene_set=gene_set)
        logging.info(
            "No taxids specified; importing all %d taxid(s) from CSDB archive",
            len(csdb_taxid_list),
        )
        taxids_requested = csdb_taxid_list
        csdb_added, csdb_skipped = import_all_csdb_from_archive(
            store, archive_path, gene_set=gene_set
        )
        added.extend(csdb_added)
        skipped.extend(csdb_skipped)
    else:
        csdb_added, csdb_skipped, csdb_failed = import_csdb_taxids(
            store,
            archive_path,
            requested,
            gene_set=gene_set,
            parents=parents or None,
        )
        added.extend(csdb_added)
        skipped.extend(csdb_skipped)
        failed.extend(csdb_failed)

    hive_path = default_hive_species_path()
    if download_csdb and not hive_path.is_file():
        download_hive_species_tsv(hive_path)
    if not hive_path.is_file():
        raise FileNotFoundError(
            f"HIVE/CoCoPUTs species TSV not found at {hive_path}. "
            "Run with network download enabled, or place o537-genbank_species.tsv "
            "in the PlasmidScreen data directory."
        )

    logging.info(
        "Importing HIVE/CoCoPUTs genomic species tables from %s",
        hive_path,
    )
    hive_added, hive_skipped, hive_failed = import_hive_species_tsv(
        store,
        hive_path,
        taxids=requested if selective else None,
        parents=parents or None,
        organelle=HIVE_ORGANELLE,
        min_cds=HIVE_MIN_CDS,
        only_missing=True,
        save_every=10000,
    )
    added.extend(hive_added)
    skipped.extend(hive_skipped)
    recovered = set(hive_added)
    failed = [t for t in failed if t not in recovered]
    failed.extend(t for t in hive_failed if t not in recovered and t not in set(added))

    store.save()

    if not (data_dir / CODON_TABLES_FILE).exists() and not added and not skipped:
        raise RuntimeError(
            f"No codon tables written to {data_dir}. "
            "Check CSDB archive / HIVE species TSV and requested taxids."
        )

    return BuildCodonReferenceResult(
        data_dir=data_dir,
        taxids_requested=taxids_requested,
        taxids_added=sorted(set(added)),
        taxids_skipped=sorted(set(skipped)),
        taxids_failed=sorted(set(failed)),
    )
