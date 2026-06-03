"""Tests for Codon Statistics Database (CSDB) import."""
from __future__ import annotations

import io
import json
import tarfile
from pathlib import Path

import pytest

from plasmidScreen.lib.codon_usage_db import CodonUsageStore
from plasmidScreen.lib.codon_usage_sources import (
    build_csdb_ancestor_map,
    import_csdb_taxids,
    index_csdb_archive,
    parse_csdb_codon_statistics,
    read_csdb_table_from_archive,
    resolve_csdb_source_taxid,
)

SAMPLE_TSV = """\
Organism = Homo sapiens
Taxonomy ID = 9606

CODON\tAmino acid\tPreferred Codon\tTotal num\tFrequency in 1000\tFraction\tRSCU
GCT\tAla\tTRUE\t100\t10\t0.25\t1.0
GCC\tAla\tFALSE\t100\t10\t0.25\t1.0
GCA\tAla\tTRUE\t100\t10\t0.25\t1.0
GCG\tAla\tFALSE\t100\t10\t0.25\t1.0
TAA\tTer\tFALSE\t1\t1\t0.33\t1.0
TAG\tTer\tFALSE\t1\t1\t0.33\t1.0
TGA\tTer\tFALSE\t1\t1\t0.33\t1.0
"""


@pytest.fixture
def csdb_archive(tmp_path: Path) -> Path:
    """Minimal CSDB tar with human (9606) and E. coli K-12 (511145)."""
    archive = tmp_path / "csdb_test.tar.gz"
    entries = {
        "9606": SAMPLE_TSV,
        "511145": SAMPLE_TSV.replace("9606", "511145").replace("Homo sapiens", "E. coli"),
    }
    with tarfile.open(archive, "w:gz") as tar:
        for taxid, body in entries.items():
            data = body.encode()
            info = tarfile.TarInfo(name=f"data/{taxid}/nuclear_codon_statistics.tsv")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return archive


def test_parse_csdb_codon_statistics() -> None:
    freqs = parse_csdb_codon_statistics(SAMPLE_TSV)
    assert freqs["GCT"] == 0.25
    assert "TAA" not in freqs
    assert len(freqs) == 4


def test_index_csdb_archive(csdb_archive: Path) -> None:
    taxids = index_csdb_archive(csdb_archive, rebuild=True)
    assert taxids == {"9606", "511145"}
    cache = csdb_archive.parent / "csdb_taxid_index.json"
    assert cache.exists()


def test_resolve_csdb_source_taxid_with_ancestor_map() -> None:
    csdb_taxids = {"511145"}
    parents = {"511145": "562", "562": "561", "561": "561"}
    ancestor_map = build_csdb_ancestor_map(csdb_taxids, parents)
    assert resolve_csdb_source_taxid("511145", csdb_taxids, ancestor_map) == "511145"
    assert resolve_csdb_source_taxid("562", csdb_taxids, ancestor_map) == "511145"


def test_read_csdb_table_from_archive(csdb_archive: Path) -> None:
    freqs = read_csdb_table_from_archive(csdb_archive, "9606")
    assert freqs["GCT"] == 0.25


def test_import_all_csdb_from_archive(tmp_path: Path, csdb_archive: Path) -> None:
    from plasmidScreen.lib.codon_usage_sources import import_all_csdb_from_archive

    store = CodonUsageStore.writable(tmp_path / "codon_usage")
    added, skipped = import_all_csdb_from_archive(store, csdb_archive)
    assert set(added) == {"9606", "511145"}
    assert skipped == []
    assert store.has_codon_table("9606")


def test_import_all_csdb_from_archive_skips_existing(
    tmp_path: Path, csdb_archive: Path
) -> None:
    from plasmidScreen.lib.codon_usage_sources import import_all_csdb_from_archive

    store = CodonUsageStore.writable(tmp_path / "codon_usage")
    store.set_codon_table("9606", {"GCT": 0.25}, source="test")
    store.save()
    added, skipped = import_all_csdb_from_archive(store, csdb_archive)
    assert added == ["511145"]
    assert skipped == ["9606"]


def test_all_csdb_taxids(csdb_archive: Path) -> None:
    from plasmidScreen.lib.codon_usage_sources import all_csdb_taxids

    assert all_csdb_taxids(csdb_archive, rebuild_index=True) == ["511145", "9606"]


def test_import_csdb_taxids(tmp_path: Path, csdb_archive: Path) -> None:
    store = CodonUsageStore.writable(tmp_path / "codon_usage")
    parents = {"511145": "562", "562": "561", "561": "561"}
    store._parents = parents
    store._dirty = True

    added, skipped, failed = import_csdb_taxids(
        store,
        csdb_archive,
        ["9606", "562", "999999"],
        parents=parents,
    )
    assert "9606" in added
    assert "562" in added
    assert "999999" in failed
    assert store.has_codon_table("562")
    assert store.get_frequencies("562")["GCT"] == 0.25
