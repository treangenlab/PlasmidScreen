"""Tests for codon usage reference store (JSON, airgapped)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from plasmidScreen.lib.codon_usage_db import (
    CodonUsageStore,
    cai_weights_from_frequencies,
    parse_taxonomy_nodes,
)
from plasmidScreen.lib.exceptions import CodonReferenceNotFoundError, MissingCodonReferenceError


def test_cai_weights_max_within_amino_acid(uniform_frequencies: dict[str, float]) -> None:
    weights = cai_weights_from_frequencies(uniform_frequencies)
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    by_aa: dict[str, list[str]] = {}
    for codon, w in weights.items():
        aa = forward.get(codon)
        if aa and aa != "*":
            by_aa.setdefault(aa, []).append(w)
    for aa, ws in by_aa.items():
        assert max(ws) == pytest.approx(1.0), f"AA {aa} should have a max weight of 1.0"


def test_store_load_and_resolve_lineage(codon_store: CodonUsageStore) -> None:
    assert codon_store.has_codon_table("562")
    ref = codon_store.resolve_reference_taxid("562")
    assert ref == "562"
    weights = codon_store.get_cai_weights_for_host("562")
    assert weights is not None
    assert "ATG" in weights


def test_missing_host_raises(codon_store: CodonUsageStore) -> None:
    with pytest.raises(MissingCodonReferenceError) as exc:
        codon_store.require_host_taxids(["99999"])
    assert "99999" in exc.value.missing_host_taxids


def test_missing_store_raises(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(CodonReferenceNotFoundError):
        CodonUsageStore.load(empty)


def test_writable_store_roundtrip(tmp_path: Path, uniform_frequencies: dict[str, float]) -> None:
    data_dir = tmp_path / "writable"
    store = CodonUsageStore.writable(data_dir)
    store.set_codon_table("123", uniform_frequencies, source="test")
    store.save()

    loaded = CodonUsageStore.load(data_dir)
    assert loaded.get_frequencies("123") == uniform_frequencies
