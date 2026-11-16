"""Tests for NCBI accession resolution and GenBank caching helpers."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.lib.genbank_fetch import (
    fetch_genbank_texts,
    resolve_ncbi_accession,
)


def test_resolve_plain_ncbi_accession() -> None:
    r = resolve_ncbi_accession("MT024802.1")
    assert r is not None
    assert r.ncbi_id == "MT024802.1"
    assert r.slice_start is None


def test_resolve_univec_slice() -> None:
    hid = "gnl|uv|X66730.1:1-2687-49|kraken:taxid|32630"
    r = resolve_ncbi_accession(hid)
    assert r is not None
    assert r.hit_id == hid
    assert r.ncbi_id == "X66730.1"
    assert r.slice_start == 1
    assert r.slice_end == 2687


def test_resolve_imgpr_skipped() -> None:
    assert (
        resolve_ncbi_accession(
            "IMGPR_plasmid_3300053491_006474|3300053491|Ga0557493_093524|kraken:taxid|36549"
        )
        is None
    )


def test_collect_ncbi_ids_from_fasta(tmp_path: Path) -> None:
    from plasmidScreen.lib.genbank_fetch import collect_ncbi_ids_from_fasta

    fa = tmp_path / "refs.fa"
    fa.write_text(
        ">MT024802.1 some plasmid\nACGT\n"
        ">gnl|uv|X66730.1:1-100-49|kraken:taxid|32630\nACGT\n"
        ">IMGPR_plasmid_1|kraken:taxid|36549\nACGT\n"
        ">MT024802.1 duplicate\nACGT\n",
        encoding="utf-8",
    )
    ids = collect_ncbi_ids_from_fasta(fa)
    assert ids == ["MT024802.1", "X66730.1"]


def test_fetch_genbank_texts_uses_cache(tmp_path: Path) -> None:
    gb = tmp_path / "X66730.1.gb"
    gb.write_text(
        "LOCUS       X66730\n"
        "ACCESSION   X66730\n"
        "VERSION     X66730.1\n"
        "FEATURES             Location/Qualifiers\n"
        "ORIGIN\n"
        "        1 a\n"
        "//\n",
        encoding="utf-8",
    )
    with patch("plasmidScreen.lib.genbank_fetch._efetch_batch") as mock_fetch:
        out = fetch_genbank_texts(["X66730.1"], cache_dir=tmp_path)
        mock_fetch.assert_not_called()
    assert "X66730.1" in out
    assert "ORIGIN" in out["X66730.1"]
