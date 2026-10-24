"""Tests for codon reference build helpers."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference


@patch("plasmidScreen.src.codon_usage.codon_usage_build.import_all_csdb_from_archive")
@patch("plasmidScreen.src.codon_usage.codon_usage_build.all_csdb_taxids")
def test_build_codon_reference_defaults_to_all_csdb(
    mock_all_taxids,
    mock_import_all,
    tmp_path: Path,
) -> None:
    mock_all_taxids.return_value = ["9606", "511145"]
    mock_import_all.return_value = (["9606", "511145"], [])
    archive = tmp_path / "csdb.tar.gz"
    archive.write_bytes(b"placeholder")

    result = build_codon_reference(
        tmp_path / "codon_usage",
        csdb_archive=archive,
        download_csdb=False,
    )
    mock_all_taxids.assert_called_once()
    mock_import_all.assert_called_once()
    assert result.taxids_requested == ["9606", "511145"]
    assert result.taxids_added == ["9606", "511145"]


@patch("plasmidScreen.src.codon_usage.codon_usage_build.import_csdb_taxids")
@patch("plasmidScreen.src.codon_usage.codon_usage_build.import_all_csdb_from_archive")
def test_build_codon_reference_with_taxids_calls_import_csdb(
    mock_import_all,
    mock_import_taxids,
    tmp_path: Path,
) -> None:
    mock_import_taxids.return_value = (["562"], [], [])
    archive = tmp_path / "csdb.tar.gz"
    archive.write_bytes(b"placeholder")

    result = build_codon_reference(
        tmp_path / "codon_usage",
        taxids=["562"],
        csdb_archive=archive,
        download_csdb=False,
    )
    mock_import_all.assert_not_called()
    mock_import_taxids.assert_called_once()
    assert result.taxids_requested == ["562"]
    assert result.taxids_added == ["562"]
