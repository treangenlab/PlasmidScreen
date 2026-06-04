"""Tests for codon reference build helpers."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.src.codon_usage.codon_usage_build import build_codon_reference


@patch("plasmidScreen.lib.codon_usage_build.import_all_csdb_from_archive")
@patch("plasmidScreen.lib.codon_usage_build.all_csdb_taxids")
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
        include_taxonomy=False,
        csdb_archive=archive,
        download_csdb=False,
    )
    mock_all_taxids.assert_called_once()
    mock_import_all.assert_called_once()
    assert result.taxids_requested == ["9606", "511145"]
    assert result.taxids_added == ["9606", "511145"]
