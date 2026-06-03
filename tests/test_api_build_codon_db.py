from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.api import build_codon_database


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_imports_all_csdb_when_no_taxids(
    mock_build, tmp_path: Path
) -> None:
    build_codon_database(output_dir=tmp_path, download_csdb=False)
    args, kwargs = mock_build.call_args
    assert args[1] is None
    assert kwargs["download_csdb"] is False


@patch("plasmidScreen.api.taxids_from_kraken_output")
@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_unions_inputs(mock_build, mock_from_kraken, tmp_path: Path) -> None:
    mock_from_kraken.return_value = {"562"}
    taxids_file = tmp_path / "taxids.txt"
    taxids_file.write_text("9606\n# comment\n511145\n", encoding="utf-8")

    build_codon_database(
        output_dir=tmp_path,
        taxids=["10090"],
        taxids_file=taxids_file,
        kraken_output=tmp_path / "kraken.out",
        download_csdb=False,
    )
    args, _kwargs = mock_build.call_args
    assert set(args[1]) == {"10090", "9606", "511145", "562"}
