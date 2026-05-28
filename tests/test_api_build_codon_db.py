from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.api import build_codon_database


@patch("plasmidScreen.api.default_reference_taxids")
@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_defaults_when_no_taxids(mock_build, mock_defaults, tmp_path: Path) -> None:
    mock_defaults.return_value = ["9606", "511145"]
    build_codon_database(output_dir=tmp_path, download_csdb=False)
    args, kwargs = mock_build.call_args
    assert list(args[1]) == ["9606", "511145"]
    assert kwargs["use_default_taxids"] is False


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
    # args[1] is the taxid_list passed to build_codon_reference
    assert set(args[1]) == {"10090", "9606", "511145", "562"}
