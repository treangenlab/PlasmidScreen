"""Unit + real-CSDB integration tests for ``build_codon_database``."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from plasmidScreen.api import build_codon_database
from plasmidScreen.lib.models import BuildCodonReferenceResult
from plasmidScreen.src.codon_usage.codon_usage_db import (
    CODON_TABLES_FILE,
    CodonUsageStore,
)

# Must match ``REAL_CSDB_TAXIDS`` in conftest session fixture.
REAL_CSDB_TAXIDS = ("511145", "562", "9606")


def _fake_result(data_dir: Path) -> BuildCodonReferenceResult:
    return BuildCodonReferenceResult(
        data_dir=data_dir,
        taxids_requested=["562"],
        taxids_added=["562"],
        taxids_skipped=[],
        taxids_failed=[],
    )


# ---------------------------------------------------------------------------
# Fast unit tests (mocked build_codon_reference)
# ---------------------------------------------------------------------------


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_imports_all_csdb_when_no_taxids(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    """No taxids / taxids_file → taxid_list is None (import entire CSDB gene set)."""
    mock_build.return_value = _fake_result(tmp_path)

    result = build_codon_database(output_dir=tmp_path, download_csdb=False)

    mock_build.assert_called_once_with(
        tmp_path,
        None,
        csdb_archive=None,
        download_csdb=False,
        gene_set="nuclear",
    )
    assert result is mock_build.return_value


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_passes_sorted_taxids(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(
        output_dir=tmp_path,
        taxids=["562", 511145, "9606"],
        download_csdb=False,
    )

    assert mock_build.call_args.args[1] == ["511145", "562", "9606"]


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_filters_empty_and_zero_taxids(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(
        output_dir=tmp_path,
        taxids=["562", "0", "", 0],
        download_csdb=False,
    )

    assert mock_build.call_args.args[1] == ["562"]


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_reads_taxids_file(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    taxids_file = tmp_path / "taxids.txt"
    taxids_file.write_text(
        "# comment line\n"
        "562\n"
        "\n"
        "9606  # trailing comment\n"
        "511145\n",
        encoding="utf-8",
    )
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(
        output_dir=tmp_path,
        taxids_file=taxids_file,
        download_csdb=False,
    )

    assert mock_build.call_args.args[1] == ["511145", "562", "9606"]


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_merges_taxids_and_file(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    taxids_file = tmp_path / "taxids.txt"
    taxids_file.write_text("9606\n562\n", encoding="utf-8")
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(
        output_dir=tmp_path,
        taxids=["562", "511145"],
        taxids_file=taxids_file,
        download_csdb=False,
    )

    assert mock_build.call_args.args[1] == ["511145", "562", "9606"]


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_passes_archive_and_gene_set(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    archive = tmp_path / "csdb.tar.gz"
    archive.write_bytes(b"placeholder")
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(
        output_dir=tmp_path,
        taxids=["562"],
        csdb_archive=archive,
        download_csdb=False,
        gene_set="ribosomal",
    )

    _, kwargs = mock_build.call_args
    assert kwargs["csdb_archive"] == archive
    assert kwargs["download_csdb"] is False
    assert kwargs["gene_set"] == "ribosomal"


@patch("plasmidScreen.api.default_codon_usage_dir")
@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_defaults_output_dir(
    mock_build: MagicMock,
    mock_default_dir: MagicMock,
    tmp_path: Path,
) -> None:
    default_dir = tmp_path / "default_codon_usage"
    mock_default_dir.return_value = default_dir
    mock_build.return_value = _fake_result(default_dir)

    build_codon_database(taxids=["562"], download_csdb=False)

    mock_default_dir.assert_called_once()
    assert mock_build.call_args.args[0] == default_dir


@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_accepts_string_output_dir(
    mock_build: MagicMock, tmp_path: Path
) -> None:
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(output_dir=str(tmp_path), download_csdb=False)

    assert mock_build.call_args.args[0] == tmp_path


@pytest.mark.parametrize("gene_set", ["nuclear", "ribosomal", "mitochondrial", "plastid"])
@patch("plasmidScreen.api.build_codon_reference")
def test_build_codon_database_accepts_all_gene_sets(
    mock_build: MagicMock, gene_set: str, tmp_path: Path
) -> None:
    mock_build.return_value = _fake_result(tmp_path)

    build_codon_database(
        output_dir=tmp_path,
        taxids=["562"],
        download_csdb=False,
        gene_set=gene_set,  # type: ignore[arg-type]
    )

    assert mock_build.call_args.kwargs["gene_set"] == gene_set


# ---------------------------------------------------------------------------
# Integration: real CSDB March 2022 archive on disk
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_build_codon_database_real_csdb_imports_requested_taxids(
    real_csdb_build,
) -> None:
    output_dir, result = real_csdb_build

    assert sorted(result.taxids_requested) == sorted(REAL_CSDB_TAXIDS)
    assert sorted(result.taxids_added) == sorted(REAL_CSDB_TAXIDS)
    assert result.taxids_failed == []
    assert result.data_dir == output_dir
    assert (output_dir / CODON_TABLES_FILE).is_file()

    store = CodonUsageStore.load(output_dir)
    for taxid in REAL_CSDB_TAXIDS:
        assert store.has_codon_table(taxid), f"missing table for {taxid}"
        freqs = store.get_frequencies(taxid)
        assert freqs is not None
        assert len(freqs) >= 50  # most nuclear CSDB tables cover all sense codons
        assert all(0.0 <= v <= 1.0 for v in freqs.values())
        assert "TAA" not in freqs  # stops excluded by parser


@pytest.mark.integration
def test_build_codon_database_real_csdb_human_table_plausible(
    real_csdb_build,
) -> None:
    """Spot-check a few well-known human synonymous fractions."""
    output_dir, _result = real_csdb_build
    freqs = CodonUsageStore.load(output_dir).get_frequencies("9606")
    assert freqs is not None
    # Ala GCT fraction from the March 2022 CSDB release (non-uniform).
    assert freqs["GCT"] == pytest.approx(0.2627, abs=1e-3)
    assert freqs["GCG"] < freqs["GCT"]


@pytest.mark.integration
def test_build_codon_database_real_csdb_taxids_file(
    real_csdb_archive: Path, tmp_path: Path
) -> None:
    taxids_file = tmp_path / "hosts.txt"
    taxids_file.write_text(
        "# integration hosts\n"
        "9606\n"
        "511145\n",
        encoding="utf-8",
    )
    out = tmp_path / "codon_from_file"

    result = build_codon_database(
        output_dir=out,
        taxids_file=taxids_file,
        csdb_archive=real_csdb_archive,
        download_csdb=False,
    )

    assert result.taxids_requested == ["511145", "9606"]
    assert set(result.taxids_added) == {"511145", "9606"}
    assert result.taxids_failed == []
    store = CodonUsageStore.load(out)
    assert store.has_codon_table("9606")
    assert store.has_codon_table("511145")


@pytest.mark.integration
def test_build_codon_database_real_csdb_skips_existing(
    real_csdb_archive: Path, real_csdb_build
) -> None:
    output_dir, _first = real_csdb_build

    second = build_codon_database(
        output_dir=output_dir,
        taxids=list(REAL_CSDB_TAXIDS),
        csdb_archive=real_csdb_archive,
        download_csdb=False,
    )

    assert second.taxids_added == []
    assert sorted(second.taxids_skipped) == sorted(REAL_CSDB_TAXIDS)
    assert second.taxids_failed == []


@pytest.mark.integration
def test_build_codon_database_real_csdb_unknown_taxid_raises(
    real_csdb_archive: Path, tmp_path: Path
) -> None:
    """A taxid absent from CSDB (and with no lineage map) yields no tables."""
    with pytest.raises(RuntimeError, match="No codon tables written"):
        build_codon_database(
            output_dir=tmp_path / "codon_unknown",
            taxids=["999999999"],
            csdb_archive=real_csdb_archive,
            download_csdb=False,
        )