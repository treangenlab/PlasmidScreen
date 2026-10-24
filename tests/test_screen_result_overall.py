"""Tests for combined engineered decisions on ScreenResult."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.lib.models import compute_engineered_overall
from plasmidScreen.src.plasmidScreen import Workflow


def test_compute_engineered_overall_kmer_only() -> None:
    assert compute_engineered_overall(
        engineered_by_kmer_scan=True,
        engineered_by_codon_cai=None,
    )
    assert not compute_engineered_overall(
        engineered_by_kmer_scan=False,
        engineered_by_codon_cai=None,
    )


def test_compute_engineered_overall_codon_only() -> None:
    assert compute_engineered_overall(
        engineered_by_kmer_scan=False,
        engineered_by_codon_cai=True,
    )
    assert not compute_engineered_overall(
        engineered_by_kmer_scan=False,
        engineered_by_codon_cai=False,
    )


def test_compute_engineered_overall_either_signal() -> None:
    assert compute_engineered_overall(
        engineered_by_kmer_scan=True,
        engineered_by_codon_cai=False,
    )


DIAMOND_LINES = [
    "read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\t\n",
]


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
@patch("plasmidScreen.src.analyze_codon_usage.compute_cai", autospec=True)
def test_overall_engineered_by_codon_cai_threshold(
    mock_cai,
    mock_diamond,
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    mock_diamond.return_value = (DIAMOND_LINES, None)
    mock_cai.return_value = 0.9  # above threshold → engineered

    wf = Workflow(
        str(sample_fasta),
        1,
        None,
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=True,
        codon_cai_engineered_threshold=0.7,
        diamond_db=str(tmp_path / "diamond.dmnd"),
    )
    result = wf.run()
    per = {r.read_id: r for r in result.per_read}

    assert per["read_natural_1"].engineered_by_kmer_scan is False
    assert per["read_natural_1"].engineered_by_codon_cai is True
    assert per["read_natural_1"].engineered_overall is True
    assert per["read_natural_1"].overall_label == "Synthetic"
    assert per["read_synthetic_1"].engineered_overall is True
    assert result.overall_synthetic_count == 2
    assert result.codon_cai_engineered_threshold == 0.7
    assert result.codon_cai_engineered_read_ids == {"read_natural_1"}
    assert result.engineered_scan.engineered_read_ids == {"read_synthetic_1"}
    assert result.engineered_read_ids == {"read_natural_1", "read_synthetic_1"}
