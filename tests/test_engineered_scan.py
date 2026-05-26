"""Tests for engineered k-mer block detection."""
from __future__ import annotations

from pathlib import Path

import pytest

from plasmidScreen.src.plasmidScreen import Workflow


@pytest.mark.parametrize(
    "kmer_field,expected",
    [
        ("562:66", False),
        ("32630:66", True),
    ],
)
def test_parse_and_run_synthetic_detection(kmer_field: str, expected: bool) -> None:
    result = Workflow.parse_and_run(kmer_field, window_size=200, threshold=25)
    assert result is expected


def test_scan_engineered_blocks_integration(
    sample_fasta: Path,
    sample_kraken_out: Path,
    tmp_path: Path,
) -> None:
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        run_kraken=False,
        run_codon_usage=False,
    )
    scan = wf.scan_engineered_blocks_kraken()
    assert scan.synthetic_count == 1
    assert scan.natural_count == 2  # natural read + unclassified (no 32630 signal)
    labels = {r.read_id: r.label for r in scan.labels}
    assert labels["read_natural_1"] == "Natural"
    assert labels["read_synthetic_1"] == "Synthetic"
    assert labels["read_unclassified"] == "Natural"
    text = report.read_text()
    assert "Natural\tread_natural_1" in text
    assert "Synthetic\tread_synthetic_1" in text
