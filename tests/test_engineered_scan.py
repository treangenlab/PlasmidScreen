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
    hit, max_eng, _target, _mapped, _cov = Workflow.parse_and_run(
        kmer_field, window_size=200, threshold=25
    )
    assert hit is expected
    assert isinstance(max_eng, int)


def test_parse_and_run_malformed_or_oversized_counts() -> None:
    assert Workflow.parse_and_run("", window_size=200, threshold=25)[0] is False
    assert Workflow.parse_and_run("562:66:extra", window_size=200, threshold=25)[0] is False
    huge = "32630:" + str(10**12)
    assert Workflow.parse_and_run(huge, window_size=200, threshold=25)[0] in (
        False,
        True,
    )


def test_parse_and_run_unclassified_zero_count_tokens() -> None:
    """Kraken unclassified reads often emit ``0:0``; must not crash or call engineered."""
    empty = (False, 0, 0, 0, 0.0)
    assert Workflow.parse_and_run("0:0", window_size=200, threshold=10) == empty
    assert Workflow.parse_and_run("0:0 0:0", window_size=200, threshold=10) == empty
    assert Workflow.parse_and_run("0:0\n", window_size=200, threshold=10) == empty
    # Zero-width tokens mixed with real evidence are ignored.
    hit, _max_eng, target, mapped, _cov = Workflow.parse_and_run(
        "562:10 0:0 32630:30", window_size=200, threshold=25
    )
    assert hit is True
    assert target == 30
    assert mapped == 40  # taxid 0 excluded from mapped total


def test_entry_logic_unclassified_0_0(tmp_path: Path) -> None:
    wf = Workflow(
        str(tmp_path / "reads.fa"),
        1,
        None,
        kraken_db=str(tmp_path / "db"),
        threads=1,
        kraken_raw_output=str(tmp_path / "k.out"),
        run_kraken=False,
        run_codon_usage=False,
    )
    (tmp_path / "reads.fa").write_text(">r1\nACGT\n")
    (tmp_path / "k.out").write_text("U\tr1\t0\t4\t0:0\n")
    kmer_max: dict = {}
    from plasmidScreen.lib.models import EngineeredScanResult

    result = EngineeredScanResult()
    t, n = wf.entry_logic("U\tr1\tunclassified (taxid 0)\t4\t0:0\n", kmer_max, result)
    assert (t, n) == (0, 0)
    assert result.natural_count == 1
    assert result.synthetic_count == 0
    assert result.labels[0].label == "Natural"


def test_workflow_completely_unmapped_read_0_0(tmp_path: Path) -> None:
    """End-to-end: a fully unclassified Kraken line with ``0:0`` must not crash."""
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">read_unmapped\nACGTACGTACGTACGT\n")
    kraken = tmp_path / "kraken.out"
    # Status U, taxid 0, empty LCA evidence ``0:0`` (breaks naive int/window paths).
    kraken.write_text("U\tread_unmapped\t0\t16\t0:0\n")
    report = tmp_path / "report.txt"

    wf = Workflow(
        str(fasta),
        1,
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(kraken),
        run_kraken=False,
        run_codon_usage=False,
    )
    scan, p_values, kmer_max_by_read = wf.scan_engineered_blocks_kraken()
    assert scan.synthetic_count == 0
    assert scan.natural_count == 1
    assert scan.labels[0].read_id == "read_unmapped"
    assert scan.labels[0].label == "Natural"
    assert len(p_values) == 1
    assert kmer_max_by_read["read_unmapped"] == 0

    result = wf.run()
    per = {r.read_id: r for r in result.per_read}
    assert per["read_unmapped"].engineered_by_kmer_scan is False
    assert per["read_unmapped"].engineered_overall is False
    assert per["read_unmapped"].overall_label == "Natural"
    assert "Natural\tread_unmapped" in report.read_text()


def test_parse_and_run_numba_unbox_fallback(monkeypatch) -> None:
    from plasmidScreen.src import plasmidScreen as ps

    def _raise_type_error(*_args, **_kwargs):
        raise TypeError(
            "can't unbox array from PyObject into native value.  "
            "The object maybe of a different type"
        )

    monkeypatch.setattr(ps, "fast_window_logic", _raise_type_error)
    hit, max_eng, *_rest = Workflow.parse_and_run(
        "32630:66", window_size=200, threshold=25
    )
    assert hit is True
    assert max_eng == 25


def test_scan_engineered_blocks_integration(
    sample_fasta: Path,
    sample_kraken_out: Path,
    tmp_path: Path,
) -> None:
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        1,
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        run_kraken=False,
        run_codon_usage=False,
    )
    scan, p_values, kmer_max_by_read = wf.scan_engineered_blocks_kraken()
    assert scan.synthetic_count == 1
    assert scan.natural_count == 2  # natural read + unclassified (no 32630 signal)
    labels = {r.read_id: r.label for r in scan.labels}
    assert labels["read_natural_1"] == "Natural"
    assert labels["read_synthetic_1"] == "Synthetic"
    assert labels["read_unclassified"] == "Natural"
    assert len(p_values) == 3
    assert set(kmer_max_by_read) >= {"read_natural_1", "read_synthetic_1"}

    # Report is written by run(), not scan alone.
    result = wf.run()
    text = result.engineered_report_path.read_text()
    assert "Natural\tread_natural_1" in text
    assert "Engineered\tread_synthetic_1\tengineered_kmer_scan" in text


def test_insignificant_pvalue_reassigns_window_synthetic_as_natural(
    sample_fasta: Path,
    tmp_path: Path,
) -> None:
    """Window threshold met but p > alpha → Natural for engineered call."""
    length = 100
    k = 35
    n_kmers = length - k + 1
    # One engineered k-mer: meets threshold=1 but is insignificant under bg=0.001.
    kraken = tmp_path / "kraken_weak.out"
    kraken.write_text(
        f"C\tread_synthetic_1\t562\t{length}\t0:1\t"
        f"32630:1 562:{max(n_kmers - 1, 1)}\n"
        f"C\tread_natural_1\t562\t{length}\t0:1\t562:{max(n_kmers, 1)}\n"
    )
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        1,
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(kraken),
        run_kraken=False,
        run_codon_usage=False,
        engineered_kmer_threshold=1,
        background_rate=0.001,
        pvalue_alpha=0.05,
    )
    result = wf.run()
    per = {r.read_id: r for r in result.per_read}
    assert per["read_synthetic_1"].p_value is not None
    assert per["read_synthetic_1"].p_value > 0.05
    assert per["read_synthetic_1"].engineered_by_kmer_scan is False
    assert per["read_synthetic_1"].engineered_overall is False
    assert per["read_synthetic_1"].overall_label == "Natural"
    assert per["read_synthetic_1"].kmer_label == "Natural"
    text = report.read_text()
    assert "Natural\tread_synthetic_1" in text
