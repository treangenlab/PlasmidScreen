"""End-to-end workflow tests (no Kraken subprocess)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from plasmidScreen.src.plasmidScreen import Workflow


DIAMOND_LINES = [
    # outfmt 6: qseqid sseqid qstart qend pident length evalue bitscore staxids sscinames
    "read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\tEscherichia coli",
]


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_workflow_run_without_kraken(
    mock_diamond,
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    mock_diamond.return_value = (DIAMOND_LINES, None)
    report = tmp_path / "report.txt"
    codon_out = tmp_path / "codon.tsv"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        write_engineered_report=True,
        codon_usage_output_path=str(codon_out),
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=True,
        diamond_db=str(tmp_path / "diamond.dmnd"),
    )
    result = wf.run()

    mock_diamond.assert_called_once()
    assert result.engineered_scan.synthetic_count == 1
    assert len(result.codon_adaptation) == 1
    assert result.codon_adaptation[0].read_id == "read_natural_1"
    assert len(result.per_read) >= 2
    per = {r.read_id: r for r in result.per_read}
    assert per["read_synthetic_1"].engineered_by_kmer_scan is True
    assert per["read_synthetic_1"].engineered_overall is True
    assert per["read_synthetic_1"].overall_label == "Synthetic"
    assert per["read_synthetic_1"].engineered_methods == ["engineered_kmer_scan"]
    assert per["read_natural_1"].engineered_by_kmer_scan is False
    assert per["read_natural_1"].engineered_overall is False
    assert per["read_natural_1"].overall_label == "Natural"
    assert result.overall_synthetic_count == 1
    assert result.overall_natural_count == 2
    assert result.engineered_read_ids == {"read_synthetic_1"}
    assert result.codon_usage_report_path is not None
    assert result.codon_usage_report_path.exists()
    assert result.engineered_report_path == report


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_workflow_does_not_write_codon_tsv_by_default(
    mock_diamond,
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    mock_diamond.return_value = (DIAMOND_LINES, None)
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        write_engineered_report=True,
        codon_usage_output_path=None,
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=True,
        diamond_db=str(tmp_path / "diamond.dmnd"),
    )
    result = wf.run()
    assert len(result.codon_adaptation) == 1
    assert result.codon_usage_report_path is None


def test_workflow_skip_codon_usage(
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        write_engineered_report=True,
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=False,
    )
    result = wf.run()
    assert result.codon_adaptation == []
    assert result.codon_usage_report_path is None


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_workflow_in_memory_no_report_files(
    mock_diamond,
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    mock_diamond.return_value = (DIAMOND_LINES, None)
    wf = Workflow(
        str(sample_fasta),
        None,
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        write_engineered_report=False,
        codon_usage_output_path=None,
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=True,
        diamond_db=str(tmp_path / "diamond.dmnd"),
    )
    result = wf.run()
    assert result.engineered_report_path is None
    assert result.codon_usage_report_path is None
    assert len(result.codon_adaptation) == 1


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_workflow_loads_precomputed_diamond_tsv(
    mock_diamond,
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    diamond_tsv = tmp_path / "diamond.tsv"
    diamond_tsv.write_text(
        "read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\t\n",
        encoding="utf-8",
    )
    mock_diamond.return_value = (
        diamond_tsv.read_text(encoding="utf-8").splitlines(),
        diamond_tsv,
    )
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        write_engineered_report=True,
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=True,
        run_diamond=False,
        diamond_output_path=str(diamond_tsv),
    )
    result = wf.run()
    mock_diamond.assert_called_once()
    assert result.diamond_output_path == diamond_tsv
    assert len(result.codon_adaptation) == 1


@patch("plasmidScreen.src.plasmidScreen.subprocess.run")
def test_run_kraken_enabled_calls_subprocess(
    mock_subprocess,
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    report = tmp_path / "report.txt"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=None,
        write_engineered_report=True,
        codon_usage_dir=str(codon_store_dir),
        run_kraken=True,
        run_codon_usage=False,
    )
    wf.run()
    mock_subprocess.assert_called_once()
    args = mock_subprocess.call_args[0][0]
    # In-memory by default: omit --output so Kraken2 writes classifications to stdout
    assert "--output" not in args
    assert "--report" not in args
