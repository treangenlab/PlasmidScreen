"""End-to-end workflow tests (no Kraken subprocess)."""
from __future__ import annotations

from pathlib import Path

from plasmidScreen.src.plasmidScreen import Workflow


def test_workflow_run_without_kraken(
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    report = tmp_path / "report.txt"
    codon_out = tmp_path / "codon.tsv"
    wf = Workflow(
        str(sample_fasta),
        str(report),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=str(sample_kraken_out),
        codon_usage_output_path=str(codon_out),
        codon_usage_dir=str(codon_store_dir),
        run_kraken=False,
        run_codon_usage=True,
    )
    result = wf.run()

    assert result.engineered_scan.synthetic_count == 1
    assert len(result.codon_adaptation) == 1
    assert result.codon_adaptation[0].read_id == "read_natural_1"
    assert len(result.per_read) >= 2
    per = {r.read_id: r for r in result.per_read}
    assert per["read_synthetic_1"].engineered_by_kmer_scan is True
    assert per["read_synthetic_1"].engineered_methods == ["engineered_kmer_scan"]
    assert per["read_natural_1"].engineered_by_kmer_scan is False
    assert result.codon_usage_report_path is not None
    assert result.codon_usage_report_path.exists()
    assert result.engineered_report_path == report
