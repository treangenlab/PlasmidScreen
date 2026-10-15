"""Tests for codon adaptation analysis."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from plasmidScreen.lib.exceptions import MissingCodonReferenceError
from plasmidScreen.lib.models import CodonAdaptationRead
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
    codon_adaptation_to_tsv_lines,
    compute_cai,
    parse_kraken_file,
    write_codon_adaptation_tsv,
    write_subset_reads,
)

DIAMOND_LINES = [
    "read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\tEscherichia coli",
]


def test_parse_kraken_file(fixtures_kraken_out: Path) -> None:
    data = parse_kraken_file(str(fixtures_kraken_out))
    assert "read_natural_1" in data
    status, taxid, length, k_info = data["read_natural_1"]
    assert status == "C"
    assert taxid == "562"
    assert length == 100
    assert "562:" in k_info


def test_compute_cai_perfect_match() -> None:
    weights = {"ATG": 1.0, "AAA": 1.0}
    assert compute_cai("ATGAAA", weights) == pytest.approx(1.0)


def test_write_subset_reads_filters_ids(sample_fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "natural_only.fasta"
    n = write_subset_reads(sample_fasta, {"read_natural_1"}, out)
    assert n == 1
    text = out.read_text()
    assert ">read_natural_1" in text
    assert "read_synthetic_1" not in text


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_analyze_codon_adaptation_natural_reads_only(
    mock_diamond,
    sample_fasta: Path,
    codon_store_dir: Path,
) -> None:
    mock_diamond.return_value = (DIAMOND_LINES, None)
    results, _path = analyze_codon_adaptation(
        sample_fasta,
        diamond_db="db.dmnd",
        codon_usage_dir=codon_store_dir,
        include_read_ids={"read_natural_1"},
    )
    assert len(results) == 1
    r = results[0]
    assert r.read_id == "read_natural_1"
    assert r.host_taxid == "562"
    assert r.reference_taxid == "562"
    assert r.cai_vs_host is not None
    assert 0.0 <= r.cai_vs_host <= 1.0
    # DIAMOND should be invoked on a temp subset, not the full input FASTA
    diamond_query = Path(mock_diamond.call_args[0][0])
    assert diamond_query != sample_fasta
    assert diamond_query.exists() is False  # cleaned up in finally


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_analyze_skips_reads_without_diamond_hits(
    mock_diamond,
    sample_fasta: Path,
    codon_store_dir: Path,
) -> None:
    mock_diamond.return_value = ([], None)
    results, _path = analyze_codon_adaptation(
        sample_fasta,
        diamond_db="db.dmnd",
        codon_usage_dir=codon_store_dir,
        include_read_ids={"read_unclassified"},
    )
    assert len(results) == 0


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_missing_reference_raises(
    mock_diamond,
    sample_fasta: Path,
    codon_store_dir: Path,
) -> None:
    mock_diamond.return_value = (
        ["read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t99999\t\n"],
        None,
    )
    with pytest.raises(MissingCodonReferenceError) as exc:
        analyze_codon_adaptation(
            sample_fasta,
            diamond_db="db.dmnd",
            codon_usage_dir=codon_store_dir,
            include_read_ids={"read_natural_1"},
        )
    assert "99999" in exc.value.missing_host_taxids


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_write_tsv_roundtrip(
    mock_diamond,
    sample_fasta: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    mock_diamond.return_value = (DIAMOND_LINES, None)
    out = tmp_path / "codon.tsv"
    path, results = write_codon_adaptation_tsv(
        out,
        sample_fasta,
        diamond_db="db.dmnd",
        codon_usage_dir=codon_store_dir,
        include_read_ids={"read_natural_1"},
    )
    assert path == str(out)
    assert len(results) == 1
    lines = out.read_text().strip().splitlines()
    assert lines[0].startswith("Read_ID\t")
    assert "read_natural_1" in lines[1]


def test_tsv_lines_match_results() -> None:
    row = CodonAdaptationRead(
        read_id="r1",
        cds_strand="+",
        cds_start=0,
        cds_end=12,
        host_taxid="562",
        reference_taxid="562",
        cds_len_bp=12,
        cai_vs_host=0.5,
    )
    lines = codon_adaptation_to_tsv_lines([row])
    assert len(lines) == 2
    assert "0.5000" in lines[1]


def test_write_results_tsv_with_cai_threshold_columns(tmp_path: Path) -> None:
    from plasmidScreen.src.analyze_codon_usage import write_codon_adaptation_results_tsv

    rows = [
        CodonAdaptationRead(
            read_id="high",
            cds_strand="+",
            cds_start=0,
            cds_end=12,
            host_taxid="562",
            reference_taxid="562",
            cds_len_bp=12,
            cai_vs_host=0.9,
        ),
        CodonAdaptationRead(
            read_id="low",
            cds_strand="+",
            cds_start=0,
            cds_end=12,
            host_taxid="562",
            reference_taxid="562",
            cds_len_bp=12,
            cai_vs_host=0.2,
        ),
    ]
    out = tmp_path / "codon_thr.tsv"
    write_codon_adaptation_results_tsv(out, rows, cai_engineered_threshold=0.6)
    text = out.read_text()
    header = text.splitlines()[0]
    assert "Codon_CAI_Threshold" in header
    assert "Engineered_By_Codon_CAI" in header
    assert "high\t" in text and "True" in text.splitlines()[1]
    assert "low\t" in text and "False" in text.splitlines()[2]
