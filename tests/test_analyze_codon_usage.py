"""Tests for codon adaptation analysis."""
from __future__ import annotations

from pathlib import Path

import pytest

from plasmidScreen.lib.exceptions import MissingCodonReferenceError
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
    codon_adaptation_to_tsv_lines,
    compute_cai,
    expand_kmer_taxids,
    parse_kraken_file,
    write_codon_adaptation_tsv,
)


def test_parse_kraken_file(fixtures_kraken_out: Path) -> None:
    data = parse_kraken_file(str(fixtures_kraken_out))
    assert "read_natural_1" in data
    status, taxid, length, k_info = data["read_natural_1"]
    assert status == "C"
    assert taxid == "562"
    assert length == 100
    assert "562:" in k_info


def test_expand_kmer_taxids_pads_short_run() -> None:
    kmers = expand_kmer_taxids("562:5", read_length=100, k=35)
    assert len(kmers) == 66
    assert kmers[0] == "562"


def test_compute_cai_perfect_match() -> None:
    weights = {"ATG": 1.0, "AAA": 1.0}
    assert compute_cai("ATGAAA", weights) == pytest.approx(1.0)


def test_analyze_codon_adaptation_natural_reads_only(
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
) -> None:
    results = analyze_codon_adaptation(
        sample_fasta,
        sample_kraken_out,
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


def test_analyze_skips_unclassified(
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
) -> None:
    results = analyze_codon_adaptation(
        sample_fasta,
        sample_kraken_out,
        codon_usage_dir=codon_store_dir,
        include_read_ids={"read_unclassified"},
    )
    assert results == []


def test_missing_reference_raises(
    sample_fasta: Path,
    tmp_path: Path,
    codon_store_dir: Path,
) -> None:
    kraken = tmp_path / "kraken_missing_host.out"
    kraken.write_text("C\tread_natural_1\t99999\t100\t0:1\t99999:66\n")
    with pytest.raises(MissingCodonReferenceError) as exc:
        analyze_codon_adaptation(
            sample_fasta,
            kraken,
            codon_usage_dir=codon_store_dir,
            include_read_ids={"read_natural_1"},
        )
    assert "99999" in exc.value.missing_host_taxids


def test_write_tsv_roundtrip(
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    out = tmp_path / "codon.tsv"
    path, results = write_codon_adaptation_tsv(
        sample_fasta,
        sample_kraken_out,
        out,
        codon_usage_dir=codon_store_dir,
        include_read_ids={"read_natural_1"},
    )
    assert path == str(out)
    assert len(results) == 1
    lines = out.read_text().strip().splitlines()
    assert lines[0].startswith("Read_ID\t")
    assert "read_natural_1" in lines[1]


def test_tsv_lines_match_results() -> None:
    from plasmidScreen.lib.models import CodonAdaptationResult

    row = CodonAdaptationResult(
        read_id="r1",
        overall_taxid="562",
        cds_strand="+",
        cds_start=0,
        cds_end=12,
        cds_taxid="562",
        host_taxid="562",
        reference_taxid="562",
        cds_len_bp=12,
        cai_vs_host=0.5,
    )
    lines = codon_adaptation_to_tsv_lines([row])
    assert len(lines) == 2
    assert "0.5000" in lines[1]
