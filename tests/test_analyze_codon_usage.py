"""Tests for codon adaptation analysis."""
from __future__ import annotations

import math
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
    parse_kraken_lines,
    write_codon_adaptation_tsv,
    write_subset_reads,
)
from plasmidScreen.src.codon_usage.codon_usage_db import (
    CodonUsageStore,
    cai_weights_from_frequencies,
)

DIAMOND_LINES = [
    "read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\tEscherichia coli",
]


def _biased_lys_frequencies() -> dict[str, float]:
    """Prefer AAA over AAG for Lys; other AAs uniform so ATG stays weight 1."""
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    by_aa: dict[str, list[str]] = {}
    for codon, aa in forward.items():
        if aa == "*":
            continue
        by_aa.setdefault(aa, []).append(codon)
    freqs: dict[str, float] = {}
    for aa, codons in by_aa.items():
        if aa == "K":
            freqs["AAA"] = 1.0
            freqs["AAG"] = 0.1
            continue
        w = 1.0 / len(codons)
        for c in codons:
            freqs[c] = w
    return freqs


def test_parse_kraken_file(fixtures_kraken_out: Path) -> None:
    data = parse_kraken_file(str(fixtures_kraken_out))
    assert "read_natural_1" in data
    status, taxid, length, k_info = data["read_natural_1"]
    assert status == "C"
    assert taxid == "562"
    assert length == 100
    assert "562:" in k_info


def test_parse_kraken_lines_skips_malformed() -> None:
    data = parse_kraken_lines(
        [
            "C\tok\t562\t100\t562:10\n",
            "too\tshort\n",
            "C\tbadlen\t562\tnot_int\t562:1\n",
        ]
    )
    assert set(data) == {"ok"}


def test_compute_cai_perfect_match() -> None:
    weights = {"ATG": 1.0, "AAA": 1.0}
    assert compute_cai("ATGAAA", weights) == pytest.approx(1.0)


def test_compute_cai_geometric_mean() -> None:
    weights = {"AAA": 1.0, "AAG": 0.25}
    # geometric mean of 1.0 and 0.25 = 0.5
    assert compute_cai("AAAAAG", weights) == pytest.approx(math.sqrt(1.0 * 0.25))


def test_compute_cai_skips_unknown_and_zero_weight_codons() -> None:
    weights = {"AAA": 1.0, "NNN": 0.0}
    # NNN ignored (zero); trailing incomplete "AA" ignored; only AAA counts
    assert compute_cai("AAANNNAA", weights) == pytest.approx(1.0)


def test_compute_cai_empty_or_unweighted_is_zero() -> None:
    assert compute_cai("", {"AAA": 1.0}) == 0.0
    assert compute_cai("ATG", {}) == 0.0
    assert compute_cai("NNN", {"AAA": 1.0}) == 0.0


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
        CodonAdaptationRead(
            read_id="missing",
            cds_strand="+",
            cds_start=0,
            cds_end=12,
            host_taxid="562",
            reference_taxid="562",
            cds_len_bp=12,
            cai_vs_host=None,
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
    # CAI missing → engineered flag NA (threshold column still present)
    assert text.splitlines()[3].endswith("\t0.6\tNA")


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_analyze_cai_higher_for_preferred_codons(
    mock_diamond,
    tmp_path: Path,
) -> None:
    """With biased Lys weights, AAA-rich CDS should score higher CAI than AAG-rich.

    Conftest ``codon_store`` is uniform (all synonym weights 1.0), so this test
    builds a small in-memory store with Lys bias instead.
    """
    freqs = _biased_lys_frequencies()
    weights = cai_weights_from_frequencies(freqs)
    assert weights["AAA"] == pytest.approx(1.0)
    assert weights["AAG"] == pytest.approx(0.1)

    store = CodonUsageStore.writable(tmp_path / "biased_store")
    store.set_codon_table("562", freqs, source="test")

    preferred = "ATG" + ("AAA" * 20)
    rare = "ATG" + ("AAG" * 20)
    fasta = tmp_path / "bias.fa"
    fasta.write_text(
        f">read_pref\n{preferred}\n>read_rare\n{rare}\n",
        encoding="utf-8",
    )
    end = len(preferred)
    mock_diamond.return_value = (
        [
            f"read_pref\tref\t1\t{end}\t99.0\t{end // 3}\t1e-5\t200\t562\t\n",
            f"read_rare\tref\t1\t{end}\t99.0\t{end // 3}\t1e-5\t200\t562\t\n",
        ],
        None,
    )
    results, _ = analyze_codon_adaptation(
        fasta,
        diamond_db="db.dmnd",
        codon_usage_store=store,
        include_read_ids={"read_pref", "read_rare"},
    )
    by_id = {r.read_id: r for r in results.labels}
    assert by_id["read_pref"].cai_vs_host > by_id["read_rare"].cai_vs_host
    assert by_id["read_pref"].cai_vs_host == pytest.approx(1.0)
    assert by_id["read_rare"].cai_vs_host == pytest.approx(
        math.exp((math.log(1.0) + 20 * math.log(0.1)) / 21)
    )


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_analyze_selects_longest_orf_interval(
    mock_diamond,
    sample_fasta: Path,
    codon_store_dir: Path,
) -> None:
    # Gap must be > merge_gap_bp (30): short ends at 20, long starts at 62
    mock_diamond.return_value = (
        [
            "read_natural_1\tshort\t1\t20\t99.0\t7\t1e-5\t50\t562\t\n",
            "read_natural_1\tlong\t62\t99\t99.0\t13\t1e-5\t200\t562\t\n",
        ],
        None,
    )
    results, _ = analyze_codon_adaptation(
        sample_fasta,
        diamond_db="db.dmnd",
        codon_usage_dir=codon_store_dir,
        include_read_ids={"read_natural_1"},
    )
    assert len(results) == 1
    assert results[0].cds_start == 61  # 0-based from qstart 62
    assert results[0].cds_end == 99
    assert results[0].cds_len_bp == 38


@patch("plasmidScreen.src.analyze_codon_usage.resolve_diamond_lines", autospec=True)
def test_analyze_empty_include_set_writes_no_query(
    mock_diamond,
    sample_fasta: Path,
    codon_store_dir: Path,
) -> None:
    results, path = analyze_codon_adaptation(
        sample_fasta,
        diamond_db="db.dmnd",
        codon_usage_dir=codon_store_dir,
        include_read_ids=set(),
    )
    assert len(results) == 0
    assert path is None
    mock_diamond.assert_not_called()
