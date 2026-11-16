"""CLI smoke / validation tests via Typer's CliRunner."""
from __future__ import annotations

import re
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from plasmidScreen.cli import app
from plasmidScreen.lib.diamond_host_taxonomy import diamond_technology_flags
from plasmidScreen.lib.types import GeneSet

runner = CliRunner()

# Rich/Typer may bold option names by inserting ANSI codes mid-token
# (e.g. ``--diamond-db`` → bold ``-`` + ``-diamond`` + ``-db``).
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def _cli_text(result) -> str:
    raw = (result.output or "") + (getattr(result, "stderr", None) or "")
    return _ANSI_RE.sub("", raw).lower()


def test_cli_help_boots() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "screen" in result.output
    assert "build" in result.output
    assert "visualize" in result.output
    assert "prefetch-genbank" in result.output


def test_screen_visualize_requires_reference_db() -> None:
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--visualize",
        ],
    )
    assert result.exit_code != 0
    assert "reference-db" in _cli_text(result)


def test_screen_map_hits_requires_reference_db() -> None:
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--map-hits",
        ],
    )
    assert result.exit_code != 0
    assert "reference-db" in _cli_text(result)


def test_screen_no_run_kraken_requires_path() -> None:
    result = runner.invoke(
        app,
        ["screen", "reads.fa", "report.txt", "kraken_db", "--no-run-kraken", "--skip-codon-usage"],
    )
    assert result.exit_code != 0
    assert "kraken-output-path" in _cli_text(result) or result.exit_code == 2


def test_screen_codon_requires_diamond_db() -> None:
    result = runner.invoke(
        app,
        ["screen", "reads.fa", "report.txt", "kraken_db"],
    )
    assert result.exit_code != 0
    assert "diamond-db" in _cli_text(result)


@patch("plasmidScreen.cli.run_screen")
def test_screen_skip_codon_usage_invokes_run_screen(mock_run: MagicMock) -> None:
    mock_run.return_value = MagicMock(overall_synthetic_count=0, per_read=[])
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--threshold",
            "15",
            "--bg",
            "0.001",
        ],
    )
    assert result.exit_code == 0, result.output
    mock_run.assert_called_once()
    kwargs = mock_run.call_args.kwargs
    assert kwargs["run_codon_usage"] is False
    assert kwargs["engineered_kmer_threshold"] == 15
    assert kwargs["background_rate"] == 0.001
    assert kwargs["kmer_size"] == 21
    assert kwargs["minimizer_size"] == 17
    assert kwargs["seq_type"] == 1  # CLI default: short read
    assert kwargs["quiet_mode"] is False


@patch("plasmidScreen.cli.run_screen")
def test_screen_passes_kmer_and_minimizer_size(mock_run: MagicMock) -> None:
    mock_run.return_value = MagicMock(overall_synthetic_count=0, per_read=[])
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--kmer-size",
            "35",
            "--minimizer-size",
            "31",
        ],
    )
    assert result.exit_code == 0, result.output
    kwargs = mock_run.call_args.kwargs
    assert kwargs["kmer_size"] == 35
    assert kwargs["minimizer_size"] == 31


@patch("plasmidScreen.cli.map_and_visualize_engineered")
@patch("plasmidScreen.cli.run_screen")
def test_screen_visualize_invokes_map_and_visualize(
    mock_run: MagicMock,
    mock_map_viz: MagicMock,
) -> None:
    screen = MagicMock(overall_synthetic_count=1, per_read=[MagicMock()])
    mock_run.return_value = screen
    mock_map_viz.return_value = (
        screen,
        MagicMock(pileup_paths=[Path("a.png")], output_dir=Path("viz_out")),
    )
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--visualize",
            "--reference-db",
            "plasmids.fa",
            "--viz-dir",
            "viz_out",
        ],
    )
    assert result.exit_code == 0, result.output
    mock_run.assert_called_once()
    mock_map_viz.assert_called_once()
    kwargs = mock_map_viz.call_args.kwargs
    assert kwargs["query_fasta"] == "reads.fa"
    assert kwargs["reference_db"] == "plasmids.fa"
    assert kwargs["viz_dir"] == "viz_out"


@patch("plasmidScreen.cli.map_engineered_reads")
@patch("plasmidScreen.cli.run_screen")
def test_screen_map_hits_without_visualize(
    mock_run: MagicMock,
    mock_map: MagicMock,
) -> None:
    screen = MagicMock(
        overall_synthetic_count=1,
        per_read=[MagicMock()],
        top_engineered_references=("ref1",),
    )
    mock_run.return_value = screen
    mock_map.return_value = screen
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--map-hits",
            "--reference-db",
            "plasmids.fa",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.output
    mock_map.assert_called_once()
    assert mock_map.call_args.kwargs["reference_db"] == "plasmids.fa"


@patch("plasmidScreen.cli.map_and_visualize_engineered")
@patch("plasmidScreen.cli.screen_result_from_engineered_report")
def test_visualize_command_invokes_pipeline(
    mock_from_report: MagicMock,
    mock_map_viz: MagicMock,
) -> None:
    rebuilt = MagicMock(overall_synthetic_count=2, per_read=[1, 2])
    mock_from_report.return_value = rebuilt
    mock_map_viz.return_value = (
        rebuilt,
        MagicMock(pileup_paths=[Path("p.png")], output_dir=Path("viz_out")),
    )
    result = runner.invoke(
        app,
        [
            "visualize",
            "reads.fa",
            "report.txt",
            "plasmids.fa",
            "viz_out",
            "--top-n",
            "3",
        ],
    )
    assert result.exit_code == 0, result.output
    mock_from_report.assert_called_once()
    mock_map_viz.assert_called_once()
    kwargs = mock_map_viz.call_args.kwargs
    assert kwargs["top_n"] == 3
    assert kwargs["reference_db"] == "plasmids.fa"
    # Default: cache-only; --fetch-genbank opts into on-demand NCBI.
    assert kwargs.get("fetch_genbank") is False


@patch("plasmidScreen.cli.prefetch_genbank_from_fasta")
def test_prefetch_genbank_command(mock_prefetch: MagicMock, tmp_path: Path) -> None:
    mock_prefetch.return_value = MagicMock(
        ids_requested=2,
        ids_cached=2,
        ids_missing=0,
        ids_fetched_this_run=1,
        cache_dir=tmp_path,
    )
    result = runner.invoke(
        app,
        [
            "prefetch-genbank",
            "plasmids.fa",
            "--cache-dir",
            str(tmp_path),
            "--ncbi-api-key",
            "KEY",
        ],
    )
    assert result.exit_code == 0, result.output
    mock_prefetch.assert_called_once()
    args, kwargs = mock_prefetch.call_args
    assert args[0] == "plasmids.fa"
    assert kwargs["cache_dir"] == str(tmp_path)
    assert kwargs["api_key"] == "KEY"


@patch("plasmidScreen.cli.run_screen")
def test_screen_seq_type_long_read(mock_run: MagicMock) -> None:
    mock_run.return_value = MagicMock(overall_synthetic_count=0, per_read=[])
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--seq-type",
            "0",
        ],
    )
    assert result.exit_code == 0, result.output
    assert mock_run.call_args.kwargs["seq_type"] == 0


@pytest.mark.parametrize(
    "seq_type,required_flags,forbidden_flags",
    [
        (
            0,
            ["--frameshift", "15", "--range-culling", "--fast", "-k", "1"],
            ["--more-sensitive", "--top"],
        ),
        (
            1,
            ["--more-sensitive", "--top", "5", "--min-orf", "10"],
            ["--frameshift", "--range-culling", "--fast"],
        ),
    ],
)
@patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run")
def test_cli_seq_type_selects_diamond_params(
    mock_diamond_run: MagicMock,
    seq_type: int,
    required_flags: list[str],
    forbidden_flags: list[str],
    sample_fasta: Path,
    sample_kraken_out: Path,
    codon_store_dir: Path,
    tmp_path: Path,
) -> None:
    """CLI --seq-type must drive the matching DIAMOND blastx flag set."""
    diamond_db = tmp_path / "protein.dmnd"
    diamond_db.write_bytes(b"")
    report = tmp_path / "report.txt"

    # One Natural-read hit so codon CAI can finish after DIAMOND returns.
    mock_diamond_run.return_value = MagicMock(
        returncode=0,
        stdout=(
            "read_natural_1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\tEscherichia coli\n"
        ),
        stderr="",
    )

    result = runner.invoke(
        app,
        [
            "screen",
            str(sample_fasta),
            str(report),
            str(tmp_path / "unused_kraken_db"),
            "--no-run-kraken",
            "--kraken-output-path",
            str(sample_kraken_out),
            "--diamond-db",
            str(diamond_db),
            "--codon-usage-dir",
            str(codon_store_dir),
            "--seq-type",
            str(seq_type),
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.output
    mock_diamond_run.assert_called_once()
    cmd = mock_diamond_run.call_args[0][0]
    assert cmd[:2] == ["diamond", "blastx"]

    for flag in required_flags:
        assert flag in cmd, f"expected {flag!r} in DIAMOND cmd for seq_type={seq_type}: {cmd}"
    for flag in forbidden_flags:
        assert flag not in cmd, f"did not expect {flag!r} for seq_type={seq_type}: {cmd}"

    tech = diamond_technology_flags(seq_type)
    assert " ".join(tech) in " ".join(cmd)


@patch("plasmidScreen.cli.build_codon_reference")
def test_build_parses_taxids(mock_build: MagicMock, tmp_path: Path) -> None:
    mock_build.return_value = MagicMock(
        taxids_added=["562"],
        taxids_skipped=[],
        taxids_failed=[],
        data_dir=tmp_path,
    )
    taxids_file = tmp_path / "taxids.txt"
    taxids_file.write_text("9606\n# comment\n511145\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "build",
            "--output-dir",
            str(tmp_path),
            "--taxids",
            "10090",
            "--taxids-file",
            str(taxids_file),
            "--no-download-csdb",
            "--csdb-archive",
            str(tmp_path / "missing.tar.gz"),
            "--gene-set",
            "nuclear",
        ],
    )
    # build_codon_reference is mocked, so missing archive is fine
    assert result.exit_code == 0, result.output
    mock_build.assert_called_once()
    args, kwargs = mock_build.call_args
    assert set(args[1]) == {"10090", "9606", "511145"}
    assert kwargs["gene_set"] == GeneSet.NUCLEAR
    assert kwargs["download_csdb"] is False
