"""Quiet-mode logging and input/error-path coverage."""
from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from plasmidScreen.cli import app
from plasmidScreen.lib.exceptions import InputError, KrakenError
from plasmidScreen.lib.logging_util import configure_logging
from plasmidScreen.lib.seqio_util import open_text_maybe_gzip
from plasmidScreen.src.plasmidScreen import Workflow

runner = CliRunner()


def test_configure_logging_quiet_raises_level() -> None:
    configure_logging(quiet=False, force=True)
    assert logging.getLogger().level == logging.INFO
    configure_logging(quiet=True, force=True)
    assert logging.getLogger().level == logging.WARNING


def test_open_text_maybe_gzip_missing_file(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="not found"):
        with open_text_maybe_gzip(tmp_path / "missing.fa"):
            pass


def test_open_text_maybe_gzip_bad_gzip(tmp_path: Path) -> None:
    bad = tmp_path / "bad.fa.gz"
    bad.write_bytes(b"not-gzip-data")
    with pytest.raises(InputError, match="Invalid gzip"):
        with open_text_maybe_gzip(bad) as handle:
            handle.read()


def test_workflow_missing_fasta(tmp_path: Path) -> None:
    wf = Workflow(
        str(tmp_path / "nope.fa"),
        1,
        str(tmp_path / "out.txt"),
        kraken_db=str(tmp_path / "db"),
        threads=1,
        kraken_raw_output=str(tmp_path / "k.out"),
        run_kraken=False,
        run_codon_usage=False,
    )
    with pytest.raises(InputError, match="sequence file not found"):
        wf.run()


def test_workflow_empty_fasta(tmp_path: Path) -> None:
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    wf = Workflow(
        str(empty),
        1,
        str(tmp_path / "out.txt"),
        kraken_db=str(tmp_path / "db"),
        threads=1,
        kraken_raw_output=str(tmp_path / "k.out"),
        run_kraken=False,
        run_codon_usage=False,
    )
    with pytest.raises(InputError, match="sequence file is empty"):
        wf.run()


def test_workflow_missing_kraken_output(tmp_path: Path, sample_fasta: Path) -> None:
    wf = Workflow(
        str(sample_fasta),
        1,
        str(tmp_path / "out.txt"),
        kraken_db=str(tmp_path / "db"),
        threads=1,
        kraken_raw_output=str(tmp_path / "missing_kraken.out"),
        run_kraken=False,
        run_codon_usage=False,
    )
    with pytest.raises(InputError, match="classifications file not found"):
        wf.run()


@patch("plasmidScreen.src.plasmidScreen.subprocess.run")
def test_run_kraken_missing_binary(
    mock_run: MagicMock,
    sample_fasta: Path,
    tmp_path: Path,
) -> None:
    mock_run.side_effect = FileNotFoundError("kraken2")
    (tmp_path / "kraken_db").mkdir()
    wf = Workflow(
        str(sample_fasta),
        1,
        str(tmp_path / "out.txt"),
        kraken_db=str(tmp_path / "kraken_db"),
        threads=1,
        kraken_raw_output=None,
        run_kraken=True,
        run_codon_usage=False,
    )
    with pytest.raises(KrakenError, match="not found"):
        wf.run()


@patch("plasmidScreen.cli.run_screen")
def test_cli_quiet_passes_quiet_mode(mock_run: MagicMock) -> None:
    mock_run.return_value = MagicMock(overall_synthetic_count=0, per_read=[])
    result = runner.invoke(
        app,
        [
            "screen",
            "reads.fa",
            "report.txt",
            "kraken_db",
            "--skip-codon-usage",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.output
    assert mock_run.call_args.kwargs["quiet_mode"] is True


@patch("plasmidScreen.cli.run_screen")
def test_cli_default_is_verbose(mock_run: MagicMock) -> None:
    mock_run.return_value = MagicMock(overall_synthetic_count=0, per_read=[])
    result = runner.invoke(
        app,
        ["screen", "reads.fa", "report.txt", "kraken_db", "--skip-codon-usage"],
    )
    assert result.exit_code == 0, result.output
    assert mock_run.call_args.kwargs["quiet_mode"] is False


@patch("plasmidScreen.cli.run_screen")
def test_cli_maps_input_error_to_exit_1(mock_run: MagicMock) -> None:
    mock_run.side_effect = InputError("Input sequence file not found: reads.fa")
    result = runner.invoke(
        app,
        ["screen", "reads.fa", "report.txt", "kraken_db", "--skip-codon-usage"],
    )
    assert result.exit_code == 1
    assert "not found" in result.output.lower() or "not found" in str(result.exception).lower()
