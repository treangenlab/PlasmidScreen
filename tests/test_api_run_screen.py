"""Tests for the public run_screen() API façade."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from plasmidScreen.api import run_screen
from plasmidScreen.lib.models import BACKGROUND_RATE


def test_run_screen_requires_diamond_db() -> None:
    with pytest.raises(ValueError, match="diamond_db"):
        run_screen(
            "reads.fa",
            1,
            "kraken_db",
            run_codon_usage=True,
            run_diamond=True,
            diamond_db=None,
        )


def test_run_screen_requires_diamond_output_when_not_running() -> None:
    with pytest.raises(ValueError, match="diamond_output_path"):
        run_screen(
            "reads.fa",
            1,
            "kraken_db",
            run_codon_usage=True,
            run_diamond=False,
            diamond_output_path=None,
        )


def test_run_screen_requires_diamond_output_for_debug_write() -> None:
    with pytest.raises(ValueError, match="diamond_output_path"):
        run_screen(
            "reads.fa",
            1,
            "kraken_db",
            run_codon_usage=True,
            diamond_db="db.dmnd",
            debug_write_diamond_output=True,
            diamond_output_path=None,
        )


def test_run_screen_rejects_bad_memory_mode() -> None:
    with pytest.raises(ValueError, match="memory_mode"):
        run_screen(
            "reads.fa",
            1,
            "kraken_db",
            run_codon_usage=False,
            memory_mode=3,
        )


@patch("plasmidScreen.api.Workflow")
def test_run_screen_delegates_defaults_to_workflow(mock_workflow_cls: MagicMock) -> None:
    mock_instance = MagicMock()
    mock_workflow_cls.return_value = mock_instance
    mock_instance.run.return_value = MagicMock(name="ScreenResult")

    result = run_screen(
        Path("reads.fa"),
        1,
        Path("kraken_db"),
        run_codon_usage=False,
    )

    mock_workflow_cls.assert_called_once()
    kwargs = mock_workflow_cls.call_args.kwargs
    # API defaults (distinct from Workflow.__init__ defaults)
    assert kwargs["seq_type"] == 1
    assert kwargs["engineered_kmer_threshold"] == 10
    assert kwargs["codon_cai_engineered_threshold"] == 0.6
    assert kwargs["window_size"] == 200
    assert kwargs["background_rate"] == BACKGROUND_RATE
    assert kwargs["mem_mode"] == 0
    assert kwargs["run_codon_usage"] is False
    mock_instance.run.assert_called_once()
    assert result is mock_instance.run.return_value


@patch("plasmidScreen.api.Workflow")
def test_run_screen_passes_background_rate(mock_workflow_cls: MagicMock) -> None:
    mock_workflow_cls.return_value.run.return_value = MagicMock()
    run_screen(
        "reads.fa",
        0,
        "kraken_db",
        run_codon_usage=False,
        background_rate=0.01,
        memory_mode=2,
    )
    kwargs = mock_workflow_cls.call_args.kwargs
    assert kwargs["seq_type"] == 0
    assert kwargs["background_rate"] == 0.01
    assert kwargs["mem_mode"] == 2
