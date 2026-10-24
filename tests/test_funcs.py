"""Unit tests for statistical helpers in plasmidScreen.lib.funcs."""
from __future__ import annotations

import numpy as np
import pytest

from plasmidScreen.lib.funcs import bh_adjust, call_synthetic_reads, synthetic_pvalue


def test_synthetic_pvalue_scalar_binomial() -> None:
    # Under bg_rate=0.5, P(X >= 3 | n=3) = 0.125
    p = synthetic_pvalue(3, 3, bg_rate=0.5)
    assert p == pytest.approx(0.125)


def test_synthetic_pvalue_zero_total_is_one() -> None:
    assert synthetic_pvalue(0, 0, bg_rate=0.01) == 1.0
    p = synthetic_pvalue([0, 5], [0, 10], bg_rate=0.01)
    assert p[0] == 1.0
    assert 0.0 <= p[1] <= 1.0


def test_synthetic_pvalue_rejects_bad_inputs() -> None:
    with pytest.raises(ValueError, match="bg_rate"):
        synthetic_pvalue(1, 10, bg_rate=0.0)
    with pytest.raises(ValueError, match="bg_rate"):
        synthetic_pvalue(1, 10, bg_rate=1.0)
    with pytest.raises(ValueError, match="0 <= k_syn"):
        synthetic_pvalue(5, 3, bg_rate=0.1)
    with pytest.raises(ValueError, match="rho"):
        synthetic_pvalue(1, 10, bg_rate=0.1, rho=0.0)


def test_bh_adjust_empty_and_monotonic() -> None:
    assert bh_adjust([]).size == 0
    p = np.array([0.01, 0.04, 0.03, 0.20])
    q = bh_adjust(p)
    assert q.shape == p.shape
    assert np.all((q >= 0.0) & (q <= 1.0))
    # Sorted q-values must be non-decreasing
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= -1e-12)


def test_bh_adjust_n_tests_inflates() -> None:
    p = np.array([0.01])
    q_default = bh_adjust(p)
    q_larger = bh_adjust(p, n_tests=100)
    assert q_larger[0] >= q_default[0]
    with pytest.raises(ValueError, match="n_tests"):
        bh_adjust([0.1, 0.2], n_tests=1)


def test_call_synthetic_reads_flags_strong_signal() -> None:
    # One read with all kmers synthetic vs tiny background → significant
    k_syn = np.array([0, 50], dtype=np.int64)
    n_total = np.array([100, 50], dtype=np.int64)
    pvals, qvals, is_syn = call_synthetic_reads(
        k_syn, n_total, bg_rate=0.0003, alpha=0.05
    )
    assert pvals.shape == (2,)
    assert qvals.shape == (2,)
    assert is_syn.dtype == bool
    assert bool(is_syn[0]) is False
    assert bool(is_syn[1]) is True
    assert np.all(is_syn == (qvals <= 0.05))
