from __future__ import annotations

import tempfile

import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binom, betabinom

"""
Per-read p-values for "more kmers mapped as synthetic than expected by chance",
with Benjamini-Hochberg (BH) FDR correction across reads.
"""


def synthetic_pvalue(k_syn, n_total, bg_rate, rho=None):
    """
    P(>= k_syn synthetic kmers out of n_total mapped kmers) under the background.

    k_syn, n_total : int or array-like, one entry per read
    bg_rate        : background fraction of mapped kmers called synthetic (0 < bg_rate < 1)
    rho            : overdispersion (0 < rho < 1) for a beta-binomial null, which
                     allows for clumped calls from overlapping kmers.
                     None -> plain binomial (assumes independent kmers).

    Returns a float for scalar input, otherwise a numpy array.
    Reads with n_total == 0 get p = 1.
    """
    k = np.asarray(k_syn, dtype=np.int64)
    n = np.asarray(n_total, dtype=np.int64)
    if np.any(k < 0) or np.any(k > n):
        raise ValueError("need 0 <= k_syn <= n_total for every read")
    if not 0.0 < bg_rate < 1.0:
        raise ValueError("bg_rate must be strictly between 0 and 1")

    if rho is None:
        p = binom.sf(k - 1, n, bg_rate)  # sf(k-1) = P(X >= k)
    else:
        if not 0.0 < rho < 1.0:
            raise ValueError("rho must be strictly between 0 and 1")
        a = bg_rate * (1.0 - rho) / rho
        b = (1.0 - bg_rate) * (1.0 - rho) / rho
        with np.errstate(invalid="ignore"):
            p = betabinom.sf(k - 1, n, a, b)

    p = np.clip(np.where(n == 0, 1.0, p), 0.0, 1.0)
    return p.item() if p.ndim == 0 else p


def bh_adjust(pvals, n_tests=None):
    """
    Benjamini-Hochberg adjusted p-values (q-values).

    pvals   : array-like of p-values
    n_tests : total number of tests (m). Defaults to len(pvals).
              Pass a larger number when pvals covers only a subset of the reads
              you tested; the missing reads are treated as having p = 1.
              This is never anti-conservative, and calls at level alpha are exact
              as long as every omitted read had p > alpha.
    """
    p = np.asarray(pvals, dtype=float)
    shape = p.shape
    p = p.ravel()
    m = p.size if n_tests is None else int(n_tests)
    if m < p.size:
        raise ValueError(f"n_tests ({m}) is smaller than the number of p-values ({p.size})")
    if p.size == 0:
        return p.reshape(shape)

    order = np.argsort(p, kind="mergesort")
    scaled = p[order] * m / np.arange(1, p.size + 1)
    q_sorted = np.minimum.accumulate(scaled[::-1])[::-1]  # enforce monotonicity
    q = np.empty_like(p)
    q[order] = np.minimum(q_sorted, 1.0)
    return q.reshape(shape)


def call_synthetic_reads(k_syn, n_total, bg_rate, n_reads_total=None, alpha=0.05, rho=None):
    """
    Score reads and call the ones with significantly more synthetic kmers than background.

    k_syn, n_total : array-like, synthetic and total mapped kmer counts per read
    bg_rate        : background rate of kmers mapped as synthetic
    n_reads_total  : total number of reads tested in the run (m for BH).
                     Defaults to len(k_syn). Set this when you only pass in a
                     subset of reads, e.g. you dropped reads with 0 synthetic
                     kmers, or kept only reads with p <= alpha while streaming.
                     Count reads with at least one mapped kmer; unmapped reads
                     weren't tested.
    alpha          : target false discovery rate
    rho            : optional beta-binomial overdispersion (see synthetic_pvalue)

    Returns (=, qvals, is_synthetic) as numpy arrays.
    """
    pvals = np.atleast_1d(synthetic_pvalue(k_syn, n_total, bg_rate, rho=rho))
    qvals = bh_adjust(pvals, n_tests=n_reads_total)
    return pvals, qvals, qvals <= alpha


def get_default_db_path(app_name: str) -> str:
    """
    Determines the platform-specific default data directory and returns
    the full path to the database file. Ensures the directory exists.

    Args:
        app_name: The name of your application (used for the folder name).

    Returns:
        The full path to the data directory.
    """
    home = Path.home()

    if sys.platform == "win32":
        base_dir = Path(os.environ.get("LOCALAPPDATA", home / "AppData" / "Local"))
        data_dir = base_dir / app_name

    elif sys.platform == "darwin":
        data_dir = home / "Library" / "Application Support" / app_name

    else:
        xdg_data = os.environ.get("XDG_DATA_HOME")
        if xdg_data:
            data_dir = Path(xdg_data) / app_name
        else:
            data_dir = home / ".local" / "share" / app_name

    try:
        data_dir.mkdir(parents=True, exist_ok=True)

        # Test write permissions explicitly (handles cases where mkdir succeeds but write fails)
        if not os.access(data_dir, os.W_OK):
            raise PermissionError

    except (PermissionError, OSError):
        # 4. Fallback Strategy: Use the OS temporary directory
        # tempfile.gettempdir() automatically resolves to /tmp on Linux/macOS
        data_dir = Path(tempfile.gettempdir()) / app_name
        data_dir.mkdir(parents=True, exist_ok=True)

        # Alert the user/logs that data persistence is volatile
        print(
            f"Warning: Primary data directory unwritable. "
            f"Falling back to temporary storage: {data_dir}",
            file=sys.stderr
        )

    return str(data_dir)
