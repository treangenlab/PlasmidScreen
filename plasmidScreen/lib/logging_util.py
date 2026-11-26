"""Logging setup shared by the CLI (and optional library callers)."""
from __future__ import annotations

import logging

from rich.logging import RichHandler

_CONFIGURED = False


def configure_logging(
        *,
        quiet: bool = False,
        verbose: bool = False,
        force: bool = False,
) -> None:
    """
    Configure PlasmidScreen logging.

    Default: INFO milestones only (screen/map/viz complete, stage starts).
    ``verbose=True``: DEBUG (commands, timings, per-file writes).
    ``quiet=True``: ERROR and above (warnings hidden; progress bars suppressed
    separately via ``quiet_mode``).

    If both ``quiet`` and ``verbose`` are True, quiet wins.
    """
    global _CONFIGURED
    if _CONFIGURED and not force:
        return

    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    root = logging.getLogger()
    # Replace handlers so re-invoking with --quiet / --verbose / force works in tests.
    for handler in list(root.handlers):
        root.removeHandler(handler)
    handler = RichHandler(rich_tracebacks=False, show_path=False)
    handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
    root.addHandler(handler)
    root.setLevel(level)
    # Keep noisy third-party loggers down (even under --verbose).
    for name in ("numba", "matplotlib", "PIL", "fontTools"):
        logging.getLogger(name).setLevel(logging.WARNING)
    _CONFIGURED = True
