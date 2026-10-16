"""Logging setup shared by the CLI (and optional library callers)."""
from __future__ import annotations

import logging

from rich.logging import RichHandler

_CONFIGURED = False


def configure_logging(*, quiet: bool = False, force: bool = False) -> None:
    """
    Configure PlasmidScreen logging.

    When ``quiet`` is False (default for CLI), INFO and above are shown.
    When ``quiet`` is True, only WARNING and above are shown (progress bars
    are also suppressed separately via ``quiet_mode``).
    """
    global _CONFIGURED
    if _CONFIGURED and not force:
        return

    level = logging.WARNING if quiet else logging.INFO
    root = logging.getLogger()
    # Replace handlers so re-invoking with --quiet / force works in tests.
    for handler in list(root.handlers):
        root.removeHandler(handler)
    handler = RichHandler(rich_tracebacks=False, show_path=False)
    handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
    root.addHandler(handler)
    root.setLevel(level)
    # Keep noisy third-party loggers down unless debugging.
    logging.getLogger("numba").setLevel(logging.WARNING)
    _CONFIGURED = True
