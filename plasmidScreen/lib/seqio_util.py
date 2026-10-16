"""Helpers for reading plain or gzip-compressed FASTA/FASTQ."""
from __future__ import annotations

import gzip
import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Literal, TextIO

from plasmidScreen.lib.exceptions import InputError

SequenceFormat = Literal["fasta", "fastq"]
logger = logging.getLogger(__name__)


def is_gzip_path(path: str | Path) -> bool:
    """True if the path looks gzip-compressed (``.gz`` suffix)."""
    return Path(path).name.lower().endswith(".gz")


def sequence_format(path: str | Path) -> SequenceFormat:
    """
    Infer ``fasta`` vs ``fastq`` from the filename, ignoring a trailing ``.gz``.

    Examples: ``reads.fa``, ``reads.fasta.gz``, ``reads.fq.gz`` → fasta/fastq.
    Unknown suffixes default to ``fasta``.
    """
    suffixes = [s.lower() for s in Path(path).suffixes]
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    if not suffixes:
        return "fasta"
    last = suffixes[-1]
    if last in (".fastq", ".fq"):
        return "fastq"
    return "fasta"


def is_fastq_path(path: str | Path) -> bool:
    return sequence_format(path) == "fastq"


@contextmanager
def open_text_maybe_gzip(path: str | Path) -> Iterator[TextIO]:
    """Open a text file, transparently decompressing ``.gz`` inputs."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"Sequence file not found: {path}")
    try:
        if is_gzip_path(path):
            # Probe the header so corrupt .gz fails before callers start parsing.
            with gzip.open(path, "rb") as probe:
                probe.read(1)
            handle: TextIO = gzip.open(path, "rt", encoding="utf-8")  # type: ignore[assignment]
        else:
            handle = open(path, "r", encoding="utf-8")
    except gzip.BadGzipFile as exc:
        logger.error("Invalid gzip content in %s", path)
        raise InputError(f"Invalid gzip file: {path}") from exc
    except OSError as exc:
        logger.error("Failed to open sequence file %s: %s", path, exc)
        raise InputError(f"Could not open sequence file {path}: {exc}") from exc
    try:
        yield handle
    finally:
        handle.close()
