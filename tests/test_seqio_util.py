"""Tests for gzip-aware FASTA/FASTQ helpers."""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from plasmidScreen.lib.seqio_util import (
    is_gzip_path,
    is_fastq_path,
    open_text_maybe_gzip,
    sequence_format,
)
from plasmidScreen.src.analyze_codon_usage import write_subset_reads


@pytest.mark.parametrize(
    "name,fmt,gzipped",
    [
        ("reads.fa", "fasta", False),
        ("reads.fasta", "fasta", False),
        ("reads.fq", "fastq", False),
        ("reads.fastq", "fastq", False),
        ("reads.fa.gz", "fasta", True),
        ("reads.fasta.gz", "fasta", True),
        ("reads.fq.gz", "fastq", True),
        ("reads.fastq.gz", "fastq", True),
    ],
)
def test_sequence_format_and_gzip_detection(name: str, fmt: str, gzipped: bool) -> None:
    path = Path(name)
    assert sequence_format(path) == fmt
    assert is_fastq_path(path) == (fmt == "fastq")
    assert is_gzip_path(path) is gzipped


def test_open_text_maybe_gzip_roundtrip(tmp_path: Path) -> None:
    plain = tmp_path / "x.fa"
    plain.write_text(">a\nACGT\n", encoding="utf-8")
    with open_text_maybe_gzip(plain) as h:
        assert h.read() == ">a\nACGT\n"

    gz = tmp_path / "x.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as out:
        out.write(">b\nTTAA\n")
    with open_text_maybe_gzip(gz) as h:
        assert h.read() == ">b\nTTAA\n"


def test_write_subset_reads_from_gzip(sample_fasta: Path, tmp_path: Path) -> None:
    gz = tmp_path / "reads.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as out:
        out.write(sample_fasta.read_text())

    dest = tmp_path / "natural.fa"
    n = write_subset_reads(gz, {"read_natural_1"}, dest)
    assert n == 1
    text = dest.read_text()
    assert ">read_natural_1" in text
    assert "read_synthetic_1" not in text
