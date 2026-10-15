"""Shared type aliases and TypedDicts for PlasmidScreen."""
from __future__ import annotations

from enum import Enum
from typing import Literal, TypedDict


class GeneSet(str, Enum):
    """CSDB gene-set categories used when building codon tables."""

    NUCLEAR = "nuclear"
    RIBOSOMAL = "ribosomal"
    MITOCHONDRIAL = "mitochondrial"
    PLASTID = "plastid"


KrakenReadInfo = tuple[str, str, int, str]
"""Kraken2 row fields: (status, taxid, read_length, kmer_info)."""


class CdsOrf(TypedDict):
    """Longest open reading frame detected in a read."""

    strand: Literal["+", "-"]
    start: int
    end: int
    length: int
    seq: str


class CodonTableEntry(TypedDict, total=False):
    """Single species entry in codon_tables.json."""

    source: str
    scientific_name: str
    frequencies: dict[str, float]
