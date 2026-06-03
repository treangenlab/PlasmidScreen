"""Shared type aliases and TypedDicts for PlasmidScreen."""
from __future__ import annotations

from typing import Literal, TypedDict

GeneSet = Literal["nuclear", "ribosomal", "mitochondrial", "plastid"]

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
