"""Shared type aliases and TypedDicts for PlasmidScreen."""
from __future__ import annotations

from typing import Literal, TypedDict

GeneSet = Literal["nuclear", "ribosomal", "mitochondrial", "plastid"]

KrakenReadInfo = tuple[str, str, int, str]
"""Kraken2 row fields: (status, taxid, read_length, kmer_info)."""

PendingCodonRead = tuple[str, str, "CdsOrf", str, str]
"""Pending codon row: (read_id, overall_taxid, cds, cds_taxid, host_taxid)."""


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
