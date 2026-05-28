"""PlasmidScreen library exceptions."""
from __future__ import annotations

class PlasmidScreenError(Exception):
    """Base error for library usage."""


class MissingCodonReferenceError(PlasmidScreenError):
    """Raised when host taxids lack a pre-built codon usage reference (airgapped runs)."""

    def __init__(self, missing_host_taxids: list[str], data_dir: str) -> None:
        self.missing_host_taxids = missing_host_taxids
        self.data_dir = data_dir
        super().__init__(
            f"No codon usage reference for host taxid(s) {missing_host_taxids} in {data_dir}. "
            "Build the reference offline with build_codon_reference() before screening."
        )


class CodonReferenceNotFoundError(PlasmidScreenError):
    """Raised when the codon usage data directory or codon_tables.json is missing."""

    def __init__(self, data_dir: str) -> None:
        self.data_dir = data_dir
        super().__init__(
            f"Codon usage reference not found at {data_dir}. "
            "Run build_codon_reference() on a networked machine first."
        )
