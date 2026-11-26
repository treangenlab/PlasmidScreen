"""PlasmidScreen library exceptions."""
from __future__ import annotations


class PlasmidScreenError(Exception):
    """Base error for library usage."""


class InputError(PlasmidScreenError):
    """Raised when required inputs are missing or invalid."""


class MissingCodonReferenceError(PlasmidScreenError):
    """Raised when no host taxid can be resolved to a codon usage reference."""

    def __init__(self, missing_host_taxids: list[str], data_dir: str) -> None:
        self.missing_host_taxids = missing_host_taxids
        self.data_dir = data_dir
        preview = missing_host_taxids[:20]
        more = len(missing_host_taxids) - len(preview)
        taxid_part = str(preview) if more <= 0 else f"{preview} (+{more} more)"
        super().__init__(
            f"No codon usage reference for host taxid(s) {taxid_part} in {data_dir}. "
            "CSDB covers a subset of NCBI taxa; ensure codon_tables.json is a full "
            "build and taxonomy_parents.json is present for lineage fallback. "
            "Rebuild offline with build_codon_reference() if the store is missing."
        )


class CodonReferenceNotFoundError(PlasmidScreenError):
    """Raised when the codon usage data directory or codon_tables.json is missing."""

    def __init__(self, data_dir: str) -> None:
        self.data_dir = data_dir
        super().__init__(
            f"Codon usage reference not found at {data_dir}. "
            "Run build_codon_reference() on a networked machine first."
        )


class KrakenError(PlasmidScreenError):
    """Raised when the Kraken2 subprocess fails or is unavailable."""

    def __init__(self, returncode: int, stderr: str | None = None) -> None:
        self.returncode = returncode
        self.stderr = stderr or ""
        detail = f" (stderr: {self.stderr.strip()})" if self.stderr.strip() else ""
        super().__init__(f"Kraken2 failed with exit code {returncode}{detail}")


class DiamondError(PlasmidScreenError):
    """Raised when the DIAMOND subprocess fails or is unavailable."""

    def __init__(self, returncode: int, stderr: str | None = None) -> None:
        self.returncode = returncode
        self.stderr = stderr or ""
        detail = f" (stderr: {self.stderr.strip()})" if self.stderr.strip() else ""
        super().__init__(f"DIAMOND failed with exit code {returncode}{detail}")


# Perseus disabled — re-enable with lib/perseus.py
# class PerseusError(PlasmidScreenError):
#     """Raised when the Perseus subprocess fails or is unavailable."""
#
#     def __init__(self, returncode: int, stderr: str | None = None) -> None:
#         self.returncode = returncode
#         self.stderr = stderr or ""
#         detail = f" (stderr: {self.stderr.strip()})" if self.stderr.strip() else ""
#         super().__init__(f"Perseus failed with exit code {returncode}{detail}")
