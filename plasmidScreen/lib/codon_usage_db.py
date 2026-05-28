"""
Read-only codon usage reference store for airgapped CAI scoring.

Populate with build_codon_reference() (see codon_usage_build.py) before screening.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Optional

from Bio.Data import CodonTable
from plasmidScreen.lib.exceptions import CodonReferenceNotFoundError, MissingCodonReferenceError
from plasmidScreen.lib.funcs import get_default_db_path
from plasmidScreen.lib.types import CodonTableEntry

CODON_TABLES_FILE = "codon_tables.json"
TAXONOMY_PARENTS_FILE = "taxonomy_parents.json"

_CODON_TO_AA = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
for _stop in ("TAA", "TAG", "TGA"):
    _CODON_TO_AA[_stop] = "*"


def default_codon_usage_dir(app_name: str = "PlasmidScreen") -> Path:
    return Path(get_default_db_path(app_name)) / "codon_usage"


def default_codon_usage_db_path(app_name: str = "PlasmidScreen") -> Path:
    """Backward-compatible alias: returns the codon usage data directory."""
    return default_codon_usage_dir(app_name)


def cai_weights_from_frequencies(frequencies: dict[str, float]) -> dict[str, float]:
    """Sharp & Li relative adaptiveness weights per codon (0–1 within each amino acid)."""
    by_aa: dict[str, dict[str, float]] = {}
    for codon, freq in frequencies.items():
        aa = _CODON_TO_AA.get(codon)
        if aa is None or aa == "*":
            continue
        by_aa.setdefault(aa, {})[codon] = max(freq, 0.0)

    weights: dict[str, float] = {}
    for _aa, codon_freqs in by_aa.items():
        if not codon_freqs:
            continue
        max_f = max(codon_freqs.values()) or 1.0
        for codon, freq in codon_freqs.items():
            weights[codon] = freq / max_f if max_f > 0 else 0.0
    return weights


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open() as f:
        return json.load(f)


def _save_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def parse_taxonomy_nodes(nodes_dmp: Path) -> dict[str, str]:
    """Parse NCBI nodes.dmp into taxid -> parent_taxid."""
    parents: dict[str, str] = {}
    with nodes_dmp.open() as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 2:
                continue
            taxid, parent = parts[0], parts[1]
            if taxid == parent:
                continue
            parents[taxid] = parent
    return parents


def taxids_from_kraken_output(kraken_path: str | Path) -> set[str]:
    """Collect unique classified taxids from Kraken2 output."""
    taxids: set[str] = set()
    with open(kraken_path, "r", buffering=1024 * 1024) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            taxid = parts[2]
            if taxid not in ("0", "U", "unclassified"):
                taxids.add(taxid)
    return taxids


class CodonUsageStore:
    """Read-only JSON codon usage store for runtime (airgapped).

    Use CodonUsageStore.writable() during offline builds only.
    """

    data_dir: Path
    tables_path: Path
    taxonomy_path: Path
    _writable: bool
    _tables: dict[str, CodonTableEntry]
    _parents: dict[str, str]
    _dirty: bool

    def __init__(
        self,
        data_dir: str | Path,
        *,
        create: bool = False,
        _tables: dict[str, CodonTableEntry] | None = None,
        _parents: dict[str, str] | None = None,
    ) -> None:
        self.data_dir = Path(data_dir)
        self.tables_path = self.data_dir / CODON_TABLES_FILE
        self.taxonomy_path = self.data_dir / TAXONOMY_PARENTS_FILE
        self._writable = create

        if create:
            self.data_dir.mkdir(parents=True, exist_ok=True)
            self._tables = _tables if _tables is not None else _load_json(self.tables_path)
            self._parents = _parents if _parents is not None else _load_json(self.taxonomy_path)
            self._dirty = False
        else:
            if not self.tables_path.exists():
                raise CodonReferenceNotFoundError(str(self.data_dir))
            self._tables = _load_json(self.tables_path)
            self._parents = _load_json(self.taxonomy_path)

    @classmethod
    def load(cls, data_dir: str | Path) -> "CodonUsageStore":
        """Load a read-only store (raises if not built)."""
        return cls(data_dir, create=False)

    @classmethod
    def writable(cls, data_dir: str | Path) -> "CodonUsageStore":
        """Open for offline building (call save() when done)."""
        return cls(data_dir, create=True)

    def save(self) -> None:
        if not self._writable:
            raise RuntimeError("Cannot save a read-only CodonUsageStore")
        if getattr(self, "_dirty", False):
            _save_json(self.tables_path, self._tables)
            _save_json(self.taxonomy_path, self._parents)
            self._dirty = False

    def has_taxonomy(self) -> bool:
        return bool(self._parents)

    def taxonomy_parents(self) -> dict[str, str]:
        """Copy of taxid -> parent_taxid lineage map (empty if taxonomy not loaded)."""
        return dict(self._parents)

    def has_codon_table(self, taxid: str | int) -> bool:
        return str(taxid) in self._tables

    def get_frequencies(self, taxid: str | int) -> Optional[dict[str, float]]:
        entry = self._tables.get(str(taxid))
        if not entry:
            return None
        return entry.get("frequencies")

    def get_cai_weights(self, taxid: str | int) -> Optional[dict[str, float]]:
        freqs = self.get_frequencies(taxid)
        if freqs is None:
            return None
        return cai_weights_from_frequencies(freqs)

    def parent_taxid(self, taxid: str | int) -> Optional[str]:
        return self._parents.get(str(taxid))

    def resolve_reference_taxid(self, taxid: str | int) -> Optional[str]:
        if str(taxid) in ("0", ""):
            return None
        current = str(taxid)
        seen: set[str] = set()
        while current and current not in seen:
            if self.has_codon_table(current):
                return current
            seen.add(current)
            parent = self.parent_taxid(current)
            if parent is None or parent == current:
                break
            current = parent
        return None

    def get_cai_weights_for_host(
        self, host_taxid: str | int
    ) -> tuple[Optional[dict[str, float]], Optional[str]]:
        ref = self.resolve_reference_taxid(host_taxid)
        if ref is None:
            return None, None
        return self.get_cai_weights(ref), ref

    def missing_host_taxids(self, host_taxids: Iterable[str | int]) -> list[str]:
        """Host taxids with no resolvable reference table in this store."""
        missing = []
        for taxid in host_taxids:
            if str(taxid) in ("0", ""):
                continue
            if self.resolve_reference_taxid(taxid) is None:
                missing.append(str(taxid))
        return sorted(set(missing))

    def require_host_taxids(self, host_taxids: Iterable[str | int]) -> None:
        """Raise MissingCodonReferenceError if any host lacks a reference table."""
        missing = self.missing_host_taxids(host_taxids)
        if missing:
            raise MissingCodonReferenceError(missing, str(self.data_dir))

    def set_codon_table(
        self,
        taxid: str | int,
        frequencies: dict[str, float],
        *,
        scientific_name: str | None = None,
        source: str = "csdb",
    ) -> None:
        if not self._writable:
            raise RuntimeError("Cannot modify a read-only CodonUsageStore")
        entry: dict[str, Any] = {"source": source, "frequencies": frequencies}
        if scientific_name:
            entry["scientific_name"] = scientific_name
        self._tables[str(taxid)] = entry
        self._dirty = True

    def load_taxonomy_from_nodes(self, nodes_dmp: Path) -> int:
        if not self._writable:
            raise RuntimeError("Cannot modify a read-only CodonUsageStore")
        self._parents = parse_taxonomy_nodes(nodes_dmp)
        self._dirty = True
        return len(self._parents)


# Backward-compatible alias
CodonUsageDB = CodonUsageStore
