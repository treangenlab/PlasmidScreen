"""Offline codon reference builder (network access only during build)."""
from __future__ import annotations

import logging
import re
import tarfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable

from plasmidScreen.lib.codon_usage_db import (
    CODON_TABLES_FILE,
    CodonUsageStore,
    parse_taxonomy_nodes,
)
from plasmidScreen.lib.models import BuildCodonReferenceResult

KAZUSA_URL = (
    "http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi"
    "?aa=1&style=N&species={taxid}"
)
KAZUSA_CODON_RE = re.compile(r"([ATGCU]{3})\s+([A-Z*])\s+([\d.]+)")
NCBI_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

from Bio.Data import CodonTable

_CODON_TO_AA = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
for _stop in ("TAA", "TAG", "TGA"):
    _CODON_TO_AA[_stop] = "*"


def _rna_to_dna(codon: str) -> str:
    return codon.replace("U", "T").replace("u", "t")


def fetch_kazusa_frequencies(taxid: str | int, timeout: int = 30) -> dict[str, float]:
    """Download codon relative frequencies from Kazusa (build-time only)."""
    url = KAZUSA_URL.format(taxid=int(taxid))
    try:
        with urllib.request.urlopen(url, timeout=timeout) as handle:
            html = handle.read().decode(errors="replace")
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Failed to fetch Kazusa codon table for taxid {taxid}: {exc}") from exc

    if "not found" in html.lower() and "<title>" in html.lower():
        raise RuntimeError(f"No Kazusa codon usage table for taxonomy ID {taxid}")

    frequencies: dict[str, float] = {}
    for codon, _aa, rel_freq in KAZUSA_CODON_RE.findall(html):
        dna_codon = _rna_to_dna(codon)
        if dna_codon in _CODON_TO_AA and _CODON_TO_AA[dna_codon] != "*":
            frequencies[dna_codon] = float(rel_freq)

    if not frequencies:
        raise RuntimeError(f"Could not parse Kazusa codon table for taxid {taxid}")

    return frequencies


def download_ncbi_taxdump(dest_dir: Path) -> Path:
    """Download and extract nodes.dmp from NCBI taxdump."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    archive = dest_dir / "taxdump.tar.gz"
    if not archive.exists():
        logging.info("Downloading NCBI taxdump from %s", NCBI_TAXDUMP_URL)
        urllib.request.urlretrieve(NCBI_TAXDUMP_URL, archive)

    nodes_path = dest_dir / "nodes.dmp"
    if not nodes_path.exists():
        logging.info("Extracting nodes.dmp from taxdump archive")
        with tarfile.open(archive, "r:gz") as tar:
            try:
                tar.extract("nodes.dmp", path=dest_dir, filter="data")
            except TypeError:
                tar.extract("nodes.dmp", path=dest_dir)

    return nodes_path


def build_codon_reference(
    data_dir: str | Path,
    taxids: Iterable[str | int],
    *,
    include_taxonomy: bool = True,
    taxdump_dir: str | Path | None = None,
    fetch_timeout: int = 30,
) -> BuildCodonReferenceResult:
    """
    Build codon_tables.json (and optional taxonomy_parents.json) for airgapped use.

    Must be run on a machine with network access before screening.
    """
    data_dir = Path(data_dir)
    taxid_list = sorted({str(t) for t in taxids if str(t) not in ("0", "")})

    added: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []

    store = CodonUsageStore.writable(data_dir)

    if include_taxonomy and not store.has_taxonomy():
        tdir = Path(taxdump_dir) if taxdump_dir else data_dir.parent / "taxdump"
        nodes = download_ncbi_taxdump(tdir)
        count = store.load_taxonomy_from_nodes(nodes)
        store.save()
        logging.info("Loaded %d taxonomy parent links", count)

    for taxid in taxid_list:
        if store.has_codon_table(taxid):
            skipped.append(taxid)
            continue
        logging.info("Fetching Kazusa codon usage for taxid %s", taxid)
        try:
            frequencies = fetch_kazusa_frequencies(taxid, timeout=fetch_timeout)
            store.set_codon_table(taxid, frequencies, source="kazusa")
            added.append(taxid)
        except RuntimeError as exc:
            logging.warning("Could not add taxid %s: %s", taxid, exc)
            failed.append(taxid)

    store.save()

    if not (data_dir / CODON_TABLES_FILE).exists() and not added and not skipped:
        raise RuntimeError(f"No codon tables written to {data_dir}")

    return BuildCodonReferenceResult(
        data_dir=data_dir,
        taxids_requested=taxid_list,
        taxids_added=added,
        taxids_skipped=skipped,
        taxids_failed=failed,
    )
