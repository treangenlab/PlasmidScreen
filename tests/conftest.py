"""Shared pytest fixtures."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from plasmidScreen.src.codon_usage.codon_usage_db import CodonUsageStore
from Bio.Data import CodonTable

FIXTURES = Path(__file__).parent / "fixtures"

_CODON_TO_AA = CodonTable.unambiguous_dna_by_id[1].forward_table.copy()
for _stop in ("TAA", "TAG", "TGA"):
    _CODON_TO_AA[_stop] = "*"


def _uniform_frequencies() -> dict[str, float]:
    """Equal relative frequency per sense codon (within each AA family)."""
    by_aa: dict[str, list[str]] = {}
    for codon, aa in _CODON_TO_AA.items():
        if aa == "*":
            continue
        by_aa.setdefault(aa, []).append(codon)
    freqs: dict[str, float] = {}
    for codons in by_aa.values():
        w = 1.0 / len(codons)
        for c in codons:
            freqs[c] = w
    return freqs


@pytest.fixture
def uniform_frequencies() -> dict[str, float]:
    return _uniform_frequencies()


@pytest.fixture
def codon_store_dir(tmp_path: Path, uniform_frequencies: dict[str, float]) -> Path:
    """Minimal airgapped codon reference with taxids 562 and 561."""
    data_dir = tmp_path / "codon_usage"
    data_dir.mkdir()
    tables = {
        "562": {"source": "test", "frequencies": uniform_frequencies},
        "561": {"source": "test", "frequencies": uniform_frequencies},
    }
    (data_dir / "codon_tables.json").write_text(json.dumps(tables, indent=2))
    parents = {"562": "561", "561": "561"}
    (data_dir / "taxonomy_parents.json").write_text(json.dumps(parents, indent=2))
    return data_dir


@pytest.fixture
def codon_store(codon_store_dir: Path) -> CodonUsageStore:
    return CodonUsageStore.load(codon_store_dir)


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    """Two reads with a short ORF (ATG start, no in-frame stop)."""
    fasta = tmp_path / "reads.fa"
    # 100 bp so k=35 yields 66 k-mer tokens
    seq = "ATG" + "AAA" * 32 + "TAA"  # stop at end of frame 0 segment; longest ORF still usable
    seq = ("ATG" + "AAA" * 32)[:100]
    if len(seq) < 100:
        seq = (seq + "A" * 100)[:100]
    fasta.write_text(
        ">read_natural_1\n"
        f"{seq}\n"
        ">read_synthetic_1\n"
        f"{seq}\n"
    )
    return fasta


@pytest.fixture
def sample_kraken_out(tmp_path: Path) -> Path:
    """Kraken lines: natural (562 k-mers) vs synthetic (32630 k-mers)."""
    length = 100
    k = 35
    n_kmers = length - k + 1
    natural_kmer_field = f"562:{max(n_kmers, 1)}"
    synthetic_kmer_field = f"32630:{max(n_kmers, 30)}"
    lines = [
        f"C\tread_natural_1\t562\t{length}\t0:1\t{natural_kmer_field}\n",
        f"C\tread_synthetic_1\t562\t{length}\t0:1\t{synthetic_kmer_field}\n",
        f"U\tread_unclassified\t0\t{length}\t0:1\t0:{max(n_kmers, 1)}\n",
    ]
    path = tmp_path / "kraken.out"
    path.write_text("".join(lines))
    return path


@pytest.fixture
def fixtures_kraken_out() -> Path:
    return FIXTURES / "kraken.out"


@pytest.fixture
def fixtures_fasta() -> Path:
    return FIXTURES / "reads.fa"
