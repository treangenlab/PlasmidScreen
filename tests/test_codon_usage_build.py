"""Tests for codon reference build helpers."""
from plasmidScreen.lib.codon_usage_build import default_reference_taxids


def test_default_reference_taxids_non_empty() -> None:
    taxids = default_reference_taxids()
    assert len(taxids) >= 50
    assert len(taxids) == len(set(taxids))
    assert "9606" in taxids
    assert "511145" in taxids
