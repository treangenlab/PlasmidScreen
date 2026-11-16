"""Tests for per-plasmid pileup aggregation and figures."""
from __future__ import annotations

from pathlib import Path

from plasmidScreen.lib.models import (
    CodonAdaptationResult,
    EngineeredScanResult,
    QueryMatchConfidence,
    ReadEngineeringLabel,
    ReadFlagDetail,
    SampleReferenceRank,
    ScreenResult,
)
from plasmidScreen.api import visualize_engineered_mappings
from plasmidScreen.src.visualization_code import (
    build_engineered_mapping_hits,
    create_pile_ups,
    group_hits_by_plasmid,
    load_genbank_annotations,
    render_plasmid_pileup,
)


def _hit(
    accession: str,
    *,
    query_id: str,
    target_start: int,
    target_end: int,
    target_length: int = 1000,
    match_score: float = 100.0,
    strand: str = "+",
) -> QueryMatchConfidence:
    return QueryMatchConfidence(
        accession_id=accession,
        database_source="test_db",
        match_score=match_score,
        query_id=query_id,
        target_start=target_start,
        target_end=target_end,
        target_length=target_length,
        functional_metadata={"strand": strand},
    )


def _screen_with_hits(hits_by_read: dict[str, tuple[QueryMatchConfidence, ...]]) -> ScreenResult:
    labels = []
    per_read = []
    for rid, hits in hits_by_read.items():
        labels.append(
            ReadEngineeringLabel(read_id=rid, label="Synthetic", eng_kmer_coverage=1.0)
        )
        summary = hits[0] if hits else None
        per_read.append(
            ReadFlagDetail(
                read_id=rid,
                kmer_label="Synthetic",
                engineered_by_kmer_scan=True,
                engineered_overall=True,
                overall_label="Synthetic",
                database_match_summary=summary,
                database_match_hits=hits,
            )
        )
    return ScreenResult(
        engineered_scan=EngineeredScanResult(
            labels=labels,
            synthetic_count=len(labels),
            natural_count=0,
        ),
        codon_adaptation=CodonAdaptationResult(labels=[]),
        engineered_kmer_threshold=10,
        engineered_kmer_window_size=200,
        codon_cai_engineered_threshold=0.6,
        per_read=per_read,
    )


def test_group_hits_by_plasmid_aggregates_queries() -> None:
    result = _screen_with_hits(
        {
            "read_a": (
                _hit("plasmid_X", query_id="read_a", target_start=10, target_end=200),
                _hit("plasmid_Y", query_id="read_a", target_start=0, target_end=50),
            ),
            "read_b": (
                _hit("plasmid_X", query_id="read_b", target_start=150, target_end=400),
            ),
        }
    )
    grouped = group_hits_by_plasmid(result.per_read_engineered)
    assert set(grouped) == {"plasmid_X", "plasmid_Y"}
    assert {h.query_id for h in grouped["plasmid_X"]} == {"read_a", "read_b"}
    assert len(grouped["plasmid_Y"]) == 1


def test_build_engineered_mapping_hits_mirrors_pileup() -> None:
    result = _screen_with_hits(
        {
            "read_a": (
                _hit("plasmid_X", query_id="read_a", target_start=10, target_end=200),
            ),
            "read_b": (
                _hit("plasmid_X", query_id="read_b", target_start=150, target_end=400),
            ),
        }
    )
    summaries = build_engineered_mapping_hits(result.per_read_engineered)
    assert len(summaries) == 1
    summary = summaries[0]
    assert summary.accession_id == "plasmid_X"
    assert summary.n_queries == 2
    assert summary.n_alignments == 2
    assert summary.query_ids == ("read_a", "read_b")
    assert summary.target_length == 1000
    assert summary.covered_bp == 390  # 10–400
    assert summary.max_coverage == 2  # overlap 150–200
    assert summary.alignments[0].target_start == 10
    assert summary.alignments[0].target_end == 200
    assert summary.alignments[1].query_id == "read_b"


def test_create_pile_ups_attaches_engineered_mapping_hits(tmp_path: Path) -> None:
    result = _screen_with_hits(
        {
            "read_a": (_hit("plasmid_X", query_id="read_a", target_start=0, target_end=100),),
            "read_b": (_hit("plasmid_Y", query_id="read_b", target_start=20, target_end=80),),
        }
    )
    assert result.engineered_mapping_hits == ()
    create_pile_ups(result, tmp_path, fetch_genbank=False)
    assert {m.accession_id for m in result.engineered_mapping_hits} == {
        "plasmid_X",
        "plasmid_Y",
    }
    by_acc = {m.accession_id: m for m in result.engineered_mapping_hits}
    assert by_acc["plasmid_X"].alignments[0].query_id == "read_a"
    assert by_acc["plasmid_X"].alignments[0].target_start == 0
    assert by_acc["plasmid_X"].alignments[0].target_end == 100


def test_visualize_api_respects_top_refs_and_writes_sidecar(tmp_path: Path) -> None:
    result = _screen_with_hits(
        {
            "read_a": (_hit("plasmid_X", query_id="read_a", target_start=0, target_end=100),),
            "read_b": (_hit("plasmid_Y", query_id="read_b", target_start=20, target_end=80),),
        }
    )
    result.top_engineered_references = (
        SampleReferenceRank(accession_id="plasmid_X", n_reads=1, sum_bitscore=100.0),
    )
    viz = visualize_engineered_mappings(
        result, tmp_path / "viz", fetch_genbank=False, write_json=True
    )
    assert viz.report_json_path is not None and viz.report_json_path.is_file()
    assert {m.accession_id for m in viz.engineered_mapping_hits} == {"plasmid_X"}
    assert all(p.name.startswith("pileup_plasmid_X") for p in viz.pileup_paths)
    assert result.engineered_mapping_hits == viz.engineered_mapping_hits


def test_render_plasmid_pileup_builds_figure() -> None:
    fig = render_plasmid_pileup(
        "plasmid_X",
        [
            _hit("plasmid_X", query_id="read_a", target_start=0, target_end=250),
            _hit(
                "plasmid_X",
                query_id="read_b",
                target_start=700,
                target_end=900,
                strand="-",
            ),
        ],
    )
    assert fig._suptitle is not None
    assert "plasmid_X" in fig._suptitle.get_text()
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_create_pile_ups_writes_one_png_per_plasmid(tmp_path: Path) -> None:
    result = _screen_with_hits(
        {
            "read_a": (_hit("plasmid_X", query_id="read_a", target_start=0, target_end=100),),
            "read_b": (_hit("plasmid_Y", query_id="read_b", target_start=20, target_end=80),),
        }
    )
    written = create_pile_ups(result, tmp_path, fetch_genbank=False)
    names = sorted(p.name for p in written)
    assert names == ["pileup_plasmid_X.png", "pileup_plasmid_Y.png"]
    assert (tmp_path / "pileup_plasmid_X.png").stat().st_size > 0


def test_compound_location_explodes_into_parts() -> None:
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    from plasmidScreen.src.visualization_code import annotations_from_features

    record = SeqRecord(Seq("A" * 3000), id="wrap")
    # Origin-spanning gene like Mob on pBBR1.
    feature = SeqFeature(
        CompoundLocation(
            [
                FeatureLocation(2496, 2687, strand=-1),
                FeatureLocation(0, 799, strand=-1),
            ],
            "order",
        ),
        type="gene",
        qualifiers={"gene": ["Mob"], "product": ["mobilization protein"]},
    )
    anns = annotations_from_features([feature])
    assert len(anns) == 2
    spans = {(a.start, a.end) for a in anns}
    assert spans == {(0, 799), (2496, 2687)}
    assert all(a.label.startswith("Mob") for a in anns)


def test_load_genbank_annotations_and_render(tmp_path: Path) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("A" * 1000), id="plasmid_X", name="plasmid_X")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(100, 400),
            type="CDS",
            qualifiers={"gene": ["bla"]},
        ),
        SeqFeature(
            FeatureLocation(500, 700, strand=-1),
            type="promoter",
            qualifiers={"note": ["P1"]},
        ),
    ]
    gb = tmp_path / "plasmid_X.gb"
    with gb.open("w") as fh:
        SeqIO.write(record, fh, "genbank")

    ann = load_genbank_annotations(gb)
    assert "plasmid_X" in ann
    assert any(a.label == "bla" for a in ann["plasmid_X"])

    fig = render_plasmid_pileup(
        "plasmid_X",
        [_hit("plasmid_X", query_id="read_a", target_start=120, target_end=350)],
        annotations=ann["plasmid_X"],
    )
    assert fig._suptitle is not None
    title = fig._suptitle.get_text()
    assert "plasmid_X" in title
    assert "feature" in title
    import matplotlib.pyplot as plt

    plt.close(fig)

    result = _screen_with_hits(
        {"read_a": (_hit("plasmid_X", query_id="read_a", target_start=120, target_end=350),)}
    )
    written = create_pile_ups(
        result, tmp_path / "out", genbank_path=gb, fetch_genbank=False
    )
    assert written and written[0].stat().st_size > 0
    assert len(result.engineered_mapping_hits) == 1
    mapping = result.engineered_mapping_hits[0]
    assert any(f.label == "bla" for f in mapping.features)
    assert "bla" in mapping.alignments[0].overlapping_features
