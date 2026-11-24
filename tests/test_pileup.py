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


def test_group_hits_by_plasmid_uses_primary_only() -> None:
    result = _screen_with_hits(
        {
            "read_a": (
                # Primary (first) is plasmid_X; secondary HSP must not appear.
                _hit("plasmid_X", query_id="read_a", target_start=10, target_end=200),
                _hit("plasmid_Y", query_id="read_a", target_start=0, target_end=50),
                _hit(
                    "plasmid_X",
                    query_id="read_a",
                    target_start=500,
                    target_end=600,
                    match_score=10.0,
                ),
            ),
            "read_b": (
                _hit("plasmid_X", query_id="read_b", target_start=150, target_end=400),
            ),
        }
    )
    grouped = group_hits_by_plasmid(result.per_read_engineered)
    assert set(grouped) == {"plasmid_X"}
    assert {h.query_id for h in grouped["plasmid_X"]} == {"read_a", "read_b"}
    assert len(grouped["plasmid_X"]) == 2
    assert grouped["plasmid_X"][0].target_start == 10
    assert grouped["plasmid_X"][0].target_end == 200


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


def test_visualize_api_respects_top_refs_and_attaches_visualization(
    tmp_path: Path,
) -> None:
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
        result, tmp_path / "viz", fetch_genbank=False
    )
    assert not (tmp_path / "viz" / "visualization_report.json").exists()
    assert {m.accession_id for m in viz.engineered_mapping_hits} == {"plasmid_X"}
    assert all(p.name.startswith("pileup_plasmid_X") for p in viz.pileup_paths)
    assert result.engineered_mapping_hits == viz.engineered_mapping_hits
    assert result.visualization is viz
    by_id = {r.read_id: r for r in result.per_read}
    assert by_id["read_a"].matched_plasmid == "plasmid_X"
    assert by_id["read_a"].match_target_start == 0
    assert by_id["read_a"].match_target_end == 100


def test_render_plasmid_pileup_builds_figure() -> None:
    from plasmidScreen.src.visualization_code import PlasmidAnnotation

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
        annotations=(
            PlasmidAnnotation(
                start=50, end=200, label="bla", feature_type="CDS", strand=1
            ),
            PlasmidAnnotation(
                start=300, end=400, label="ori", feature_type="rep_origin", strand=1
            ),
        ),
    )
    assert fig._suptitle is not None
    assert "plasmid_X" in fig._suptitle.get_text()
    legend = fig.axes[0].get_legend()
    assert legend is not None
    labels = [t.get_text() for t in legend.get_texts()]
    assert "CDS / gene" in labels
    assert "origin of replication" in labels
    assert "query alignment" in labels
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


def test_dedupe_annotations_collapses_repeated_aliases() -> None:
    from plasmidScreen.src.visualization_code import (
        PlasmidAnnotation,
        _dedupe_annotations,
    )

    base = (
        PlasmidAnnotation(start=10, end=100, label="bla", feature_type="gene", strand=1),
        PlasmidAnnotation(start=10, end=100, label="bla", feature_type="CDS", strand=1),
        PlasmidAnnotation(
            start=200, end=210, label="promoter", feature_type="promoter", strand=1
        ),
        # Near-duplicate RBS (1 bp difference) should collapse.
        PlasmidAnnotation(
            start=1904,
            end=1909,
            label="ribosome binding site",
            feature_type="regulatory",
            strand=1,
        ),
        PlasmidAnnotation(
            start=1904,
            end=1910,
            label="ribosome binding site",
            feature_type="regulatory",
            strand=1,
        ),
    )
    # Simulate concatenating the same record under many GenBank accession aliases.
    duplicated = base * 14
    deduped = _dedupe_annotations(duplicated)
    assert len(deduped) == 3  # gene+CDS, promoter, one RBS
    types = {a.feature_type for a in deduped}
    assert "promoter" in types
    assert "regulatory" in types
    assert ("gene" in types) or ("CDS" in types)
    rbs = [a for a in deduped if a.feature_type == "regulatory"]
    assert len(rbs) == 1
    assert rbs[0].end - rbs[0].start == 6  # kept the longer 1904–1910


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
