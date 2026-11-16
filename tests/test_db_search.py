"""Tests for blastn reference-hit selection and ScreenResult enrichment."""
from __future__ import annotations

from dataclasses import replace as dc_replace
from pathlib import Path
from unittest.mock import patch

from plasmidScreen.lib.db_search import (
    BlastnReferenceSearchEngine,
    HitFilterConfig,
    ReferenceSearchConfig,
    build_query_matches,
    enrich_screen_result,
    filter_significant_hits,
    parse_blast_outfmt6_line,
    rank_sample_references,
    write_engineered_query_fasta,
)
from plasmidScreen.lib.models import (
    NO_DATABASE_MATCH_NOTE,
    CodonAdaptationResult,
    EngineeredScanResult,
    QueryMatchConfidence,
    ReadEngineeringLabel,
    ReadFlagDetail,
    ScreenResult,
)


def _hit(
    accession: str,
    *,
    bitscore: float,
    evalue: float | None = 1e-20,
    pident: float = 90.0,
    query_id: str = "read1",
    query_start: int | None = None,
    query_end: int | None = None,
) -> QueryMatchConfidence:
    return QueryMatchConfidence(
        accession_id=accession,
        database_source="test_db",
        match_score=bitscore,
        sequence_identity=pident,
        query_id=query_id,
        evalue=evalue,
        bitscore=bitscore,
        functional_metadata={"alignment_tool": "blastn", "strand": "+"},
        query_start=query_start,
        query_end=query_end,
    )


def _blast_line(
    qid: str,
    sid: str,
    *,
    pident: float = 95.0,
    length: int = 100,
    qstart: int = 1,
    qend: int = 100,
    sstart: int = 10,
    send: int = 109,
    evalue: float = 1e-40,
    bitscore: float = 180.0,
) -> str:
    return "\t".join(
        [
            qid,
            sid,
            str(pident),
            str(length),
            "0",
            "0",
            str(qstart),
            str(qend),
            str(sstart),
            str(send),
            str(evalue),
            str(bitscore),
        ]
    )


def test_filter_significant_hits_by_bitscore() -> None:
    hits = [
        _hit("A", bitscore=40),
        _hit("B", bitscore=100),
        _hit("C", bitscore=50),
    ]
    kept = filter_significant_hits(hits, HitFilterConfig(min_bitscore=50.0))
    assert [h.accession_id for h in kept] == ["B", "C"]


def test_rank_sample_references_by_read_coverage() -> None:
    hits_by_query = {
        "r1": [_hit("ACC_A", bitscore=200, query_id="r1"), _hit("ACC_B", bitscore=80, query_id="r1")],
        "r2": [_hit("ACC_A", bitscore=150, query_id="r2")],
        "r3": [_hit("ACC_C", bitscore=300, query_id="r3")],
        "r4": [_hit("ACC_B", bitscore=90, query_id="r4"), _hit("ACC_C", bitscore=100, query_id="r4")],
    }
    ranks = rank_sample_references(hits_by_query, top_n=2)
    # All three cover 2 reads; bitscore sum breaks ties: C=400, A=350, B=170.
    assert [r.accession_id for r in ranks] == ["ACC_C", "ACC_A"]
    assert ranks[0].n_reads == 2
    assert ranks[0].sum_bitscore == 400.0


def test_build_query_matches_keeps_only_top_references() -> None:
    # TOP1 covers 3 reads, TOP2 covers 2, LOW covers 1 → top_n=2 keeps TOP1/TOP2.
    hits_by_query = {
        "r1": [
            _hit("TOP1", bitscore=200, query_id="r1", query_start=0, query_end=50),
            _hit("LOW", bitscore=190, query_id="r1", query_start=0, query_end=40),
        ],
        "r2": [
            _hit("TOP1", bitscore=180, query_id="r2", query_start=0, query_end=50),
            _hit("TOP2", bitscore=170, query_id="r2", query_start=0, query_end=50),
        ],
        "r3": [
            _hit("TOP1", bitscore=160, query_id="r3", query_start=0, query_end=50),
            _hit("TOP2", bitscore=155, query_id="r3", query_start=0, query_end=50),
        ],
    }
    result = build_query_matches(
        hits_by_query,
        HitFilterConfig(min_bitscore=50.0, top_n_references=2),
    )
    assert [r.accession_id for r in result.top_references] == ["TOP1", "TOP2"]
    assert {h.accession_id for h in result.matches["r1"].hits} == {"TOP1"}
    assert {h.accession_id for h in result.matches["r3"].hits} == {"TOP1", "TOP2"}
    assert "LOW" not in {
        h.accession_id for m in result.matches.values() for h in m.hits
    }


def test_parse_blast_outfmt6_line_coords() -> None:
    line = _blast_line(
        "read1",
        "REF_ACC",
        qstart=5,
        qend=104,
        sstart=20,
        send=10,  # reverse subject
        bitscore=210.0,
    )
    hit = parse_blast_outfmt6_line(line, database_source="refs:v1")
    assert hit is not None
    assert hit.accession_id == "REF_ACC"
    assert hit.query_id == "read1"
    assert hit.query_start == 4  # 0-based
    assert hit.query_end == 104
    assert hit.target_start == 9
    assert hit.target_end == 20
    assert hit.bitscore == 210.0
    assert hit.match_score == 210.0
    assert hit.functional_metadata["strand"] == "-"
    assert hit.functional_metadata["alignment_tool"] == "blastn"


def _minimal_screen_result(
    read_ids: list[str], *, engineered: set[str] | None = None
) -> ScreenResult:
    engineered = engineered or set()
    labels = []
    per_read = []
    for rid in read_ids:
        is_eng = rid in engineered
        label = "Synthetic" if is_eng else "Natural"
        labels.append(
            ReadEngineeringLabel(
                read_id=rid, label=label, eng_kmer_coverage=1.0 if is_eng else 0.0
            )
        )
        per_read.append(
            ReadFlagDetail(
                read_id=rid,
                kmer_label=label,  # type: ignore[arg-type]
                engineered_by_kmer_scan=is_eng,
                engineered_overall=is_eng,
                overall_label=label,  # type: ignore[arg-type]
            )
        )
    return ScreenResult(
        engineered_scan=EngineeredScanResult(
            labels=labels,
            synthetic_count=sum(1 for r in labels if r.label == "Synthetic"),
            natural_count=sum(1 for r in labels if r.label == "Natural"),
        ),
        codon_adaptation=CodonAdaptationResult(labels=[]),
        engineered_kmer_threshold=10,
        engineered_kmer_window_size=200,
        codon_cai_engineered_threshold=0.6,
        per_read=per_read,
    )


def test_enrich_screen_result_tags_no_match() -> None:
    result = _minimal_screen_result(
        ["read_a", "read_b"], engineered={"read_a", "read_b"}
    )
    hit_a = dc_replace(
        _hit("ACC1", bitscore=100, query_id="read_a"),
        query_start=0,
        query_end=100,
        target_start=10,
        target_end=110,
        target_length=1000,
    )
    enriched = enrich_screen_result(
        result,
        {"read_a": [hit_a]},
        database_source="test_db",
        annotate_only_ids={"read_a", "read_b"},
    )

    by_id = {r.read_id: r for r in enriched.per_read}
    assert by_id["read_a"].has_database_match is True
    assert by_id["read_a"].database_match_summary is not None
    assert by_id["read_a"].database_match_summary.accession_id == "ACC1"
    assert by_id["read_b"].has_database_match is False
    assert by_id["read_b"].database_match_note == NO_DATABASE_MATCH_NOTE
    assert by_id["read_b"].engineered_overall is True
    assert len(enriched.top_engineered_references) == 1
    assert enriched.top_engineered_references[0].accession_id == "ACC1"
    # Mapping summaries are viz-only (not filled by enrichment).
    assert enriched.engineered_mapping_hits == ()


def test_write_engineered_query_fasta(sample_fasta: Path, tmp_path: Path) -> None:
    out = write_engineered_query_fasta(
        {"read_synthetic_1"}, sample_fasta, output_path=tmp_path / "eng.fa"
    )
    assert out is not None
    text = out.read_text(encoding="utf-8")
    assert ">read_synthetic_1" in text
    assert "read_natural_1" not in text


def test_blast_engine_search_from_blast_lines(tmp_path: Path) -> None:
    lines = [
        _blast_line("read1", "ACC_HIT", bitscore=200.0),
        _blast_line("read1", "ACC_WEAK", bitscore=20.0),
        _blast_line("read2", "ACC_OTHER", bitscore=150.0),
    ]
    config = ReferenceSearchConfig(
        reference_db=tmp_path / "unused.fa",
        database_source="fixture_db",
        filters=HitFilterConfig(min_bitscore=50.0, top_n_references=5),
        run_blast=False,
    )
    engine = BlastnReferenceSearchEngine(config)
    out = engine.search_from_blast_lines(
        lines,
        query_ids=["read1", "read2", "read3"],
        query_lengths={"read1": 100, "read2": 100, "read3": 50},
    )

    assert list(out.matches.keys()) == ["read1", "read2", "read3"]
    assert len(out.matches["read1"].hits) == 1
    assert out.matches["read1"].hits[0].accession_id == "ACC_HIT"
    assert out.matches["read2"].hits[0].accession_id == "ACC_OTHER"
    assert out.matches["read3"].hits == ()
    assert {r.accession_id for r in out.top_references} == {"ACC_HIT", "ACC_OTHER"}


@patch("plasmidScreen.lib.db_search.subprocess.run")
@patch("plasmidScreen.lib.db_search.read_query_lengths", autospec=True)
def test_blast_engine_search_batch_invokes_blastn(
    mock_lengths, mock_run, tmp_path: Path, sample_fasta: Path
) -> None:
    mock_lengths.return_value = {"read_synthetic_1": 100}
    mock_run.return_value.stdout = (
        _blast_line("read_synthetic_1", "REF1", bitscore=220.0) + "\n"
    )
    mock_run.return_value.returncode = 0
    ref = tmp_path / "ref.fa"
    ref.write_text(">REF1\nACGT\n", encoding="utf-8")
    config = ReferenceSearchConfig(
        reference_db=ref,
        database_source="plasmid_refs",
        threads=2,
        filters=HitFilterConfig(min_bitscore=50.0, top_n_references=3),
        run_blast=True,
    )
    engine = BlastnReferenceSearchEngine(config)
    out = engine.search_batch(
        sample_fasta, query_ids=["read_synthetic_1", "read_natural_1"]
    )
    mock_run.assert_called_once()
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "blastn"
    assert "-subject" in cmd
    assert "-outfmt" in cmd
    assert out.matches["read_synthetic_1"].hits[0].accession_id == "REF1"
    assert out.matches["read_natural_1"].hits == ()
