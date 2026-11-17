"""Example PlasmidScreen library usage: screen → map-hits → optional viz.

After ``map_engineered_reads``, each ``ReadFlagDetail`` carries the same primary
BLAST fields written into the engineered report TSV:

    matched_plasmid, match_target_start/end, match_sequence_identity,
    match_score, overlapping_features
"""
from __future__ import annotations

from pathlib import Path

from plasmidScreen.api import (
    ScreenResult,
    map_engineered_reads,
    run_screen,
    visualize_engineered_mappings,
)

# Paths — adjust for your machine.
ROOT = Path(__file__).resolve().parents[1]
reads = ROOT / "example" / "data" / "example_pKT0174.fasta"
# Alternate: Path("/dodo/fmq1/pKT0174.fasta")
kraken_db = Path("/dodo/fmq1/PlasmidScreenDBNoAddgene")
# diamond_db = Path("/dodo/dbs/uniref_march_2025_with_tax.dmnd")
reference_db = ROOT / "Downloads" / "pKT0174_viz_ref.fasta"
report_path = ROOT / "Downloads" / "pKT0174_validate" / "api_example_report.txt"
viz_dir = ROOT / "Downloads" / "pKT0174_validate" / "api_example_viz"

# Dummy shape of one engineered read after map-hits (illustrative only).
DUMMY_READ_AFTER_MAP = {
    "read_id": "pKT0174",
    "overall_label": "Synthetic",
    "engineered_overall": True,
    "engineered_methods": ["engineered_kmer_scan"],
    "p_value": 3.1e-57,
    "eng_kmer_coverage": 0.0148,
    # Filled by map_engineered_reads (--map-hits):
    "matched_plasmid": "PZ096693.1",
    "match_target_start": 3473,
    "match_target_end": 6150,
    "match_sequence_identity": 99.888,
    "match_score": 4813.0,
    "overlapping_features": ("promoter", "ColE1 origin", "bla"),
}


def _print_screen_result(result: ScreenResult, *, title: str) -> None:
    print(f"\n=== {title} ===")
    print(f"engineered_report_path: {result.engineered_report_path}")
    print(
        f"counts: engineered={result.overall_synthetic_count} "
        f"natural={result.overall_natural_count}"
    )
    print(f"top_engineered_references: {result.top_engineered_references}")
    print("per_read:")
    for r in result.per_read:
        print(
            {
                "read_id": r.read_id,
                "overall_label": r.overall_label,
                "engineered_overall": r.engineered_overall,
                "engineered_methods": r.engineered_methods,
                "p_value": r.p_value,
                "eng_kmer_coverage": r.eng_kmer_coverage,
                "matched_plasmid": r.matched_plasmid,
                "match_target_start": r.match_target_start,
                "match_target_end": r.match_target_end,
                "match_sequence_identity": r.match_sequence_identity,
                "match_score": r.match_score,
                "overlapping_features": r.overlapping_features,
            }
        )


if __name__ == "__main__":
    print("Dummy ReadFlagDetail after map-hits (expected shape):")
    print(DUMMY_READ_AFTER_MAP)

    # 1) Screen (k-mer + optional codon). Match columns are still empty here.
    report_path.parent.mkdir(parents=True, exist_ok=True)
    screen_result: ScreenResult = run_screen(
        fasta_file=reads,
        seq_type=1,
        engineered_kmer_threshold=10,
        kraken_db=str(kraken_db),
        engineered_report_path=str(report_path),
        threads=8,
        quiet_mode=False,
        run_codon_usage=False,  # set True + diamond_db path to enable CAI
        # diamond_db="/dodo/dbs/uniref_march_2025_with_tax.dmnd",
        # codon_cai_threshold=0.6,
        memory_mode=0,
        background_rate=0.001,
    )
    _print_screen_result(screen_result, title="After run_screen (before map)")

    # 2) Map-hits: BLAST + GenBank overlap → fills report + ScreenResult fields.
    screen_result = map_engineered_reads(
        screen_result,
        query_fasta=reads,
        reference_db=reference_db,
        threads=8,
        top_n=5,
    )
    _print_screen_result(screen_result, title="After map_engineered_reads")
    if report_path.is_file():
        print("\nEngineered report TSV:")
        print(report_path.read_text(encoding="utf-8"))

    # 3) Optional viz: pileup PNGs only (report columns already filled by map).
    viz = visualize_engineered_mappings(
        screen_result,
        viz_dir,
        fetch_genbank=False,
    )
    print("\n=== After visualize_engineered_mappings ===")
    print(f"pileup_paths: {[str(p) for p in viz.pileup_paths]}")
    assert screen_result.visualization is viz
    print(f"visualization.output_dir: {viz.output_dir}")
