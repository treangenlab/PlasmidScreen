from __future__ import annotations

from plasmidScreen.lib.diamond_host_taxonomy import (
    DIAMOND_OUTFMT,
    hits_to_orf_intervals,
    infer_orfs_and_host_taxids,
    parse_diamond_tsv,
    resolve_diamond_lines,
    run_diamond_blastx,
)


def _diamond_line(
    qseqid: str,
    sseqid: str,
    qstart: str,
    qend: str,
    pident: str,
    length: str,
    evalue: str,
    bitscore: str,
    staxids: str,
    sscinames: str = "",
) -> str:
    return "\t".join(
        [
            qseqid,
            sseqid,
            qstart,
            qend,
            pident,
            length,
            evalue,
            bitscore,
            staxids,
            sscinames,
        ]
    )


def test_run_diamond_blastx_command_args() -> None:
    from unittest.mock import patch

    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.return_value.stdout = ""
        mock_run.return_value.returncode = 0
        run_diamond_blastx("reads.fa", "db.dmnd", threads=2)
    cmd = mock_run.call_args[0][0]
    assert "--block-size" in cmd
    idx = cmd.index("-f")
    assert cmd[idx + 1 : idx + 1 + len(DIAMOND_OUTFMT)] == list(DIAMOND_OUTFMT)
    assert "qstart" in DIAMOND_OUTFMT
    assert "qend" in DIAMOND_OUTFMT
    assert "-o" not in cmd


def test_run_diamond_blastx_writes_output_file(tmp_path: Path) -> None:
    from unittest.mock import patch

    out = tmp_path / "diamond.tsv"
    out.write_text(
        "read1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\t\n",
        encoding="utf-8",
    )
    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        lines = run_diamond_blastx("reads.fa", "db.dmnd", output_path=out)
    cmd = mock_run.call_args[0][0]
    assert "-o" in cmd
    assert str(out) in cmd
    assert len(lines) == 1


def test_resolve_diamond_lines_loads_file(tmp_path: Path) -> None:
    from plasmidScreen.lib.diamond_host_taxonomy import resolve_diamond_lines

    path = tmp_path / "diamond.tsv"
    path.write_text("read1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\t\n", encoding="utf-8")
    lines, used = resolve_diamond_lines(
        "reads.fa", None, run_diamond=False, output_path=path
    )
    assert used == path
    assert len(lines) == 1


def test_parse_diamond_tsv_staxids_and_coordinates() -> None:
    lines = [
        _diamond_line(
            "read1",
            "protA",
            "10",
            "50",
            "98.5",
            "120",
            "1e-50",
            "200",
            "562",
            "Escherichia coli",
        ),
        _diamond_line("read1", "protB", "55", "90", "90.0", "80", "1e-20", "150", "562", ""),
        _diamond_line("read2", "protC", "1", "60", "85.0", "60", "1e-10", "100", "9606", ""),
    ]
    hits = parse_diamond_tsv(lines)
    assert len(hits["read1"]) == 2
    assert hits["read1"][0].bitscore == 200.0
    assert hits["read1"][0].staxids == ("562",)
    assert hits["read1"][0].query_start == 9
    assert hits["read1"][0].query_end == 50
    assert hits["read1"][1].staxids == ("562",)
    assert hits["read2"][0].staxids == ("9606",)


def test_hits_to_orf_intervals_merge() -> None:
    lines = [
        _diamond_line("read1", "a", "10", "50", "99", "40", "1e-50", "200", "562", ""),
        _diamond_line("read1", "b", "55", "90", "98", "35", "1e-40", "180", "562", ""),
        _diamond_line("read1", "c", "200", "260", "90", "60", "1e-10", "80", "9606", ""),
    ]
    hits = parse_diamond_tsv(lines)["read1"]
    merged = hits_to_orf_intervals("read1", hits, merge_gap_bp=10, min_orf_len_bp=30)
    assert len(merged) == 2
    assert merged[0].start == 9
    assert merged[0].end == 90
    assert merged[1].start == 199
    assert merged[1].end == 260


def test_parse_diamond_tsv_skips_malformed_rows() -> None:
    lines = [
        "read1\tprotA\tbad\t50\t99\t33\t1e-5\t200\t562",
        _diamond_line("read1", "protA", "1", "99", "99", "33", "1e-5", "200", "562", ""),
    ]
    hits = parse_diamond_tsv(lines)
    assert len(hits["read1"]) == 1


def test_infer_orfs_and_host_taxids_majority_host() -> None:
    lines = [
        _diamond_line("read1", "a", "1", "120", "99", "120", "1e-50", "200", "562", ""),
        _diamond_line("read1", "b", "1", "100", "98", "100", "1e-40", "180", "562", ""),
        _diamond_line("read1", "c", "300", "360", "90", "50", "1e-10", "50", "9606", ""),
    ]
    orfs_by_read, host_by_read = infer_orfs_and_host_taxids(lines)
    assert host_by_read["read1"] == "562"
    assert len(orfs_by_read["read1"]) >= 1


def test_infer_orfs_and_host_taxids() -> None:
    lines = [
        _diamond_line("read1", "a", "1", "99", "99", "33", "1e-5", "200", "562", ""),
    ]
    orfs_by_read, host_by_read = infer_orfs_and_host_taxids(
        lines, min_orf_len_bp=30
    )
    assert len(orfs_by_read["read1"]) == 1
    assert orfs_by_read["read1"][0].length_bp == 99
    assert host_by_read["read1"] == "562"
