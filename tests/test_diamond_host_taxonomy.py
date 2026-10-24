from __future__ import annotations

from pathlib import Path

from plasmidScreen.lib.diamond_host_taxonomy import (
    DIAMOND_OUTFMT,
    diamond_technology_flags,
    hits_to_orf_intervals,
    infer_orfs_and_host_taxids,
    parse_diamond_tsv,
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


def _touch_query_and_db(tmp_path: Path) -> tuple[Path, Path]:
    reads = tmp_path / "reads.fa"
    reads.write_text(">a\nACGT\n", encoding="utf-8")
    db = tmp_path / "db.dmnd"
    db.write_bytes(b"")
    return reads, db


def test_run_diamond_blastx_command_args(tmp_path: Path) -> None:
    from unittest.mock import patch

    reads, db = _touch_query_and_db(tmp_path)
    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.return_value.stdout = ""
        mock_run.return_value.returncode = 0
        run_diamond_blastx(reads, db, threads=2)
    cmd = mock_run.call_args[0][0]
    assert "--block-size" in cmd
    idx = cmd.index("-f")
    assert cmd[idx + 1 : idx + 1 + len(DIAMOND_OUTFMT)] == list(DIAMOND_OUTFMT)
    assert "qstart" in DIAMOND_OUTFMT
    assert "qend" in DIAMOND_OUTFMT
    assert "-o" not in cmd
    # Default seq technology is short-read (1)
    assert "--more-sensitive" in cmd
    assert "--frameshift" not in cmd


def test_run_diamond_blastx_seq_type_flags(tmp_path: Path) -> None:
    from unittest.mock import patch

    reads, db = _touch_query_and_db(tmp_path)
    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.return_value.stdout = ""
        mock_run.return_value.returncode = 0
        run_diamond_blastx(reads, db, sequencing_technology=0)
        long_cmd = mock_run.call_args[0][0]
        run_diamond_blastx(reads, db, sequencing_technology=1)
        short_cmd = mock_run.call_args[0][0]

    assert " ".join(diamond_technology_flags(0)) in " ".join(long_cmd)
    assert " ".join(diamond_technology_flags(1)) in " ".join(short_cmd)
    assert "--frameshift" in long_cmd and "--frameshift" not in short_cmd
    assert "--more-sensitive" in short_cmd and "--more-sensitive" not in long_cmd


def test_run_diamond_blastx_mem_mode_block_sizes(tmp_path: Path) -> None:
    from unittest.mock import patch

    reads, db = _touch_query_and_db(tmp_path)
    expected = {0: ("4", "3"), 1: ("6", "1"), 2: ("12", "1")}
    for mode, (block, chunks) in expected.items():
        with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
            mock_run.return_value.stdout = ""
            mock_run.return_value.returncode = 0
            run_diamond_blastx(reads, db, threads=1, mem_mode=mode)
        cmd = mock_run.call_args[0][0]
        assert cmd[cmd.index("--block-size") + 1] == block
        assert cmd[cmd.index("--index-chunks") + 1] == chunks


def test_run_diamond_blastx_invalid_mem_mode() -> None:
    import pytest

    with pytest.raises(ValueError, match="memory_mode"):
        run_diamond_blastx("reads.fa", "db.dmnd", threads=1, mem_mode=9)


def test_run_diamond_blastx_writes_output_file(tmp_path: Path) -> None:
    from unittest.mock import patch

    reads, db = _touch_query_and_db(tmp_path)
    out = tmp_path / "diamond.tsv"
    out.write_text(
        "read1\tref\t1\t99\t99.0\t33\t1e-5\t200\t562\t\n",
        encoding="utf-8",
    )
    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        lines = run_diamond_blastx(reads, db, output_path=out)
    cmd = mock_run.call_args[0][0]
    assert "-o" in cmd
    assert str(out) in cmd
    assert len(lines) == 1


def test_run_diamond_blastx_wraps_subprocess_failure(tmp_path: Path) -> None:
    import subprocess
    from unittest.mock import patch

    import pytest

    from plasmidScreen.lib.exceptions import DiamondError

    reads, db = _touch_query_and_db(tmp_path)
    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.side_effect = subprocess.CalledProcessError(
            1, "diamond", stderr="boom"
        )
        with pytest.raises(DiamondError, match="exit code 1"):
            run_diamond_blastx(reads, db)


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
