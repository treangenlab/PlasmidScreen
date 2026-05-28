from __future__ import annotations

from plasmidScreen.lib.diamond_host_taxonomy import (
    hits_to_orf_intervals,
    infer_orfs_and_host_taxids,
    lca_taxid,
    parse_diamond_tsv,
    run_diamond_blastx,
)


def test_run_diamond_blastx_command_args() -> None:
    from unittest.mock import patch

    with patch("plasmidScreen.lib.diamond_host_taxonomy.subprocess.run") as mock_run:
        mock_run.return_value.stdout = ""
        mock_run.return_value.returncode = 0
        run_diamond_blastx("reads.fa", "db.dmnd", threads=2)
    cmd = mock_run.call_args[0][0]
    assert "--block-size" in cmd
    idx = cmd.index("--block-size")
    assert cmd[idx + 1] == "200"
    assert cmd[idx + 2] == "--index-chunks"


def test_parse_diamond_tsv_and_merge_intervals() -> None:
    lines = [
        "read1\t10\t50\t100\t562",
        "read1\t55\t90\t90\t562",
        "read1\t200\t260\t80\t9606",
    ]
    hits = parse_diamond_tsv(lines)["read1"]
    merged = hits_to_orf_intervals("read1", hits, merge_gap_bp=10, min_orf_len_bp=30)
    assert len(merged) == 2
    (s1, e1, _h1) = merged[0]
    (s2, e2, _h2) = merged[1]
    assert (s1, e1) == (10, 90)
    assert (s2, e2) == (200, 260)


def test_lca_taxid() -> None:
    # Simple taxonomy chain: 562 -> 561 -> 543
    parents = {"562": "561", "561": "543", "543": "543"}
    assert lca_taxid(["562", "562"], parents) == "562"
    assert lca_taxid(["562", "561"], parents) == "561"


def test_infer_orfs_and_host_taxids_majority() -> None:
    parents = {"562": "561", "561": "543", "543": "543", "9606": "9605", "9605": "314295", "314295": "314295"}
    lines = [
        "read1\t1\t120\t200\t562",
        "read1\t130\t260\t180\t562",
        "read1\t300\t360\t50\t9606",
    ]
    orfs_by_read, host_by_read = infer_orfs_and_host_taxids(lines, taxonomy_parents=parents)
    assert "read1" in orfs_by_read
    assert host_by_read["read1"] == "562"
