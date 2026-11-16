"""Engineered-call TSV report read/write helpers."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

from plasmidScreen.lib.models import ReadFlagDetail, ScreenResult

logger = logging.getLogger(__name__)

REPORT_COLUMNS = (
    "Label",
    "Read_ID",
    "Methods",
    "P-Value",
    "Eng-Coverage Percentage",
    "Plasmid",
    "Target_Start",
    "Target_End",
    "Sequence_Identity",
    "Match_Score",
    "Overlapping_Features",
)

_NA = "NA"


def match_fields_for_read(detail: ReadFlagDetail) -> tuple[str, str, str, str, str, str]:
    """
    Format primary BLAST-hit columns for one report row.

    Reads ``matched_plasmid`` / coords / identity / score from
    ``database_match_summary`` and ``overlapping_features`` from the detail
    (filled by :func:`map_engineered_reads`).
    """
    if detail.database_match_summary is None:
        return (_NA, _NA, _NA, _NA, _NA, _NA)

    plasmid = detail.matched_plasmid or _NA
    start = (
        _NA if detail.match_target_start is None else str(detail.match_target_start)
    )
    end = _NA if detail.match_target_end is None else str(detail.match_target_end)
    identity = (
        _NA
        if detail.match_sequence_identity is None
        else f"{detail.match_sequence_identity:.3f}".rstrip("0").rstrip(".")
    )
    score = (
        _NA if detail.match_score is None else f"{detail.match_score:g}"
    )
    features = (
        ",".join(detail.overlapping_features)
        if detail.overlapping_features
        else _NA
    )
    return plasmid, start, end, identity, score, features


def format_engineered_report_row(detail: ReadFlagDetail) -> str:
    """One TSV data line (no trailing newline) for ``detail``."""
    label = "Engineered" if detail.engineered_overall else "Natural"
    methods = detail.engineered_methods
    methods_cell = ",".join(methods) if methods else _NA
    p_value = str(detail.p_value)
    coverage = (
        detail.eng_kmer_coverage if detail.eng_kmer_coverage is not None else 0.0
    )
    plasmid, start, end, identity, score, features = match_fields_for_read(detail)
    cells = (
        label,
        detail.read_id,
        methods_cell,
        p_value,
        f"{coverage:.2%}",
        plasmid,
        start,
        end,
        identity,
        score,
        features,
    )
    return "\t".join(cells)


def write_engineered_report(
        path: str | Path,
        per_read: Iterable[ReadFlagDetail],
) -> Path:
    """
    Write the engineered-call TSV (including optional primary BLAST columns).

    Mapping columns (including ``Overlapping_Features``) are ``NA`` until
    :func:`map_engineered_reads` has attached hits and GenBank annotations.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(REPORT_COLUMNS)]
    for detail in per_read:
        lines.append(format_engineered_report_row(detail))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("Wrote engineered report to %s", out)
    return out


def rewrite_engineered_report_if_present(result: ScreenResult) -> Path | None:
    """Rewrite ``result.engineered_report_path`` from current per-read state."""
    path = result.engineered_report_path
    if path is None:
        return None
    return write_engineered_report(path, result.per_read)
