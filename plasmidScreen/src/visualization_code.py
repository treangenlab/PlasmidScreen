"""Reference pileup figures for engineered-read database matches."""
from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from Bio import SeqIO

from plasmidScreen.lib.models import (
    EngineeredMappingHit,
    MappedQueryInterval,
    QueryMatchConfidence,
    ReadFlagDetail,
    ReferenceFeatureInterval,
    ScreenResult,
)

logger = logging.getLogger(__name__)

_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")

# Soft plasmid-map palette (dna_features_viewer).
_FEATURE_COLORS = {
    "CDS": "#ffcc88",
    "gene": "#ffcc88",
    "promoter": "#f4a261",
    "terminator": "#e76f51",
    "rep_origin": "#cdb4db",
    "misc_feature": "#bde0fe",
    "regulatory": "#ffe08a",
    "protein_bind": "#b5e48c",
    "misc_binding": "#90e0ef",
    "primer_bind": "#a2d2ff",
    "oriT": "#cdb4db",
}
_DEFAULT_FEATURE_COLOR = "#ced4da"
_HIT_PLUS = "#4c78a8"
_HIT_MINUS = "#e45756"
_TEXT = "#1c1917"
_MUTED = "#6b7280"
# Above this many alignments, draw a coverage aggregate instead of one row each.
_PER_QUERY_MAX = 12
_ANNOTATION_TYPES = (
    "CDS",
    "gene",
    "promoter",
    "terminator",
    "rep_origin",
    "misc_feature",
    "regulatory",
    "protein_bind",
    "misc_binding",
    "primer_bind",
    "oriT",
)
# Prefer these when gene + CDS share a span.
_PRIMARY_FEATURE_TYPES = {"gene", "CDS", "promoter", "terminator", "rep_origin", "regulatory"}


@dataclass(frozen=True)
class PlasmidAnnotation:
    """One interval annotation on a reference plasmid (typically from GenBank)."""

    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    label: str
    feature_type: str
    strand: int = 1  # +1 / -1 / 0
    product: str | None = None


def group_hits_by_plasmid(
        reads: Iterable[ReadFlagDetail],
) -> dict[str, list[QueryMatchConfidence]]:
    """
    Aggregate primary alignments by reference accession.

    Uses ``database_match_summary`` only (same hit as the engineered report
    Plasmid / coords / identity / score / Overlapping_Features columns), not
    secondary HSPs in ``database_match_hits``.
    """
    by_plasmid: dict[str, list[QueryMatchConfidence]] = defaultdict(list)
    for read in reads:
        summary = read.database_match_summary
        if summary is None:
            continue
        by_plasmid[summary.accession_id].append(summary)
    return dict(by_plasmid)


def build_engineered_mapping_hits(
        reads: Iterable[ReadFlagDetail],
        *,
        annotations_by_accession: Mapping[str, Sequence[PlasmidAnnotation]] | None = None,
) -> tuple[EngineeredMappingHit, ...]:
    """
    Build per-plasmid engineered mapping hits that mirror pileup figure content.

    Each entry lists query alignment intervals on the reference, optional
    GenBank features, and coverage aggregates (covered bp / mean / max depth).
    """
    hits_by_plasmid = group_hits_by_plasmid(reads)
    summaries: list[EngineeredMappingHit] = []
    for accession_id, hits in sorted(hits_by_plasmid.items()):
        annots = _dedupe_annotations(
            _annotations_for_accession(accession_id, annotations_by_accession)
        )
        ref_len = _reference_length(hits, annots)
        ordered = sorted(
            hits,
            key=lambda h: (
                h.target_start if h.target_start is not None else 10**18,
                h.query_id or "",
                -(h.match_score or 0.0),
            ),
        )
        alignments: list[MappedQueryInterval] = []
        spans: list[tuple[int, int]] = []
        query_ids: list[str] = []
        seen_queries: set[str] = set()
        for hit in ordered:
            interval = _hit_interval(hit)
            if interval is None:
                continue
            start, end = interval
            spans.append((start, end))
            qid = hit.query_id or "?"
            if qid not in seen_queries:
                seen_queries.add(qid)
                query_ids.append(qid)
            alignments.append(
                MappedQueryInterval(
                    query_id=qid,
                    target_start=start,
                    target_end=end,
                    sequence_identity=hit.sequence_identity,
                    match_score=hit.match_score,
                    overlapping_features=tuple(
                        _genes_covered_by_span(annots, start, end)
                    ),
                )
            )
        if not alignments:
            continue
        cov = _coverage_profile(spans, ref_len=ref_len)
        covered_bp = sum(1 for depth in cov if depth > 0)
        mean_coverage = (sum(cov) / ref_len) if ref_len else 0.0
        max_coverage = max(cov) if cov else 0
        features = tuple(
            ReferenceFeatureInterval(
                label=_feature_display_label(ann) or ann.label,
                feature_type=ann.feature_type,
                start=ann.start,
                end=ann.end,
                strand=ann.strand,
            )
            for ann in annots
        )
        summaries.append(
            EngineeredMappingHit(
                accession_id=accession_id,
                target_length=ref_len,
                n_queries=len(query_ids),
                n_alignments=len(alignments),
                query_ids=tuple(query_ids),
                alignments=tuple(alignments),
                features=features,
                covered_bp=covered_bp,
                mean_coverage=mean_coverage,
                max_coverage=max_coverage,
            )
        )
    return tuple(summaries)


def _qualifier(feature, *keys: str) -> str | None:
    for key in keys:
        values = feature.qualifiers.get(key)
        if values:
            return str(values[0])
    return None


def _feature_label(feature) -> str:
    ftype = feature.type
    gene = _qualifier(feature, "gene", "locus_tag", "label")
    product = _qualifier(feature, "product", "note")
    bound = _qualifier(feature, "bound_moiety")
    reg_class = _qualifier(feature, "regulatory_class")

    if ftype in {"protein_bind", "misc_binding"} and bound:
        return bound if len(bound) <= 42 else bound[:39] + "…"
    if ftype == "regulatory":
        pretty = (reg_class or "regulatory").replace("_", " ")
        # GenBank often puts the associated gene in note, not gene=.
        tag = gene or (product if product and len(product) <= 16 else None)
        if tag and tag.lower() not in pretty.lower():
            return f"{pretty} ({tag})"
        return pretty
    if gene and product and product.lower() not in {gene.lower(), "n/a"}:
        short = product if len(product) <= 36 else product[:33] + "…"
        return f"{gene} ({short})"
    if gene:
        return gene
    if product:
        return product if len(product) <= 42 else product[:39] + "…"
    return ftype


def _with_richer_label(
        primary: PlasmidAnnotation,
        other: PlasmidAnnotation,
) -> PlasmidAnnotation:
    """Prefer gene name + CDS product/note when both occupy the same span."""
    label = primary.label
    product = primary.product or other.product
    if "(" not in label and other.label and other.label != label:
        # e.g. gene "Mob" + CDS "Mob (involved in mobilization)"
        if other.label.startswith(label):
            label = other.label
        elif product and product.lower() not in label.lower():
            short = product if len(product) <= 36 else product[:33] + "…"
            label = f"{label} ({short})"
    elif product and "(" not in label and product.lower() not in label.lower():
        short = product if len(product) <= 36 else product[:33] + "…"
        label = f"{label} ({short})"
    if label == primary.label and product == primary.product:
        return primary
    return PlasmidAnnotation(
        start=primary.start,
        end=primary.end,
        label=label,
        feature_type=primary.feature_type,
        strand=primary.strand,
        product=product,
    )


def _feature_parts(feature) -> list[tuple[int, int, int]]:
    """
    Expand compound / join / order locations into discrete intervals.

    Critical for origin-spanning genes (e.g. Mob on pBBR1) so we do not draw one
    bar across the whole plasmid.
    """
    strand = int(feature.location.strand or 0)
    parts = getattr(feature.location, "parts", None) or [feature.location]
    out: list[tuple[int, int, int]] = []
    for part in parts:
        start = int(part.start)
        end = int(part.end)
        if end <= start:
            continue
        part_strand = int(getattr(part, "strand", None) or strand or 0)
        out.append((start, end, part_strand))
    return out


def annotations_from_features(
        features,
        *,
        feature_types: Sequence[str] = _ANNOTATION_TYPES,
) -> list[PlasmidAnnotation]:
    """Convert BioPython SeqFeature objects into per-interval annotations."""
    wanted = {t.lower() for t in feature_types}
    annotations: list[PlasmidAnnotation] = []
    for feature in features:
        ftype = feature.type
        if ftype.lower() not in wanted or ftype == "source":
            continue
        label = _feature_label(feature)
        product = _qualifier(feature, "product", "note")
        for start, end, strand in _feature_parts(feature):
            annotations.append(
                PlasmidAnnotation(
                    start=start,
                    end=end,
                    label=label,
                    feature_type=ftype,
                    strand=strand,
                    product=product,
                )
            )
    return annotations


def load_genbank_annotations(
        genbank_path: str | Path,
        *,
        feature_types: Sequence[str] = _ANNOTATION_TYPES,
) -> dict[str, list[PlasmidAnnotation]]:
    """
    Parse GenBank/GenPept records into accession -> annotation intervals.

    Coordinates are converted to 0-based half-open to match blastn hit coords.
    Compound locations are exploded into parts. Records are indexed by ``id``,
    ``name``, and accession tokens so pileup accessions can resolve flexibly.
    """
    path = Path(genbank_path)
    by_key: dict[str, list[PlasmidAnnotation]] = defaultdict(list)

    for record in SeqIO.parse(str(path), "genbank"):
        annotations = annotations_from_features(
            record.features, feature_types=feature_types
        )
        if not annotations:
            continue
        keys = {record.id, record.name}
        for xref in record.annotations.get("accessions", []) or []:
            keys.add(str(xref))
        for key in list(keys):
            if key and "." in key:
                keys.add(key.split(".", 1)[0])
        for key in keys:
            if key:
                by_key[key].extend(annotations)

    logger.debug(
        "Loaded GenBank annotations from %s for %d accession key(s)",
        path,
        len(by_key),
    )
    return dict(by_key)


def _reference_length(
        hits: Sequence[QueryMatchConfidence],
        annotations: Sequence[PlasmidAnnotation] = (),
) -> int:
    lengths = [h.target_length for h in hits if h.target_length and h.target_length > 0]
    max_end = max((h.target_end or 0) for h in hits) if hits else 0
    if annotations:
        max_end = max(max_end, max(a.end for a in annotations))
    if lengths:
        return max(max(lengths), max_end)
    return max(max_end, 1)


def _hit_interval(hit: QueryMatchConfidence) -> tuple[int, int] | None:
    if hit.target_start is None or hit.target_end is None:
        return None
    start = min(hit.target_start, hit.target_end)
    end = max(hit.target_start, hit.target_end)
    if end <= start:
        end = start + 1
    return start, end


def _strand_of(hit: QueryMatchConfidence) -> str:
    strand = hit.functional_metadata.get("strand")
    return strand if strand in {"+", "-"} else "+"


def _dedupe_annotations(
        annotations: Sequence[PlasmidAnnotation],
) -> list[PlasmidAnnotation]:
    """Drop redundant CDS when a gene occupies the same span; keep other types."""
    by_span: dict[tuple[int, int, int], PlasmidAnnotation] = {}
    extras: list[PlasmidAnnotation] = []
    for ann in annotations:
        key = (ann.start, ann.end, ann.strand)
        if ann.feature_type in {"gene", "CDS"}:
            existing = by_span.get(key)
            if existing is None:
                by_span[key] = ann
            elif existing.feature_type == "CDS" and ann.feature_type == "gene":
                by_span[key] = _with_richer_label(ann, existing)
            elif existing.feature_type == "gene" and ann.feature_type == "CDS":
                by_span[key] = _with_richer_label(existing, ann)
            elif existing.feature_type == ann.feature_type:
                if len(ann.label) > len(existing.label):
                    by_span[key] = ann
        else:
            extras.append(ann)
    return sorted(list(by_span.values()) + extras, key=lambda a: (a.start, a.end, a.label))


def _short_accession(accession_id: str) -> str:
    if accession_id.startswith("gnl|uv|"):
        parts = accession_id.split("|")
        if len(parts) >= 3:
            return parts[2]
    return accession_id if len(accession_id) <= 48 else accession_id[:45] + "…"


def _coord_label(start: int, end: int) -> str:
    """1-based inclusive display coords for biologists."""
    return f"{start + 1:,}–{end:,}"


def _import_dna_features_viewer():
    try:
        from dna_features_viewer import GraphicFeature, GraphicRecord
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Pileup figures require dna_features_viewer. "
            "Install with: pip install 'plasmidScreen[viz]' "
            "or: pip install dna_features_viewer"
        ) from exc
    return GraphicFeature, GraphicRecord


def _feature_display_label(ann: PlasmidAnnotation) -> str | None:
    """Label genes / CDS / oris; leave tiny sites unlabeled."""
    if ann.feature_type not in {
        "gene",
        "CDS",
        "promoter",
        "terminator",
        "rep_origin",
        "regulatory",
        "oriT",
    }:
        return None
    label = ann.label
    # Prefer short gene names over long product notes.
    if "(" in label:
        label = label.split("(", 1)[0].strip()
    if len(label) > 28:
        label = label[:25] + "…"
    return label or None


def _hit_window(
        hits: Sequence[QueryMatchConfidence],
        ref_len: int,
) -> tuple[int, int] | None:
    spans = [_hit_interval(h) for h in hits]
    spans = [s for s in spans if s is not None]
    if not spans:
        return None
    start = min(s[0] for s in spans)
    end = max(s[1] for s in spans)
    # If the hit already covers most of the plasmid, show the full map.
    if (end - start) >= 0.7 * ref_len:
        return None
    pad = max(500, int((end - start) * 0.8))
    return max(0, start - pad), min(ref_len, end + pad)


def _genes_covered_by_span(
        annots: Sequence[PlasmidAnnotation],
        start: int,
        end: int,
) -> list[str]:
    names: list[str] = []
    for ann in annots:
        label = _feature_display_label(ann)
        if not label:
            continue
        overlap = max(0, min(ann.end, end) - max(ann.start, start))
        if overlap >= 0.45 * max(ann.end - ann.start, 1):
            names.append(label)
    return names


def _coverage_profile(
        spans: Sequence[tuple[int, int]],
        *,
        ref_len: int,
) -> list[int]:
    """Per-base count of how many query intervals cover each position."""
    if ref_len <= 0:
        return []
    delta = [0] * (ref_len + 1)
    for start, end in spans:
        s = max(0, min(start, ref_len))
        e = max(0, min(end, ref_len))
        if e <= s:
            continue
        delta[s] += 1
        delta[e] -= 1
    cov: list[int] = []
    running = 0
    for i in range(ref_len):
        running += delta[i]
        cov.append(running)
    return cov


def render_plasmid_pileup(
        accession_id: str,
        hits: Sequence[QueryMatchConfidence],
        *,
        annotations: Sequence[PlasmidAnnotation] = (),
        dpi: int = 180,
):
    """
    Draw GenBank features (top) and one row per primary query alignment (bottom).

    Keeping alignments off the annotation map avoids the crowded stacked-arrow
    look when many reads hit the same plasmid. Returns a Matplotlib figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    GraphicFeature, GraphicRecord = _import_dna_features_viewer()

    ordered = sorted(
        hits,
        key=lambda h: (
            h.target_start if h.target_start is not None else 10**18,
            h.query_id or "",
            -(h.match_score or 0.0),
        ),
    )
    annots = _dedupe_annotations(annotations)
    ref_len = _reference_length(ordered, annots)
    short = _short_accession(accession_id)
    # Multi-hit maps stay full-length so distant loci remain comparable.
    window = None if len(ordered) > 1 else _hit_window(ordered, ref_len)
    x0, x1 = (window if window is not None else (0, ref_len))

    gene_features = []
    for ann in annots:
        label = _feature_display_label(ann)
        if label is None:
            continue
        strand = ann.strand if ann.strand in (-1, 1) else 0
        gene_features.append(
            GraphicFeature(
                start=ann.start,
                end=ann.end,
                strand=strand,
                color=_FEATURE_COLORS.get(ann.feature_type, _DEFAULT_FEATURE_COLOR),
                label=label,
                thickness=12,
                box_linewidth=0.6,
                label_link_color="#9ca3af",
            )
        )

    gene_record = GraphicRecord(
        sequence_length=ref_len,
        features=gene_features,
        feature_level_height=0.95,
        labels_spacing=20,
        plots_indexing="genbank",
    )
    if window is not None:
        gene_record = gene_record.crop(window)

    hit_rows: list[tuple[str, int, int, str]] = []
    for hit in ordered:
        interval = _hit_interval(hit)
        if interval is None:
            continue
        start, end = interval
        qid = hit.query_id or "?"
        short_q = qid if len(qid) <= 28 else qid[:25] + "…"
        covered = _genes_covered_by_span(annots, start, end)
        note_parts: list[str] = []
        if covered:
            note_parts.append(covered[0])
        if hit.sequence_identity is not None:
            note_parts.append(f"{hit.sequence_identity:.0f}%")
        note = " · ".join(note_parts)
        hit_rows.append((short_q, start, end, note))

    aggregate = len(hit_rows) > _PER_QUERY_MAX
    n_hits = max(len(hit_rows), 1)
    gene_h = 2.4 if gene_features else 1.2
    if aggregate:
        pile_h = 2.0
    else:
        pile_h = max(1.4, 0.55 * n_hits + 0.6)
    fig_h = gene_h + pile_h + 0.55
    fig, (ax_genes, ax_pile) = plt.subplots(
        2,
        1,
        figsize=(13.0, fig_h),
        dpi=dpi,
        sharex=True,
        gridspec_kw={"height_ratios": [gene_h, pile_h], "hspace": 0.18},
    )
    fig.patch.set_facecolor("white")

    ax_genes.set_facecolor("white")
    gene_record.plot(
        ax=ax_genes,
        with_ruler=True,
        annotate_inline=False,
        elevate_outline_annotations=False,
        max_label_length=32,
        max_line_length=22,
    )
    ax_genes.set_xlim(x0, x1)
    # Crop excess whitespace DFV leaves above the feature lanes.
    y0, y1 = ax_genes.get_ylim()
    span = y1 - y0 if y1 != y0 else 1.0
    ax_genes.set_ylim(y0 - 0.05 * span, y1 - 0.12 * span)
    for spine in ("top", "right"):
        ax_genes.spines[spine].set_visible(False)
    ax_genes.set_xlabel("")

    ax_pile.set_facecolor("#fafafa")
    if not hit_rows:
        ax_pile.text(
            0.5,
            0.5,
            "No alignment intervals",
            transform=ax_pile.transAxes,
            ha="center",
            va="center",
            color=_MUTED,
        )
        ax_pile.set_yticks([])
        ax_pile.set_title("Queries", fontsize=10, loc="left", color=_MUTED, pad=6)
    elif aggregate:
        # Coverage depth + merged hotspot bars (no per-read IDs).
        spans = [(s, e) for _, s, e, _ in hit_rows]
        cov = _coverage_profile(spans, ref_len=ref_len)
        xs = list(range(ref_len))
        ax_pile.fill_between(
            xs,
            cov,
            step="pre",
            color=_HIT_PLUS,
            alpha=0.35,
            linewidth=0,
            zorder=2,
        )
        ax_pile.step(
            xs,
            cov,
            where="post",
            color=_HIT_PLUS,
            linewidth=1.2,
            zorder=3,
        )
        max_cov = max(cov) if cov else 1
        ax_pile.set_ylim(0, max_cov * 1.12)
        ax_pile.set_ylabel("# overlapping queries", fontsize=9, color=_MUTED)
        ax_pile.set_title(
            f"Query coverage ({len(hit_rows)} reads aggregated)",
            fontsize=10,
            loc="left",
            color=_MUTED,
            pad=6,
        )
    else:
        y_labels: list[str] = []
        for i, (qid, start, end, note) in enumerate(hit_rows):
            y = n_hits - 1 - i
            width = max(end - start, 1)
            ax_pile.add_patch(
                FancyBboxPatch(
                    (start, y - 0.28),
                    width,
                    0.56,
                    boxstyle="round,pad=0.01,rounding_size=2.0",
                    linewidth=0,
                    facecolor=_HIT_PLUS,
                    alpha=0.9,
                    zorder=2,
                )
            )
            y_labels.append(qid)
            if note:
                ax_pile.text(
                    end + 0.008 * max(x1 - x0, 1),
                    y,
                    note,
                    ha="left",
                    va="center",
                    fontsize=8,
                    color=_MUTED,
                    zorder=4,
                    clip_on=False,
                )
        ax_pile.set_yticks(list(range(n_hits)))
        ax_pile.set_yticklabels(list(reversed(y_labels)), fontsize=9)
        ax_pile.set_ylim(-0.7, n_hits - 0.3)
        ax_pile.set_title("Queries", fontsize=10, loc="left", color=_MUTED, pad=6)

    ax_pile.set_xlim(x0, x1)
    ax_pile.set_xlabel("Reference position (bp)", fontsize=10, color=_TEXT)
    ax_pile.grid(axis="x", color="#e5e7eb", linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax_pile.spines[spine].set_visible(False)
    ax_pile.spines["left"].set_color("#d1d5db")
    ax_pile.spines["bottom"].set_color("#d1d5db")

    if window is not None:
        region = (
            f"match region {window[0] + 1:,}–{window[1]:,} "
            f"of {ref_len:,} bp"
        )
    else:
        region = f"{ref_len:,} bp"
    title = f"{short}   ·   {len(ordered)} alignment(s)"
    if annots:
        title += f"   ·   {len(annots)} features"
    title += f"   ·   {region}"
    fig.suptitle(title, fontsize=12, color=_TEXT, fontweight="semibold", y=0.99)
    fig._plasmid_accession = accession_id  # type: ignore[attr-defined]
    fig.subplots_adjust(left=0.12, right=0.98, top=0.93, bottom=0.09, hspace=0.20)
    return fig



def _safe_filename(accession_id: str) -> str:
    cleaned = _SAFE_NAME.sub("_", accession_id).strip("._")
    return cleaned or "plasmid"


def _annotation_lookup_keys(accession_id: str) -> list[str]:
    """Candidate keys for matching GenBank annotations to a reference accession."""
    keys: list[str] = []
    seen: set[str] = set()

    def _add(value: str | None) -> None:
        if not value or value in seen:
            return
        seen.add(value)
        keys.append(value)

    _add(accession_id)
    head = accession_id.split("|", 1)[0]
    _add(head)
    _add(accession_id.split(".", 1)[0])
    _add(head.split(".", 1)[0])

    parts = accession_id.split("|")
    if len(parts) >= 3 and parts[0] == "gnl" and parts[1] == "uv":
        token = parts[2]
        _add(token)
        acc = token.split(":", 1)[0]
        _add(acc)
        _add(acc.split(".", 1)[0])

    return keys


def _annotations_for_accession(
        accession_id: str,
        annotations_by_accession: Mapping[str, Sequence[PlasmidAnnotation]] | None,
) -> list[PlasmidAnnotation]:
    if not annotations_by_accession:
        return []
    for key in _annotation_lookup_keys(accession_id):
        if key in annotations_by_accession:
            return list(annotations_by_accession[key])
    return []


def write_plasmid_pileups(
        hits_by_plasmid: Mapping[str, Sequence[QueryMatchConfidence]],
        output_dir: str | Path,
        *,
        annotations_by_accession: Mapping[str, Sequence[PlasmidAnnotation]] | None = None,
        dpi: int = 180,
) -> list[Path]:
    """Write one PNG pileup figure per plasmid accession."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for accession_id, hits in sorted(hits_by_plasmid.items()):
        annots = _annotations_for_accession(accession_id, annotations_by_accession)
        fig = render_plasmid_pileup(
            accession_id,
            hits,
            annotations=annots,
            dpi=dpi,
        )
        path = out_dir / f"pileup_{_safe_filename(accession_id)}.png"
        import matplotlib.pyplot as plt

        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        written.append(path)
        logger.debug(
            "Wrote pileup for %s (%d hits, %d features) -> %s",
            accession_id,
            len(hits),
            len(annots),
            path,
        )
    return written


def _load_annotations_for_accessions(
        accession_ids: Iterable[str],
        *,
        genbank_path: str | Path | None = None,
        annotations_by_accession: Mapping[str, Sequence[PlasmidAnnotation]] | None = None,
        fetch_genbank: bool = False,
        genbank_cache_dir: str | Path | None = None,
        ncbi_api_key: str | None = None,
) -> dict[str, list[PlasmidAnnotation]]:
    """Merge optional GenBank file + local cache (+ optional NCBI fetch)."""
    ann_map: dict[str, list[PlasmidAnnotation]] = {}
    if genbank_path is not None:
        ann_map.update(
            {key: list(vals) for key, vals in load_genbank_annotations(genbank_path).items()}
        )
    if annotations_by_accession:
        for key, vals in annotations_by_accession.items():
            ann_map[key] = list(vals)

    from plasmidScreen.lib.genbank_fetch import annotations_for_hit_ids

    cached_or_fetched = annotations_for_hit_ids(
        accession_ids,
        cache_dir=genbank_cache_dir,
        api_key=ncbi_api_key,
        fetch=fetch_genbank,
    )
    for key, vals in cached_or_fetched.items():
        ann_map.setdefault(key, list(vals))
    return ann_map


def annotate_screen_result_overlapping_features(
        screen_result: ScreenResult,
        *,
        genbank_path: str | Path | None = None,
        annotations_by_accession: Mapping[str, Sequence[PlasmidAnnotation]] | None = None,
        fetch_genbank: bool = False,
        genbank_cache_dir: str | Path | None = None,
        ncbi_api_key: str | None = None,
) -> ScreenResult:
    """
    Fill ``ReadFlagDetail.overlapping_features`` from the primary BLAST interval.

    Uses the local GenBank cache by default (same source as pileup tracks).
    Call after :func:`map_engineered_reads` enrichment so report columns include
    features without requiring visualization.
    """
    accessions = {
        detail.matched_plasmid
        for detail in screen_result.per_read
        if detail.matched_plasmid
    }
    if not accessions:
        return screen_result

    ann_map = _load_annotations_for_accessions(
        accessions,
        genbank_path=genbank_path,
        annotations_by_accession=annotations_by_accession,
        fetch_genbank=fetch_genbank,
        genbank_cache_dir=genbank_cache_dir,
        ncbi_api_key=ncbi_api_key,
    )
    if not ann_map:
        logger.debug(
            "No GenBank annotations available for overlapping features "
            "(%d accession(s)); leave Overlapping_Features empty.",
            len(accessions),
        )
        return screen_result

    updated: list[ReadFlagDetail] = []
    changed = False
    for detail in screen_result.per_read:
        summary = detail.database_match_summary
        if summary is None or summary.target_start is None or summary.target_end is None:
            updated.append(detail)
            continue
        start = min(summary.target_start, summary.target_end)
        end = max(summary.target_start, summary.target_end)
        annots = _dedupe_annotations(
            _annotations_for_accession(summary.accession_id, ann_map)
        )
        feats = tuple(_genes_covered_by_span(annots, start, end))
        if feats != detail.overlapping_features:
            updated.append(replace(detail, overlapping_features=feats))
            changed = True
        else:
            updated.append(detail)
    if changed:
        screen_result.per_read = updated
    return screen_result


def create_pile_ups(
        screen_result: ScreenResult,
        visualize_path: str | Path,
        *,
        genbank_path: str | Path | None = None,
        annotations_by_accession: Mapping[str, Sequence[PlasmidAnnotation]] | None = None,
        dpi: int = 180,
        engineered_only: bool = True,
        fetch_genbank: bool = False,
        genbank_cache_dir: str | Path | None = None,
        ncbi_api_key: str | None = None,
) -> list[Path]:
    """
    Write one PNG pileup per plasmid using each read's primary BLAST hit.

    Intervals match the engineered report mapping columns (primary accession /
    coords). Optional ``genbank_path`` / ``annotations_by_accession`` overlays
    reference features above the alignment tracks. By default, NCBI features
    are read only from the local GenBank cache (fill it with
    ``prefetch_genbank_from_fasta`` / ``plasmidscreen prefetch-genbank``). Set
    ``fetch_genbank=True`` only if you intentionally want on-demand NCBI
    downloads during viz.
    """
    reads = (
        screen_result.per_read_engineered
        if engineered_only
        else screen_result.per_read
    )
    hits_by_plasmid = group_hits_by_plasmid(reads)
    # Prefer the sample-wide BLAST top references when enrichment already ranked them.
    top_refs = {
        r.accession_id for r in screen_result.top_engineered_references
    }
    if top_refs:
        hits_by_plasmid = {
            acc: hits
            for acc, hits in hits_by_plasmid.items()
            if acc in top_refs
        }
    if not hits_by_plasmid:
        logger.info("No primary database matches available; skipping pileups.")
        return []

    ann_map = _load_annotations_for_accessions(
        hits_by_plasmid.keys(),
        genbank_path=genbank_path,
        annotations_by_accession=annotations_by_accession,
        fetch_genbank=fetch_genbank,
        genbank_cache_dir=genbank_cache_dir,
        ncbi_api_key=ncbi_api_key,
    )

    # Attach mapping summaries used by figures / ScreenResult.visualization.
    screen_result.engineered_mapping_hits = build_engineered_mapping_hits(
        reads,
        annotations_by_accession=ann_map or None,
    )
    if top_refs:
        screen_result.engineered_mapping_hits = tuple(
            m
            for m in screen_result.engineered_mapping_hits
            if m.accession_id in top_refs
        )

    return write_plasmid_pileups(
        hits_by_plasmid,
        visualize_path,
        annotations_by_accession=ann_map or None,
        dpi=dpi,
    )
