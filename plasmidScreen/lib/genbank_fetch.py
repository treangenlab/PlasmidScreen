"""
On-demand NCBI GenBank downloads for pileup reference annotations.

Only resolves IDs that look like NCBI nucleotide accessions (plain ``MT024802.1``
or embedded UniVec tokens like ``gnl|uv|X66730.1:1-2687-49|…``). IMGPR / custom
local IDs are skipped. Requests are batched and paced to stay within NCBI's
unauthenticated eutils guidance (~≤3 requests/sec).
"""
from __future__ import annotations

import logging
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Iterable, Sequence

from Bio import SeqIO

from plasmidScreen.lib.funcs import get_default_db_path

logger = logging.getLogger(__name__)

# Plain NCBI nucleotide accessions (RefSeq / GenBank style).
_NCBI_ACC = re.compile(
    r"^(?P<acc>[A-Z]{1,2}_?[0-9]{5,}(?:\.[0-9]+)?)$"
)
# UniVec / Kraken-tagged: gnl|uv|ACCESSION:start-end…|kraken:taxid|…
_UV_SLICE = re.compile(
    r"^gnl\|uv\|(?P<acc>[A-Z]{1,2}_?[0-9]{5,}(?:\.[0-9]+)?)"
    r"(?::(?P<start>\d+)-(?P<end>\d+)(?:-\d+)?)?",
    re.IGNORECASE,
)
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# Stay under ~3 req/s without an API key; batching keeps total requests tiny.
_MIN_REQUEST_INTERVAL_S = 0.4
_MAX_IDS_PER_REQUEST = 50
_USER_AGENT = "PlasmidScreen/1.0 (genbank annotation cache; polite eutils client)"


@dataclass(frozen=True)
class ResolvedAccession:
    """NCBI accession derived from a reference hit ID, with optional local slice."""

    hit_id: str
    ncbi_id: str
    slice_start: int | None = None  # 1-based inclusive on the NCBI record
    slice_end: int | None = None  # 1-based inclusive on the NCBI record


def default_genbank_cache_dir() -> Path:
    """``~/.local/share/PlasmidScreen/genbank_cache`` (or platform equivalent)."""
    return Path(get_default_db_path("PlasmidScreen")) / "genbank_cache"


def resolve_ncbi_accession(hit_id: str) -> ResolvedAccession | None:
    """
    Extract a downloadable NCBI accession from a BLAST / Kraken-style hit ID.

    Returns ``None`` for IMGPR and other non-NCBI local names.
    """
    raw = (hit_id or "").strip().split()[0]
    if not raw:
        return None

    m = _UV_SLICE.match(raw)
    if m:
        start = int(m.group("start")) if m.group("start") else None
        end = int(m.group("end")) if m.group("end") else None
        return ResolvedAccession(
            hit_id=hit_id,
            ncbi_id=m.group("acc"),
            slice_start=start,
            slice_end=end,
        )

    head = raw.split("|", 1)[0]
    m = _NCBI_ACC.match(head)
    if m:
        return ResolvedAccession(hit_id=hit_id, ncbi_id=m.group("acc"))

    return None


def _cache_paths(cache_dir: Path, ncbi_id: str) -> tuple[Path, Path]:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", ncbi_id)
    return cache_dir / f"{safe}.gb", cache_dir / f"{safe}.missing"


def _read_cached_gb(path: Path) -> str | None:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    return path.read_text(encoding="utf-8")


def _mark_missing(cache_dir: Path, ncbi_id: str) -> None:
    gb, missing = _cache_paths(cache_dir, ncbi_id)
    missing.write_text(f"not found or empty for {ncbi_id}\n", encoding="utf-8")
    if gb.exists():
        try:
            gb.unlink()
        except OSError:
            pass


def _is_known_missing(cache_dir: Path, ncbi_id: str) -> bool:
    _, missing = _cache_paths(cache_dir, ncbi_id)
    return missing.is_file()


def _efetch_batch(ncbi_ids: Sequence[str], *, api_key: str | None = None) -> str:
    params: dict[str, str] = {
        "db": "nuccore",
        "id": ",".join(ncbi_ids),
        "rettype": "gb",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key
    url = f"{_EFETCH}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _split_gb_records(text: str) -> dict[str, str]:
    """Map VERSION / ACCESSION id -> full GenBank record text."""
    out: dict[str, str] = {}
    if not text.strip():
        return out
    chunks = re.split(r"(?=\nLOCUS\s)", "\n" + text.lstrip())
    for chunk in chunks:
        body = chunk.lstrip("\n")
        if not body.startswith("LOCUS"):
            continue
        version = None
        accession = None
        for line in body.splitlines()[:40]:
            if line.startswith("VERSION"):
                parts = line.split()
                if len(parts) >= 2:
                    version = parts[1]
            elif line.startswith("ACCESSION"):
                parts = line.split()
                if len(parts) >= 2:
                    accession = parts[1]
        key = version or accession
        if key:
            out[key] = body if body.endswith("\n") else body + "\n"
            # Also index bare accession without version.
            bare = key.split(".", 1)[0]
            out.setdefault(bare, out[key])
            if accession:
                out.setdefault(accession, out[key])
    return out


def fetch_genbank_texts(
        ncbi_ids: Iterable[str],
        *,
        cache_dir: str | Path | None = None,
        api_key: str | None = None,
        min_interval_s: float = _MIN_REQUEST_INTERVAL_S,
) -> dict[str, str]:
    """
    Return ``{ncbi_id: genbank_text}`` for IDs that can be fetched or are cached.

    Missing / failed IDs are recorded under ``*.missing`` and omitted from the
    return map so later runs do not re-hammer NCBI.
    """
    cache = Path(cache_dir) if cache_dir else default_genbank_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)

    wanted: list[str] = []
    seen: set[str] = set()
    for raw in ncbi_ids:
        nid = (raw or "").strip()
        if not nid or nid in seen:
            continue
        seen.add(nid)
        wanted.append(nid)

    found: dict[str, str] = {}
    to_fetch: list[str] = []
    for nid in wanted:
        gb_path, _ = _cache_paths(cache, nid)
        if _is_known_missing(cache, nid):
            continue
        cached = _read_cached_gb(gb_path)
        if cached is not None:
            found[nid] = cached
            continue
        to_fetch.append(nid)

    last_req = 0.0
    for i in range(0, len(to_fetch), _MAX_IDS_PER_REQUEST):
        batch = to_fetch[i : i + _MAX_IDS_PER_REQUEST]
        wait = min_interval_s - (time.monotonic() - last_req)
        if wait > 0:
            time.sleep(wait)
        try:
            logger.info(
                "Fetching %d GenBank record(s) from NCBI eutils (cache=%s)",
                len(batch),
                cache,
            )
            text = _efetch_batch(batch, api_key=api_key)
            last_req = time.monotonic()
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            logger.warning("NCBI efetch failed for batch %s: %s", batch, exc)
            continue

        by_id = _split_gb_records(text)
        for nid in batch:
            record = by_id.get(nid) or by_id.get(nid.split(".", 1)[0])
            gb_path, _ = _cache_paths(cache, nid)
            if record and "ORIGIN" in record:
                gb_path.write_text(record, encoding="utf-8")
                found[nid] = record
            else:
                logger.info("No GenBank record for %s; caching as missing", nid)
                _mark_missing(cache, nid)

    return found


def _annotations_from_gb_text(text: str):
    """Parse GenBank text into accession -> PlasmidAnnotation lists."""
    from collections import defaultdict

    from plasmidScreen.src.visualization_code import annotations_from_features

    by_key: dict[str, list] = defaultdict(list)
    for record in SeqIO.parse(StringIO(text), "genbank"):
        annotations = annotations_from_features(record.features)
        if not annotations:
            continue
        keys = {record.id, record.name, record.id.split(".", 1)[0]}
        for xref in record.annotations.get("accessions", []) or []:
            keys.add(str(xref))
            keys.add(str(xref).split(".", 1)[0])
        for key in keys:
            if key:
                by_key[key].extend(annotations)
    return dict(by_key)


def _slice_annotations(
        annotations: Sequence,
        *,
        slice_start: int | None,
        slice_end: int | None,
):
    """Project NCBI feature coords onto a UniVec-style local slice (0-based)."""
    from plasmidScreen.src.visualization_code import PlasmidAnnotation

    if slice_start is None or slice_end is None:
        return list(annotations)
    # SeqIO locations are 0-based half-open relative to the full record.
    origin = slice_start - 1
    local_len = max(0, slice_end - slice_start + 1)
    out = []
    for ann in annotations:
        start = max(ann.start, origin) - origin
        end = min(ann.end, origin + local_len) - origin
        if end <= start:
            continue
        out.append(
            PlasmidAnnotation(
                start=start,
                end=end,
                label=ann.label,
                feature_type=ann.feature_type,
                strand=ann.strand,
                product=getattr(ann, "product", None),
            )
        )
    return out


def collect_ncbi_ids_from_fasta(fasta_path: str | Path) -> list[str]:
    """Unique NCBI accessions resolvable from FASTA headers (sorted)."""
    path = Path(fasta_path)
    found: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            hid = line[1:].strip().split()[0]
            resolved = resolve_ncbi_accession(hid)
            if resolved is not None:
                found.add(resolved.ncbi_id)
    return sorted(found)


@dataclass(frozen=True)
class PrefetchResult:
    """Outcome of a bulk GenBank cache fill."""

    cache_dir: Path
    ids_requested: int
    ids_cached: int
    ids_missing: int
    ids_fetched_this_run: int


def prefetch_genbank_from_fasta(
        fasta_path: str | Path,
        *,
        cache_dir: str | Path | None = None,
        api_key: str | None = None,
        min_interval_s: float = _MIN_REQUEST_INTERVAL_S,
) -> PrefetchResult:
    """
    Download GenBank for every NCBI-resolvable ID in ``fasta_path`` into the cache.

    Already-cached and previously-missing IDs are skipped. Safe to re-run.
    """
    cache = Path(cache_dir) if cache_dir else default_genbank_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    ids = collect_ncbi_ids_from_fasta(fasta_path)
    logger.info(
        "Prefetching GenBank for %d unique NCBI accession(s) from %s → %s",
        len(ids),
        fasta_path,
        cache,
    )

    before_cached = 0
    before_missing = 0
    to_fetch: list[str] = []
    for nid in ids:
        gb_path, missing_path = _cache_paths(cache, nid)
        if _read_cached_gb(gb_path) is not None:
            before_cached += 1
        elif missing_path.is_file():
            before_missing += 1
        else:
            to_fetch.append(nid)

    fetched = fetch_genbank_texts(
        to_fetch,
        cache_dir=cache,
        api_key=api_key,
        min_interval_s=min_interval_s,
    )

    after_cached = 0
    after_missing = 0
    for nid in ids:
        gb_path, missing_path = _cache_paths(cache, nid)
        if _read_cached_gb(gb_path) is not None:
            after_cached += 1
        elif missing_path.is_file():
            after_missing += 1

    return PrefetchResult(
        cache_dir=cache,
        ids_requested=len(ids),
        ids_cached=after_cached,
        ids_missing=after_missing,
        ids_fetched_this_run=len(fetched),
    )


def annotations_for_hit_ids(
        hit_ids: Iterable[str],
        *,
        cache_dir: str | Path | None = None,
        api_key: str | None = None,
        fetch: bool = False,
) -> dict[str, list]:
    """
    Build ``hit_id -> annotations`` for reference accessions that resolve to NCBI.

    Annotations are keyed by the **original hit ID** (BLAST subject / accession) so
    pileup lookup can use the DB accession verbatim. By default only the on-disk
    cache is used (airgapped-safe after ``prefetch_genbank_from_fasta``). Set
    ``fetch=True`` to download missing records from NCBI during this call.
    """
    resolved: list[ResolvedAccession] = []
    for hid in hit_ids:
        r = resolve_ncbi_accession(hid)
        if r is not None:
            resolved.append(r)

    if not resolved:
        logger.info("No NCBI-resolvable reference IDs among pileup accessions.")
        return {}

    by_ncbi = {r.ncbi_id for r in resolved}
    texts: dict[str, str] = {}
    if fetch:
        texts = fetch_genbank_texts(
            by_ncbi, cache_dir=cache_dir, api_key=api_key
        )
    else:
        cache = Path(cache_dir) if cache_dir else default_genbank_cache_dir()
        for nid in by_ncbi:
            gb_path, _ = _cache_paths(cache, nid)
            cached = _read_cached_gb(gb_path)
            if cached:
                texts[nid] = cached

    parsed_by_ncbi: dict[str, list] = {}
    for nid, text in texts.items():
        ann_map = _annotations_from_gb_text(text)
        anns = ann_map.get(nid) or ann_map.get(nid.split(".", 1)[0])
        if anns is None and ann_map:
            anns = next(iter(ann_map.values()))
        parsed_by_ncbi[nid] = list(anns or [])

    out: dict[str, list] = {}
    for r in resolved:
        anns = parsed_by_ncbi.get(r.ncbi_id) or parsed_by_ncbi.get(
            r.ncbi_id.split(".", 1)[0]
        )
        if not anns:
            continue
        out[r.hit_id] = _slice_annotations(
            anns, slice_start=r.slice_start, slice_end=r.slice_end
        )
        out.setdefault(r.ncbi_id, list(anns))
        out.setdefault(r.ncbi_id.split(".", 1)[0], list(anns))

    logger.info(
        "Prepared GenBank annotations for %d/%d resolvable hit accession(s)",
        sum(1 for r in resolved if r.hit_id in out),
        len(resolved),
    )
    return out
