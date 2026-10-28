"""Perseus (Kraken2 post-processor) integration for engineered taxid detection."""
from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

from plasmidScreen.lib.exceptions import InputError, PerseusError

logger = logging.getLogger(__name__)

# Bump when Perseus wiring changes; must appear in screen logs if the new build is loaded.
PERSEUS_BRIDGE_BUILD = "2026-10-27c"

# NCBI "synthetic construct" — same taxid used by the engineered k-mer scan.
ENGINEERED_TAXID = 32630

PERSEUS_TAXID_COLUMN = "perseus_taxid"
SEQUENCE_ID_COLUMN = "sequence_id"

# Kraken2 --use-names puts e.g. "Escherichia coli (taxid 562)" in column 3.
_TAXID_IN_NAME_RE = re.compile(r"\(taxid\s+(\d+)\)", re.IGNORECASE)


def kraken_read_id(line: str) -> str | None:
    """Return the sequence ID (column 2) from a Kraken2 classification line."""
    parts = line.rstrip("\n\r").split("\t")
    if len(parts) < 2:
        return None
    return parts[1]


def kraken_line_status(line: str) -> str | None:
    """Return the Kraken2 classification status (column 1: ``C`` or ``U``)."""
    parts = line.rstrip("\n\r").split("\t", 1)
    if not parts or not parts[0]:
        return None
    return parts[0]


def count_kraken_classified(lines: Iterable[str]) -> tuple[int, int]:
    """Return ``(n_classified, n_unclassified)`` from Kraken2 lines."""
    n_c = 0
    n_u = 0
    for line in lines:
        status = kraken_line_status(line)
        if status == "C":
            n_c += 1
        elif status == "U":
            n_u += 1
    return n_c, n_u


def filter_kraken_lines(
        lines: Iterable[str],
        read_ids: set[str],
) -> list[str]:
    """Keep Kraken classification lines whose read ID is in ``read_ids``."""
    if not read_ids:
        return []
    out: list[str] = []
    for line in lines:
        rid = kraken_read_id(line)
        if rid is not None and rid in read_ids:
            out.append(line if line.endswith("\n") else line + "\n")
    return out


def sanitize_read_id_for_perseus(read_id: str) -> str:
    """
    Strip characters Perseus mis-parses in sequence IDs.

    Perseus treats ``|`` as MESS ``id|taxid`` syntax (``row.ID.split("|")[1]``).
    Colons are removed so IDs cannot collide with ``taxid:count`` k-mer tokens.
    """
    return read_id.replace("|", "_").replace(":", "_")


def extract_taxid_from_kraken_taxonomy(taxonomy_field: str) -> str:
    """
    Return a bare taxid string from a Kraken column-3 value.

    Handles ``--use-names`` fields like ``Escherichia coli (taxid 562)`` and
    plain numeric taxids. Returns the original field if no taxid is found.
    """
    text = taxonomy_field.strip()
    if text.isdigit():
        return text
    match = _TAXID_IN_NAME_RE.search(text)
    if match:
        return match.group(1)
    return text


def rewrite_kraken_ids_for_perseus(
        lines: Iterable[str],
) -> tuple[list[str], dict[str, str]]:
    """
    Rewrite Kraken lines for Perseus compatibility.

    - Column 2 (read ID): replace ``|`` and ``:`` (Perseus MESS / token pitfalls)
    - Column 3 (taxonomy): coerce ``--use-names`` text to a bare taxid when possible

    Returns ``(rewritten_lines, sanitized_id -> original_id)``. Collision-safe:
    if two originals sanitize to the same string, a ``__N`` suffix is added.
    """
    rewritten: list[str] = []
    sanitized_to_original: dict[str, str] = {}
    used: set[str] = set()
    n_id_changed = 0
    n_tax_simplified = 0

    for line in lines:
        raw = line if line.endswith("\n") else line + "\n"
        parts = raw.rstrip("\n\r").split("\t")
        if len(parts) < 2:
            rewritten.append(raw)
            continue
        original = parts[1]
        base = sanitize_read_id_for_perseus(original)
        sanitized = base
        if sanitized in used and sanitized_to_original.get(sanitized) != original:
            n = 2
            while f"{base}__{n}" in used:
                n += 1
            sanitized = f"{base}__{n}"
        if sanitized != original:
            n_id_changed += 1
        used.add(sanitized)
        sanitized_to_original[sanitized] = original
        parts[1] = sanitized

        if len(parts) >= 3:
            simplified = extract_taxid_from_kraken_taxonomy(parts[2])
            if simplified != parts[2]:
                n_tax_simplified += 1
                parts[2] = simplified

        rewritten.append("\t".join(parts) + "\n")

    if n_id_changed:
        logger.info(
            "Sanitized %d Kraken read ID(s) for Perseus (removed '|' and/or ':').",
            n_id_changed,
        )
    if n_tax_simplified:
        logger.info(
            "Simplified %d Kraken taxonomy field(s) to bare taxid for Perseus "
            "(from --use-names style).",
            n_tax_simplified,
        )
    return rewritten, sanitized_to_original


def remap_sanitized_taxids(
        taxids_by_sanitized: dict[str, int],
        sanitized_to_original: dict[str, str],
) -> dict[str, int]:
    """Map Perseus sequence_id keys back to original read IDs."""
    if not sanitized_to_original:
        return taxids_by_sanitized
    return {
        sanitized_to_original.get(sid, sid): taxid
        for sid, taxid in taxids_by_sanitized.items()
    }


def _assert_extract_shards(shards_dir: Path) -> None:
    """
    Perseus filter crashes with ``manifest has no outputs`` when extract wrote
    zero shards. Fail earlier with a clearer message.

    Also recovers when ``.pt`` shard files exist on disk but were never recorded
    in the Manager.list-backed manifest (empty ``outputs`` after a long extract).
    """
    manifest = shards_dir / "manifest.json"
    if not manifest.is_file():
        raise PerseusError(
            1,
            f"Perseus extract did not write shard manifest: {manifest}",
        )
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PerseusError(
            1,
            f"Could not read Perseus shard manifest {manifest}: {exc}",
        ) from exc
    outputs = data.get("outputs") or []
    if outputs:
        return

    pt_files = sorted(p.name for p in shards_dir.glob("*.pt"))
    if pt_files:
        data["outputs"] = pt_files
        try:
            manifest.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        except OSError as exc:
            raise PerseusError(
                1,
                f"Perseus manifest empty but found {len(pt_files)} .pt shard(s); "
                f"failed to rewrite {manifest}: {exc}",
            ) from exc
        logger.warning(
            "Perseus manifest had empty outputs; recovered %d .pt shard(s) on disk "
            "into %s (Manager.list registration miss).",
            len(pt_files),
            manifest,
        )
        return

    raise PerseusError(
        1,
        "Perseus extract produced no feature shards (empty manifest, no .pt files). "
        "Perseus only builds features for Kraken2 lines with status 'C' "
        "(classified), and skips reads whose taxids cannot be resolved in "
        "the ETE3 DB from `perseus setup`. Check that Natural reads include "
        "classified lines with a k-mer LCA string, and that --perseus-db is "
        f"a complete taxonomy DB. Manifest: {manifest}",
    )


def parse_perseus_output(
        path: str | Path | None = None,
        *,
        lines: Iterable[str] | None = None,
        engineered_taxid: int = ENGINEERED_TAXID,
) -> dict[str, int]:
    """
    Parse Perseus filter TSV into ``{sequence_id: perseus_taxid}``.

    Accepts a file path and/or in-memory lines. The header must include
    ``sequence_id`` and ``perseus_taxid`` (Perseus v1.2.0 layout).
    """
    if lines is None:
        if path is None:
            raise ValueError("Provide path or lines to parse Perseus output.")
        try:
            text = Path(path).read_text(encoding="utf-8")
        except OSError as exc:
            raise InputError(f"Could not read Perseus output {path}: {exc}") from exc
        lines = text.splitlines()

    row_iter = iter(lines)
    header_line: str | None = None
    for raw in row_iter:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        header_line = line
        break
    if header_line is None:
        return {}

    headers = header_line.split("\t")
    try:
        id_idx = headers.index(SEQUENCE_ID_COLUMN)
        tax_idx = headers.index(PERSEUS_TAXID_COLUMN)
    except ValueError as exc:
        raise InputError(
            "Perseus output missing required columns "
            f"'{SEQUENCE_ID_COLUMN}' and/or '{PERSEUS_TAXID_COLUMN}'. "
            f"Found: {headers}"
        ) from exc

    result: dict[str, int] = {}
    for raw in row_iter:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) <= max(id_idx, tax_idx):
            continue
        read_id = cols[id_idx].strip()
        tax_raw = cols[tax_idx].strip()
        if not read_id or not tax_raw or tax_raw.upper() == "NA":
            continue
        try:
            result[read_id] = int(float(tax_raw))
        except ValueError:
            logger.warning(
                "Skipping Perseus row for %s with non-integer taxid %r",
                read_id,
                tax_raw,
            )
    _ = engineered_taxid  # documented for callers; filtering is separate
    return result


def engineered_read_ids_from_perseus(
        taxids_by_read: dict[str, int],
        *,
        engineered_taxid: int = ENGINEERED_TAXID,
) -> set[str]:
    """Return read IDs whose Perseus taxid equals the engineered synthetic taxid."""
    return {
        rid for rid, taxid in taxids_by_read.items() if taxid == engineered_taxid
    }


def _log_perseus_stream(step: str, stream_name: str, text: str | None) -> None:
    if not text or not text.strip():
        return
    # Perseus logs progress to stderr; keep the last chunk readable.
    lines = text.strip().splitlines()
    tail = lines[-40:] if len(lines) > 40 else lines
    logger.info("Perseus %s %s (last %d line(s)):", step, stream_name, len(tail))
    for line in tail:
        logger.info("  %s", line)


def _run_perseus_cmd(cmd: list[str], *, step: str) -> None:
    logger.info("Running Perseus %s: %s", step, " ".join(cmd))
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise PerseusError(
            127,
            "perseus executable not found; is Perseus installed and on PATH?",
        ) from exc
    except subprocess.CalledProcessError as exc:
        _log_perseus_stream(step, "stderr", exc.stderr)
        _log_perseus_stream(step, "stdout", exc.stdout)
        raise PerseusError(exc.returncode, exc.stderr) from exc
    _log_perseus_stream(step, "stderr", proc.stderr)
    _log_perseus_stream(step, "stdout", proc.stdout)


def run_perseus(
        kraken_file: str | Path,
        perseus_db: str | Path,
        output_path: str | Path,
        *,
        work_dir: str | Path | None = None,
        extract_threads: int = 0,
) -> Path:
    """
    Run ``perseus extract`` then ``perseus filter`` on a Kraken2 classifications file.

    ``perseus_db`` is the ETE3 taxonomy directory created by ``perseus setup``.
    ``extract_threads`` is passed to ``perseus extract --threads`` (``0`` = Perseus
    auto / all available CPUs, same as Perseus CLI default).

    Returns the path to the filtered Perseus TSV.
    """
    kraken_file = Path(kraken_file)
    perseus_db = Path(perseus_db)
    output_path = Path(output_path)

    if not kraken_file.is_file():
        raise InputError(f"Kraken classifications for Perseus not found: {kraken_file}")
    if kraken_file.stat().st_size == 0:
        raise InputError(f"Kraken classifications for Perseus is empty: {kraken_file}")
    if not perseus_db.is_dir():
        raise InputError(
            f"Perseus taxonomy database directory not found: {perseus_db}. "
            "Run `perseus setup <db_path>` first."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    threads = max(0, int(extract_threads))

    def _extract_and_filter(root: Path) -> Path:
        shards_dir = root / "shards"
        shards_dir.mkdir(parents=True, exist_ok=True)
        _run_perseus_cmd(
            [
                "perseus",
                "extract",
                str(kraken_file),
                str(shards_dir),
                str(perseus_db),
                "--threads",
                str(threads),
            ],
            step="extract",
        )
        _assert_extract_shards(shards_dir)
        _run_perseus_cmd(
            [
                "perseus",
                "filter",
                str(shards_dir),
                str(kraken_file),
                str(output_path),
                str(perseus_db),
            ],
            step="filter",
        )
        if not output_path.is_file():
            raise PerseusError(
                1,
                f"Perseus filter did not produce output at {output_path}",
            )
        return output_path

    if work_dir is not None:
        root = Path(work_dir)
        root.mkdir(parents=True, exist_ok=True)
        logger.info("Perseus work directory: %s", root)
        return _extract_and_filter(root)

    # Durable temp: keep on failure so empty-manifest / extract issues can be inspected.
    root = Path(tempfile.mkdtemp(prefix="plasmidscreen_perseus_"))
    logger.info("Perseus work directory: %s", root)
    try:
        result = _extract_and_filter(root)
    except Exception:
        logger.error(
            "Perseus failed; work directory preserved for debugging: %s "
            "(contains natural_kraken.out / shards/manifest.json when present)",
            root,
        )
        raise
    shutil.rmtree(root, ignore_errors=True)
    return result


def resolve_perseus_taxids(
        *,
        kraken_lines: list[str],
        natural_read_ids: set[str],
        perseus_db: str | Path | None,
        run_perseus_enabled: bool = True,
        perseus_output_path: str | Path | None = None,
        debug_write_perseus_output: bool = False,
        work_dir: str | Path | None = None,
        extract_threads: int = 0,
) -> tuple[dict[str, int], Path | None]:
    """
    Obtain Perseus taxids for Natural reads.

    When ``run_perseus_enabled`` is True, writes a Natural-only Kraken subset and
    runs Perseus. When False, loads a precomputed Perseus TSV from
    ``perseus_output_path``.

    ``extract_threads`` is forwarded to ``perseus extract --threads`` (0 = auto).

    Returns ``(taxids_by_read, output_path_or_None)``.
    """
    if not natural_read_ids:
        logger.info("No Natural reads; skipping Perseus.")
        return {}, None

    if not run_perseus_enabled:
        if perseus_output_path is None:
            raise ValueError(
                "perseus_output_path is required when run_perseus=False."
            )
        path = Path(perseus_output_path)
        if not path.is_file():
            raise InputError(f"Perseus output not found: {path}")
        logger.info("Loading Perseus output from %s", path)
        return parse_perseus_output(path), path

    if perseus_db is None:
        raise ValueError("perseus_db is required when run_perseus=True.")
    if debug_write_perseus_output and perseus_output_path is None:
        raise ValueError(
            "perseus_output_path is required when debug_write_perseus_output=True."
        )

    subset = filter_kraken_lines(kraken_lines, natural_read_ids)
    if not subset:
        logger.info(
            "No Kraken lines matched %d Natural read ID(s); skipping Perseus.",
            len(natural_read_ids),
        )
        return {}, None

    n_classified, n_unclassified = count_kraken_classified(subset)
    logger.info(
        "Perseus prep build=%s module=%s: %d Natural line(s), %d classified (C), "
        "%d unclassified (U), extract_threads=%s",
        PERSEUS_BRIDGE_BUILD,
        Path(__file__).resolve(),
        len(subset),
        n_classified,
        n_unclassified,
        extract_threads,
    )
    if n_classified == 0:
        logger.warning(
            "Skipping Perseus: none of %d Natural Kraken line(s) are classified "
            "(status 'C'). Perseus extract only processes classified reads and "
            "would produce an empty shard manifest.",
            len(subset),
        )
        return {}, None

    subset, id_map = rewrite_kraken_ids_for_perseus(subset)

    def _run_in_work(root: Path, out_path: Path) -> Path:
        root.mkdir(parents=True, exist_ok=True)
        kraken_path = root / "natural_kraken.out"
        kraken_path.write_text("".join(subset), encoding="utf-8")
        logger.info("Wrote Perseus Kraken subset (%d line(s)) → %s", len(subset), kraken_path)
        try:
            return run_perseus(
                kraken_path,
                perseus_db,
                out_path,
                work_dir=root,
                extract_threads=extract_threads,
            )
        except Exception:
            logger.error(
                "Perseus failed; work directory preserved for debugging: %s",
                root,
            )
            raise

    if debug_write_perseus_output:
        out = Path(perseus_output_path)  # type: ignore[arg-type]
        work = Path(work_dir) if work_dir is not None else out.parent
        saved = _run_in_work(work, out)
        taxids = remap_sanitized_taxids(parse_perseus_output(saved), id_map)
        return taxids, saved

    if work_dir is not None:
        root = Path(work_dir)
        out = Path(perseus_output_path) if perseus_output_path else root / "perseus.tsv"
        saved = _run_in_work(root, out)
        taxids = remap_sanitized_taxids(parse_perseus_output(saved), id_map)
        return taxids, Path(perseus_output_path) if perseus_output_path else saved

    root = Path(tempfile.mkdtemp(prefix="plasmidscreen_perseus_"))
    out = Path(perseus_output_path) if perseus_output_path else root / "perseus.tsv"
    try:
        saved = _run_in_work(root, out)
        taxids = remap_sanitized_taxids(parse_perseus_output(saved), id_map)
    except Exception:
        raise
    else:
        shutil.rmtree(root, ignore_errors=True)
        return taxids, Path(perseus_output_path) if perseus_output_path else None
