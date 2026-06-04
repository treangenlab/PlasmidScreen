"""Engineered k-mer screening workflow and Kraken integration."""
from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

try:
    from numba import jit
except ModuleNotFoundError:  # pragma: no cover
    def jit(*_args: object, **_kwargs: object):  # type: ignore[misc]
        def _wrap(fn: object) -> object:
            return fn

        return _wrap

from plasmidScreen.src.codon_usage.codon_usage_db import default_codon_usage_dir
from plasmidScreen.lib.models import (
    CodonAdaptationResult,
    EngineeredScanResult,
    ReadFlagDetail,
    ReadEngineeringLabel,
    ScreenResult,
    compute_engineered_overall,
)
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
    parse_kraken_lines,
    write_codon_adaptation_results_tsv,
)


@jit(nopython=True)
def fast_window_logic(
    tids,
    counts,
    window_size,
    threshold,
):
    target_tid = 32630
    max_kmers = window_size - 21 + 1

    eng_count = 0
    non_eng_count = 0
    total_count = 0
    max_eng_seen = 0

    for i in range(len(tids)):
        tid = tids[i]
        count = counts[i]

        for _ in range(count):
            if total_count < max_kmers:
                if tid == target_tid:
                    eng_count += 1
                else:
                    non_eng_count += 1
                total_count += 1
            else:
                if tid == target_tid:
                    if eng_count + 1 < max_kmers:
                        eng_count += 1
                    else:
                        eng_count = max_kmers

                    if non_eng_count - 1 > 0:
                        non_eng_count -= 1
                    else:
                        non_eng_count = 0
                else:
                    if non_eng_count + 1 < max_kmers:
                        non_eng_count += 1
                    else:
                        non_eng_count = max_kmers

                    if eng_count - 1 > 0:
                        eng_count -= 1
                    else:
                        eng_count = 0

            if eng_count > max_eng_seen:
                max_eng_seen = eng_count

            if eng_count >= threshold:
                return True, max_eng_seen

    return False, max_eng_seen


def engineered_scan_to_tsv_lines(
    labels: list[ReadEngineeringLabel],
    *,
    threshold: int,
    window_size: int,
    kmer_max_by_read: dict[str, int] | None = None,
) -> list[str]:
    """Format engineered k-mer scan labels as TSV lines (header + rows)."""
    header = (
        "Label\tRead_ID\tMethods\tEngineeredKmerMaxInWindow\t"
        "KmerThreshold\tWindowSize"
    )
    lines = [header]
    for lbl in labels:
        methods = "engineered_kmer_scan" if lbl.label == "Synthetic" else ""
        max_eng = (kmer_max_by_read or {}).get(lbl.read_id, 0)
        lines.append(
            f"{lbl.label}\t{lbl.read_id}\t{methods}\t{max_eng}\t"
            f"{threshold}\t{window_size}"
        )
    return lines


def write_engineered_report_tsv(
    output_path: str | Path,
    labels: list[ReadEngineeringLabel],
    *,
    threshold: int,
    window_size: int,
    kmer_max_by_read: dict[str, int] | None = None,
) -> str:
    """Write engineered k-mer scan report TSV; returns output path."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = engineered_scan_to_tsv_lines(
        labels,
        threshold=threshold,
        window_size=window_size,
        kmer_max_by_read=kmer_max_by_read,
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(out)


class Workflow:
    """
    End-to-end screening: Kraken2 engineered k-mer scan and optional DIAMOND codon CAI.

    Kraken output is held in memory by default. Codon analysis uses DIAMOND blastx for
    ORF intervals and host taxids, restricted to reads labeled Natural by the k-mer scan.
    """

    fasta_file: Path
    report_output_path: Path | None
    write_engineered_report: bool
    kraken_db: Path
    kraken_output_path: Path | None
    debug_write_kraken_output: bool
    debug_write_kraken_report: bool
    _kraken_lines: list[str] | None
    _kraken_data: dict[str, tuple[str, str, int, str]] | None
    threshold: int
    window_size: int
    max_threads: int
    codon_usage_output_path: Path | None
    codon_usage_dir: Path
    run_kraken_enabled: bool
    run_codon_usage: bool
    codon_cai_engineered_threshold: float | None
    diamond_db: Path | None
    diamond_threads: int
    diamond_extra_args: list[str] | None
    diamond_output_path: Path | None
    debug_write_diamond_output: bool
    run_diamond_enabled: bool
    _diamond_output_saved: Path | None

    def __init__(
        self,
        fasta_file: str | Path,
        output_report_path: str | Path | None,
        kraken_db: str | Path,
        threads: int,
        kraken_raw_output: str | Path | None,
        *,
        write_engineered_report: bool = False,
        window_size: int = 200,
        engineered_kmer_threshold: int = 25,
        codon_usage_output_path: str | Path | None = None,
        codon_usage_dir: str | Path | None = None,
        run_kraken: bool = True,
        debug_write_kraken_output: bool = False,
        debug_write_kraken_report: bool = False,
        run_codon_usage: bool = True,
        codon_cai_engineered_threshold: float | None = None,
        diamond_db: str | Path | None = None,
        diamond_threads: int = 4,
        diamond_extra_args: list[str] | None = None,
        diamond_output_path: str | Path | None = None,
        debug_write_diamond_output: bool = False,
        run_diamond: bool = True,
    ) -> None:
        self.fasta_file = Path(fasta_file)
        self.report_output_path = (
            Path(output_report_path) if output_report_path else None
        )
        self.write_engineered_report = write_engineered_report
        self.kraken_db = Path(kraken_db)
        self.kraken_output_path = Path(kraken_raw_output) if kraken_raw_output else None
        self.threshold = engineered_kmer_threshold
        self.window_size = window_size
        self.max_threads = threads
        self.codon_usage_output_path = (
            Path(codon_usage_output_path) if codon_usage_output_path else None
        )
        self.codon_usage_dir = (
            Path(codon_usage_dir) if codon_usage_dir else default_codon_usage_dir()
        )
        self.run_kraken_enabled = run_kraken
        self.debug_write_kraken_output = debug_write_kraken_output
        self.debug_write_kraken_report = debug_write_kraken_report
        self.run_codon_usage = run_codon_usage
        self.codon_cai_engineered_threshold = codon_cai_engineered_threshold
        self._kraken_lines = None
        self._kraken_data = None
        self.diamond_db = Path(diamond_db) if diamond_db else None
        self.diamond_threads = diamond_threads
        self.diamond_extra_args = diamond_extra_args
        self.diamond_output_path = (
            Path(diamond_output_path) if diamond_output_path else None
        )
        self.debug_write_diamond_output = debug_write_diamond_output
        self.run_diamond_enabled = run_diamond
        self._diamond_output_saved = None

    def _ensure_kraken_in_memory(self) -> None:
        if self._kraken_lines is not None and self._kraken_data is not None:
            return
        if self.kraken_output_path and self.kraken_output_path.exists():
            lines = self.kraken_output_path.read_text(encoding="utf-8").splitlines(True)
            self._kraken_lines = lines
            self._kraken_data = parse_kraken_lines(lines)

    def run_kraken(self) -> None:
        try:
            logging.info("Running Kraken on %s", self.fasta_file)
            if self.debug_write_kraken_output and self.kraken_output_path is None:
                raise ValueError(
                    "debug_write_kraken_output=True requires kraken_output_path to be set."
                )
            cmd = [
                "kraken2",
                "--db",
                str(self.kraken_db),
                "--report-minimizer-data",
                "--use-names",
                str(self.fasta_file),
                "--threads",
                str(self.max_threads),
            ]

            # Raw Kraken output: in-memory by default (Kraken2 writes classifications to stdout
            # when --output is omitted); write to kraken_output_path only in debug mode.
            if self.debug_write_kraken_output:
                cmd += ["--output", str(self.kraken_output_path)]

            # Kraken report file: debug-only.
            if self.debug_write_kraken_report:
                report_file_path = self.kraken_db.parent.joinpath("report.txt")
                cmd += ["--report", str(report_file_path)]

            logging.info("Running the command: %s", " ".join(cmd))

            if self.debug_write_kraken_output:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                self._ensure_kraken_in_memory()
            else:
                proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
                lines = proc.stdout.splitlines(True)
                self._kraken_lines = lines
                self._kraken_data = parse_kraken_lines(lines)
        except subprocess.CalledProcessError as e:
            logging.error(
                "Kraken2 failed (Code %s). Error details: %s", e.returncode, e.stderr
            )
            sys.exit(1)
        logging.info("All files processed.")

    @staticmethod
    def parse_and_run(
        kmer_pos_info: str, window_size: int, threshold: int
    ) -> tuple[bool, int]:
        raw_data = kmer_pos_info.replace("|:|", "").split()
        if not raw_data:
            return False, 0

        window_size = int(window_size)
        threshold = int(threshold)
        max_kmers = max(window_size - 21 + 1, 1)

        n = len(raw_data)
        tids = np.zeros(n, dtype=np.int64)
        counts = np.zeros(n, dtype=np.int64)

        try:
            for i, item in enumerate(raw_data):
                t, c = item.split(":", 1)
                if t == "A":
                    tids[i] = -1
                else:
                    tids[i] = int(t)
                counts[i] = min(int(c), max_kmers)
        except (ValueError, IndexError, OverflowError):
            return False, 0

        tids = np.ascontiguousarray(tids, dtype=np.int64)
        counts = np.ascontiguousarray(counts, dtype=np.int64)
        try:
            return fast_window_logic(tids, counts, window_size, threshold)
        except TypeError:
            # Older Numba builds can fail to unbox int32-annotated arrays; int64 + fallback.
            return Workflow._fast_window_logic_python(
                tids, counts, window_size, threshold
            )

    @staticmethod
    def _fast_window_logic_python(
        tids: NDArray[np.int64],
        counts: NDArray[np.int64],
        window_size: int,
        threshold: int,
    ) -> tuple[bool, int]:
        """Pure-Python fallback when Numba cannot compile/unbox inputs."""
        target_tid = 32630
        max_kmers = window_size - 21 + 1

        eng_count = 0
        non_eng_count = 0
        total_count = 0
        max_eng_seen = 0

        for i in range(len(tids)):
            tid = int(tids[i])
            count = int(counts[i])

            for _ in range(count):
                if total_count < max_kmers:
                    if tid == target_tid:
                        eng_count += 1
                    else:
                        non_eng_count += 1
                    total_count += 1
                else:
                    if tid == target_tid:
                        if eng_count + 1 < max_kmers:
                            eng_count += 1
                        else:
                            eng_count = max_kmers

                        if non_eng_count - 1 > 0:
                            non_eng_count -= 1
                        else:
                            non_eng_count = 0
                    else:
                        if non_eng_count + 1 < max_kmers:
                            non_eng_count += 1
                        else:
                            non_eng_count = max_kmers

                        if eng_count - 1 > 0:
                            eng_count -= 1
                        else:
                            eng_count = 0

                if eng_count > max_eng_seen:
                    max_eng_seen = eng_count

                if eng_count >= threshold:
                    return True, max_eng_seen

        return False, max_eng_seen

    def scan_engineered_blocks_kraken(self) -> EngineeredScanResult:
        result = EngineeredScanResult()
        self._ensure_kraken_in_memory()
        kmer_max_by_read: dict[str, int] = {}
        entries = self._kraken_lines
        for entry in tqdm(entries, total=len(entries)):
            categories = entry.split("\t")
            if len(categories) < 2:
                continue
            synthetic_boolean, max_eng = self.parse_and_run(
                categories[-1], self.window_size, self.threshold
            )
            read_id = categories[1]
            kmer_max_by_read[read_id] = max_eng
            if synthetic_boolean:
                result.synthetic_count += 1
                label = "Synthetic"
            else:
                result.natural_count += 1
                label = "Natural"
            result.labels.append(ReadEngineeringLabel(read_id=read_id, label=label))

        if self.write_engineered_report:
            if self.report_output_path is None:
                raise ValueError(
                    "report_output_path is required when write_engineered_report=True."
                )
            write_engineered_report_tsv(
                self.report_output_path,
                result.labels,
                threshold=self.threshold,
                window_size=self.window_size,
                kmer_max_by_read=kmer_max_by_read,
            )

        logging.info(
            "Engineered k-mer scan complete: %s/%s synthetic reads.",
            result.synthetic_count,
            result.synthetic_count + result.natural_count,
        )
        return result

    def run_codon_adaptation(self, natural_read_ids: set[str]) -> list[CodonAdaptationResult]:
        if not natural_read_ids:
            logging.info("No reads labeled Natural; skipping codon usage analysis.")
            return []

        if self.run_diamond_enabled and self.diamond_db is None:
            raise ValueError(
                "DIAMOND host-taxonomy is required for codon CAI. "
                "Provide diamond_db to Workflow (or --diamond-db in CLI)."
            )
        if not self.run_diamond_enabled and self.diamond_output_path is None:
            raise ValueError(
                "diamond_output_path is required when run_diamond=False "
                "(precomputed DIAMOND TSV)."
            )
        if self.debug_write_diamond_output and self.diamond_output_path is None:
            raise ValueError(
                "diamond_output_path is required when debug_write_diamond_output=True."
            )
        results, diamond_path = analyze_codon_adaptation(
            self.fasta_file,
            diamond_db=self.diamond_db,
            diamond_threads=self.diamond_threads,
            diamond_extra_args=self.diamond_extra_args,
            run_diamond=self.run_diamond_enabled,
            diamond_output_path=self.diamond_output_path,
            debug_write_diamond_output=self.debug_write_diamond_output,
            codon_usage_dir=self.codon_usage_dir,
            include_read_ids=natural_read_ids,
        )
        self._diamond_output_saved = diamond_path
        return results

    def run(self) -> ScreenResult:
        if self.run_kraken_enabled:
            self.run_kraken()
        engineered_scan = self.scan_engineered_blocks_kraken()

        codon_results: list[CodonAdaptationResult] = []
        codon_path: Path | None = None

        if self.run_codon_usage and engineered_scan.natural_read_ids:
            if engineered_scan.any_synthetic:
                logging.info(
                    "Engineered reads detected; running codon usage only on reads labeled Natural."
                )

            logging.info(
                "Running codon usage analysis on %d Natural reads (reference: %s)",
                len(engineered_scan.natural_read_ids),
                self.codon_usage_dir,
            )
            codon_results = self.run_codon_adaptation(engineered_scan.natural_read_ids)
            if self.codon_usage_output_path is not None:
                codon_path_str = write_codon_adaptation_results_tsv(
                    self.codon_usage_output_path,
                    codon_results,
                    cai_engineered_threshold=self.codon_cai_engineered_threshold,
                )
                codon_path = Path(codon_path_str)

        codon_by_read: dict[str, CodonAdaptationResult] = {
            r.read_id: r for r in codon_results
        }

        kmer_max_by_read: dict[str, int] = {}
        if self._kraken_lines is not None:
            for line in self._kraken_lines:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                rid = parts[1]
                _hit, max_eng = self.parse_and_run(
                    parts[-1], self.window_size, self.threshold
                )
                kmer_max_by_read[rid] = max_eng

        per_read: list[ReadFlagDetail] = []
        for lbl in engineered_scan.labels:
            codon = codon_by_read.get(lbl.read_id)
            cai = codon.cai_vs_host if codon else None
            engineered_by_codon: bool | None = None
            if cai is not None and self.codon_cai_engineered_threshold is not None:
                engineered_by_codon = cai < self.codon_cai_engineered_threshold

            engineered_by_kmer = lbl.label == "Synthetic"
            engineered_overall = compute_engineered_overall(
                engineered_by_kmer_scan=engineered_by_kmer,
                engineered_by_codon_cai=engineered_by_codon,
            )
            overall_label: Literal["Natural", "Synthetic"] = (
                "Synthetic" if engineered_overall else "Natural"
            )

            per_read.append(
                ReadFlagDetail(
                    read_id=lbl.read_id,
                    kmer_label=lbl.label,
                    engineered_by_kmer_scan=engineered_by_kmer,
                    engineered_overall=engineered_overall,
                    overall_label=overall_label,
                    engineered_kmer_max_in_window=kmer_max_by_read.get(lbl.read_id),
                    engineered_kmer_threshold=self.threshold,
                    engineered_kmer_window_size=self.window_size,
                    cai_vs_host=cai,
                    engineered_by_codon_cai=engineered_by_codon,
                    codon_cai_threshold=self.codon_cai_engineered_threshold,
                )
            )

        engineered_path = (
            self.report_output_path
            if self.write_engineered_report and self.report_output_path
            else None
        )
        return ScreenResult(
            engineered_scan=engineered_scan,
            codon_adaptation=codon_results,
            per_read=per_read,
            engineered_report_path=engineered_path,
            codon_usage_report_path=codon_path,
            diamond_output_path=self._diamond_output_saved,
            engineered_kmer_threshold=self.threshold,
            engineered_kmer_window_size=self.window_size,
            codon_cai_engineered_threshold=self.codon_cai_engineered_threshold,
        )
