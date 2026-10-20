"""Engineered k-mer screening workflow and Kraken integration."""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Literal, List

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from plasmidScreen.lib.funcs import call_synthetic_reads
from plasmidScreen.lib.exceptions import InputError, KrakenError

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
    BACKGROUND_RATE,
)
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
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


class Workflow:
    """
    End-to-end screening: Kraken2 engineered k-mer scan and optional DIAMOND codon CAI.

    Kraken output is held in memory by default. Codon analysis uses DIAMOND blastx for
    ORF intervals and host taxids, restricted to reads labeled Natural by the k-mer scan.
    """

    def __init__(
            self,
            fasta_file: str | Path,
            seq_type: int,
            output_report_path: str | Path | None,
            kraken_db: str | Path,
            threads: int,
            kraken_raw_output: str | Path | None,
            *,
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
            diamond_output_path: str | Path | None = None,
            debug_write_diamond_output: bool = False,
            run_diamond: bool = True,
            quiet_mode: bool = False,
            mem_mode: int = 0,
            kmer_size: int = 35,
            minimizer_size: int = 31,
            background_rate: float = BACKGROUND_RATE
    ) -> None:
        self.fasta_file = Path(fasta_file)
        self.seq_type = seq_type
        self.report_output_path = (
            Path(output_report_path) if output_report_path else None
        )
        self.kraken_db = Path(kraken_db)
        self.kraken_output_path = Path(kraken_raw_output) if kraken_raw_output else None
        self.threshold = engineered_kmer_threshold
        self.window_size = window_size
        self.max_threads = threads
        self.k = kmer_size
        self.l = minimizer_size
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
        self.diamond_db = Path(diamond_db) if diamond_db else None
        self.diamond_threads = self.max_threads
        self.diamond_output_path = (
            Path(diamond_output_path) if diamond_output_path else None
        )
        self.debug_write_diamond_output = debug_write_diamond_output
        self.run_diamond_enabled = run_diamond
        self._diamond_output_saved = None
        self.quiet_mode = quiet_mode
        self.mem_mode = mem_mode
        self.background_rate: float = background_rate

    def _validate_inputs(self) -> None:
        if not self.fasta_file.is_file():
            raise InputError(f"Input sequence file not found: {self.fasta_file}")
        if self.run_kraken_enabled:
            if not self.kraken_db.is_dir():
                raise InputError(f"Kraken2 database directory not found: {self.kraken_db}")
        else:
            if self.kraken_output_path is None:
                raise InputError(
                    "kraken_output_path is required when run_kraken=False."
                )
            if not self.kraken_output_path.is_file():
                raise InputError(
                    f"Kraken classifications file not found: {self.kraken_output_path}"
                )

    def _ensure_kraken_in_memory(self) -> None:
        if self._kraken_lines is not None:
            return
        if self.kraken_output_path and self.kraken_output_path.exists():
            try:
                lines = self.kraken_output_path.read_text(encoding="utf-8").splitlines(
                    True
                )
            except OSError as exc:
                raise InputError(
                    f"Could not read Kraken classifications {self.kraken_output_path}: {exc}"
                ) from exc
            self._kraken_lines = lines
            logging.info(
                "Loaded %d Kraken classification line(s) from %s",
                len(lines),
                self.kraken_output_path,
            )
            return
        raise InputError(
            "No Kraken classifications available. Run Kraken2 or provide "
            "an existing --kraken-output-path."
        )

    def run_kraken(self) -> None:
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

        try:
            if self.debug_write_kraken_output:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                self._ensure_kraken_in_memory()
            else:
                proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
                lines = proc.stdout.splitlines(True)
                self._kraken_lines = lines
        except FileNotFoundError as e:
            logging.error("Kraken2 executable not found on PATH")
            raise KrakenError(
                127,
                "kraken2 executable not found; is Kraken2 installed and on PATH?",
            ) from e
        except subprocess.CalledProcessError as e:
            logging.error(
                "Kraken2 failed (Code %s). Error details: %s", e.returncode, e.stderr
            )
            raise KrakenError(e.returncode, e.stderr) from e
        logging.info("Kraken2 finished; classifications ready.")

    @staticmethod
    def parse_and_run(
            kmer_pos_info: str, window_size: int, threshold: int, k: int = 35, l: int = 31,
    ):
        """
        Parse a Kraken k-mer string and scan for engineered (taxid 32630) windows.

        Returns
        -------
        (hit, max_eng, target_kmers, total_mapped_kmers, eng_kmer_coverage)
        """
        empty = (False, 0, 0, 0, 0.0)
        raw_data = kmer_pos_info.replace("|:|", "").split()
        if not raw_data:
            return empty

        window_size = int(window_size)
        threshold = int(threshold)

        max_minimizers = max(window_size - k - l + 2, 1)
        target_tid = "32630"
        n = len(raw_data)
        tids = np.zeros(n)
        counts = np.zeros(n)
        target_kmers = 0
        total_kmers_including_unmapped = 0
        total_mapped_kmers = 0
        try:
            for i, item in enumerate(raw_data):
                t, c = item.split(":", 1)
                if t == "A":
                    tids[i] = -1
                else:
                    tids[i] = int(t)
                    total_mapped_kmers += int(c)
                total_kmers_including_unmapped += int(c)
                counts[i] = min(int(c), max_minimizers)
                if target_tid == t:
                    target_kmers += counts[i]

        except (ValueError, IndexError, OverflowError):
            return empty

        coverage = (
            target_kmers / total_kmers_including_unmapped
            if total_kmers_including_unmapped
            else 0.0
        )
        tids = np.ascontiguousarray(tids, dtype=np.int64)
        counts = np.ascontiguousarray(counts, dtype=np.int64)
        try:
            hit, max_eng = fast_window_logic(tids, counts, window_size, threshold)
        except TypeError:
            # Older Numba builds can fail to unbox arrays; fall back to pure Python.
            hit, max_eng = Workflow._fast_window_logic_python(
                tids, counts, window_size, threshold
            )
        return hit, max_eng, target_kmers, total_mapped_kmers, coverage

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

    def entry_logic(self, entry: str, kmer_max_by_read, result):
        categories = entry.split("\t")
        if len(categories) < 2:
            return 0.0, 0.0
        synthetic_boolean, max_eng, target_kmers, tot_kmers, eng_kmer_coverage = self.parse_and_run(
            categories[-1], self.window_size, self.threshold, self.k, self.l)
        read_id = categories[1]
        kmer_max_by_read[read_id] = max_eng
        if synthetic_boolean:
            result.synthetic_count += 1
            label = "Synthetic"
        else:
            result.natural_count += 1
            label = "Natural"
        result.labels.append(ReadEngineeringLabel(
            read_id=read_id,
            label=label,
            eng_kmer_coverage=eng_kmer_coverage,
        ))
        return target_kmers, tot_kmers

    def scan_engineered_blocks_kraken(
            self,
    ) -> tuple[EngineeredScanResult, NDArray, dict[str, int]]:
        result = EngineeredScanResult()
        self._ensure_kraken_in_memory()
        kmer_max_by_read: dict[str, int] = {}
        entries = self._kraken_lines or []
        target_kmers = np.zeros(len(entries))
        tot_kmers = np.zeros(len(entries))
        index = 0
        if not self.quiet_mode:
            for entry in tqdm(entries, total=len(entries)):
                target_kmers[index], tot_kmers[index] = self.entry_logic(
                    entry, kmer_max_by_read, result
                )
                index += 1
        else:
            for entry in entries:
                target_kmers[index], tot_kmers[index] = self.entry_logic(
                    entry, kmer_max_by_read, result
                )
                index += 1
        logging.info(
            "Engineered k-mer scan complete: %s/%s synthetic reads.",
            result.synthetic_count,
            result.synthetic_count + result.natural_count,
        )
        p_values, _qvals, _syn = call_synthetic_reads(
            target_kmers,
            tot_kmers,
            self.background_rate,
            n_reads_total=len(entries),
            alpha=0.05,
        )
        return result, p_values, kmer_max_by_read

    def run_codon_adaptation(self, natural_read_ids: set[str]) -> CodonAdaptationResult:
        if not natural_read_ids:
            logging.info("No reads labeled Natural; skipping codon usage analysis.")
            return CodonAdaptationResult(labels=[])

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
            mem_mode=self.mem_mode,
            diamond_threads=self.diamond_threads,
            run_diamond=self.run_diamond_enabled,
            diamond_output_path=self.diamond_output_path,
            debug_write_diamond_output=self.debug_write_diamond_output,
            codon_usage_dir=self.codon_usage_dir,
            include_read_ids=natural_read_ids,
            seq_type=self.seq_type,
        )
        self._diamond_output_saved = diamond_path
        return results

    def write_screen_result(self, per_read: List[ReadFlagDetail]) -> None:
        if self.report_output_path is None:
            return
        header = "Label\tRead_ID\tMethods\tP-Value\tEng-Coverage Percentage"
        try:
            self.report_output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.report_output_path, "w", encoding="utf-8") as write_obj:
                write_obj.write(header + "\n")
                for read in per_read:
                    if read.engineered_overall:
                        line = "Engineered\t"
                    else:
                        line = "Natural\t"
                    line += read.read_id + "\t"
                    if read.engineered_by_kmer_scan:
                        line += "engineered_kmer_scan"
                    elif read.engineered_by_codon_cai:
                        line += "codon_cai"
                    else:
                        line += "NA"
                    line += "\t" + str(read.p_value)
                    coverage = (
                        read.eng_kmer_coverage
                        if read.eng_kmer_coverage is not None
                        else 0.0
                    )
                    line += "\t" + f"{coverage:.2%}"
                    write_obj.write(line + "\n")
        except OSError as exc:
            logging.error("Failed to write report %s: %s", self.report_output_path, exc)
            raise InputError(
                f"Could not write engineered report {self.report_output_path}: {exc}"
            ) from exc
        logging.info("Wrote engineered report to %s", self.report_output_path)

    def run(self) -> ScreenResult:
        self._validate_inputs()
        if self.run_kraken_enabled:
            self.run_kraken()
        engineered_scan, p_values, kmer_max_by_read = self.scan_engineered_blocks_kraken()

        codon_path: Path | None = None
        codon_results = CodonAdaptationResult(labels=[])
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

        codon_by_read = {r.read_id: r for r in codon_results.labels}

        per_read: list[ReadFlagDetail] = []
        for index, lbl in enumerate(engineered_scan.labels):
            codon = codon_by_read.get(lbl.read_id)
            cai = codon.cai_vs_host if codon else None
            engineered_by_codon: bool | None = None
            if cai is not None and self.codon_cai_engineered_threshold is not None:
                # CAI above threshold flags engineered (codon optimization signal).
                engineered_by_codon = cai > self.codon_cai_engineered_threshold

            engineered_by_kmer = lbl.label == "Synthetic"
            engineered_overall = compute_engineered_overall(
                engineered_by_kmer_scan=engineered_by_kmer,
                engineered_by_codon_cai=engineered_by_codon,
            )
            overall_label: Literal["Natural", "Synthetic"] = (
                "Synthetic" if engineered_overall else "Natural"
            )
            kmer_label: Literal["Natural", "Synthetic"] = (
                "Natural" if p_values[index] > 0.05 else lbl.label
            )

            per_read.append(
                ReadFlagDetail(
                    read_id=lbl.read_id,
                    kmer_label=kmer_label,
                    eng_kmer_coverage=lbl.eng_kmer_coverage,
                    engineered_by_kmer_scan=engineered_by_kmer,
                    engineered_overall=engineered_overall,
                    overall_label=overall_label,
                    engineered_kmer_max_in_window=kmer_max_by_read.get(lbl.read_id),
                    engineered_kmer_threshold=self.threshold,
                    engineered_kmer_window_size=self.window_size,
                    cai_vs_host=cai,
                    engineered_by_codon_cai=engineered_by_codon,
                    codon_cai_threshold=self.codon_cai_engineered_threshold,
                    p_value=float(p_values[index]),
                )
            )

        if self.report_output_path is not None:
            self.write_screen_result(per_read)

        return ScreenResult(
            engineered_scan=engineered_scan,
            codon_adaptation=codon_results,
            per_read=per_read,
            engineered_report_path=self.report_output_path,
            codon_usage_report_path=codon_path,
            diamond_output_path=self._diamond_output_saved,
            engineered_kmer_threshold=self.threshold,
            engineered_kmer_window_size=self.window_size,
            codon_cai_engineered_threshold=self.codon_cai_engineered_threshold,
        )
