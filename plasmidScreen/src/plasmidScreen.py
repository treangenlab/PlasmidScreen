"""

"""
import logging
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm
import numpy as np
try:
    from numba import jit
except ModuleNotFoundError:  # pragma: no cover
    def jit(*_args, **_kwargs):  # type: ignore
        def _wrap(fn):
            return fn
        return _wrap

from plasmidScreen.lib.codon_usage_db import default_codon_usage_dir
from plasmidScreen.lib.models import (
    CodonAdaptationResult,
    EngineeredScanResult,
    ReadFlagDetail,
    ReadEngineeringLabel,
    ScreenResult,
)
from plasmidScreen.src.analyze_codon_usage import (
    analyze_codon_adaptation,
    write_codon_adaptation_tsv,
)


@jit(nopython=True)
def fast_window_logic(tids, counts, window_size: int, threshold: int):
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

    def __init__(
        self,
        fasta_file: str,
        output_report_path: str,
        kraken_db: str,
        threads: int,
        kraken_raw_output: str,
        window_size: int = 200,
        engineered_kmer_threshold: int = 15,
        codon_usage_output_path: str | None = None,
        codon_usage_dir: str | None = None,
        run_kraken: bool = True,
        run_codon_usage: bool = True,
        kmer_len: int = 35,
        codon_cai_engineered_threshold: float | None = None,
    ):
        self.fasta_file: Path = Path(fasta_file)
        self.report_output_path = Path(output_report_path)
        self.kraken_db = Path(kraken_db)
        self.kraken_output_path = Path(kraken_raw_output)
        self.threshold = engineered_kmer_threshold
        self.window_size = window_size
        self.max_threads = threads
        self.codon_usage_output_path = Path(codon_usage_output_path) if codon_usage_output_path else None
        self.codon_usage_dir = (
            Path(codon_usage_dir) if codon_usage_dir else default_codon_usage_dir()
        )
        self.run_kraken = run_kraken
        self.run_codon_usage = run_codon_usage
        self.kmer_len = kmer_len
        self.codon_cai_engineered_threshold = codon_cai_engineered_threshold

    def run_kraken(self) -> None:
        try:
            logging.info(f"Running Kraken on {self.fasta_file}")
            report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
            logging.info(
                f"Running the command: kraken2 --db {self.kraken_db} --report-minimizer-data "
                f"--report {report_file_path} --output {self.kraken_output_path} --use-names {self.fasta_file} "
                f"--threads {self.max_threads}"
            )
            subprocess.run(
                [
                    "kraken2",
                    "--db",
                    str(self.kraken_db),
                    "--report-minimizer-data",
                    "--report",
                    str(report_file_path),
                    "--output",
                    str(self.kraken_output_path),
                    "--use-names",
                    str(self.fasta_file),
                    "--threads",
                    str(self.max_threads),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Kraken2 failed (Code {e.returncode}). Error details: {e.stderr}")
            sys.exit(1)
        logging.info("All files processed.")

    @staticmethod
    def parse_and_run(kmer_pos_info: str, window_size: int, threshold: int) -> tuple[bool, int]:
        raw_data = kmer_pos_info.replace("|:|", "").split()

        n = len(raw_data)
        tids = np.zeros(n, dtype=np.int32)
        counts = np.zeros(n, dtype=np.int32)

        for i, item in enumerate(raw_data):
            t, c = item.split(":")
            if t == "A":
                tids[i] = -1
            else:
                tids[i] = int(t)
            counts[i] = int(c)

        return fast_window_logic(tids, counts, window_size, threshold)

    def scan_engineered_blocks_kraken(self) -> EngineeredScanResult:
        result = EngineeredScanResult()
        with open(self.kraken_output_path, "r") as kraken_file, open(
            self.report_output_path, "w"
        ) as report:
            report.write("Label\tRead_ID\tMethods\tEngineeredKmerMaxInWindow\tKmerThreshold\tWindowSize\n")
            entries = kraken_file.readlines()
            for entry in tqdm(entries, total=len(entries)):
                categories = entry.split("\t")
                synthetic_boolean, max_eng = self.parse_and_run(
                    categories[-1], self.window_size, self.threshold
                )
                read_id = categories[1]
                if synthetic_boolean:
                    result.synthetic_count += 1
                    label = "Synthetic"
                else:
                    result.natural_count += 1
                    label = "Natural"
                result.labels.append(ReadEngineeringLabel(read_id=read_id, label=label))
                methods = "engineered_kmer_scan" if label == "Synthetic" else ""
                report.write(
                    f"{label}\t{read_id}\t{methods}\t{max_eng}\t{self.threshold}\t{self.window_size}\n"
                )

        logging.info(
            f"Engineered k-mer scan complete: {result.synthetic_count}/"
            f"{result.synthetic_count + result.natural_count} synthetic reads."
        )
        return result

    def run_codon_adaptation(self, natural_read_ids: set[str]) -> list[CodonAdaptationResult]:
        if not natural_read_ids:
            logging.info("No reads labeled Natural; skipping codon usage analysis.")
            return []

        return analyze_codon_adaptation(
            self.fasta_file,
            self.kraken_output_path,
            codon_usage_dir=self.codon_usage_dir,
            include_read_ids=natural_read_ids,
            kmer_len=self.kmer_len,
        )

    def run(self) -> ScreenResult:
        if self.run_kraken:
            self.run_kraken()

        engineered_scan = self.scan_engineered_blocks_kraken()

        codon_results: list[CodonAdaptationResult] = []
        codon_path: Path | None = None

        if self.run_codon_usage and engineered_scan.natural_read_ids:
            if engineered_scan.any_synthetic:
                logging.info(
                    "Engineered reads detected; running codon usage only on reads labeled Natural."
                )

            codon_path = self.codon_usage_output_path
            if codon_path is None:
                codon_path = self.report_output_path.with_suffix(
                    self.report_output_path.suffix + ".codon_usage.tsv"
                )

            logging.info(
                f"Running codon usage analysis on {len(engineered_scan.natural_read_ids)} "
                f"Natural reads (reference: {self.codon_usage_dir})"
            )
            codon_path_str, codon_results = write_codon_adaptation_tsv(
                self.fasta_file,
                self.kraken_output_path,
                codon_path,
                include_read_ids=engineered_scan.natural_read_ids,
                codon_usage_dir=self.codon_usage_dir,
                kmer_len=self.kmer_len,
                cai_engineered_threshold=self.codon_cai_engineered_threshold,
            )
            codon_path = Path(codon_path_str)

        codon_by_read: dict[str, CodonAdaptationResult] = {r.read_id: r for r in codon_results}

        # k-mer evidence lookup (max engineered k-mers observed in any window)
        kmer_max_by_read: dict[str, int] = {}
        with open(self.kraken_output_path, "r") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                rid = parts[1]
                hit, max_eng = self.parse_and_run(parts[-1], self.window_size, self.threshold)
                kmer_max_by_read[rid] = max_eng

        per_read: list[ReadFlagDetail] = []
        for lbl in engineered_scan.labels:
            codon = codon_by_read.get(lbl.read_id)
            cai = codon.cai_vs_host if codon else None
            engineered_by_codon: bool | None = None
            if cai is not None and self.codon_cai_engineered_threshold is not None:
                engineered_by_codon = cai < self.codon_cai_engineered_threshold

            per_read.append(
                ReadFlagDetail(
                    read_id=lbl.read_id,
                    kmer_label=lbl.label,
                    engineered_by_kmer_scan=(lbl.label == "Synthetic"),
                    engineered_kmer_max_in_window=kmer_max_by_read.get(lbl.read_id),
                    engineered_kmer_threshold=self.threshold,
                    engineered_kmer_window_size=self.window_size,
                    cai_vs_host=cai,
                    engineered_by_codon_cai=engineered_by_codon,
                    codon_cai_threshold=self.codon_cai_engineered_threshold,
                )
            )

        return ScreenResult(
            engineered_scan=engineered_scan,
            codon_adaptation=codon_results,
            per_read=per_read,
            engineered_report_path=self.report_output_path,
            codon_usage_report_path=codon_path,
        )
