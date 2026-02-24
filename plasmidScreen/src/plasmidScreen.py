"""

"""
import logging
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm
import numpy as np
from numba import jit


@jit(nopython=True)
def fast_window_logic(tids, counts, window_size: int, threshold: int):
    #window_size = 200
    #threshold = 20

    # window_size = 500
    target_tid = 32630  # Use int for speed comparisons
    # overlap_max = 15
    max_kmers = window_size - 21 + 1

    # Local variables (Faster than dicts)
    eng_count = 0
    non_eng_count = 0
    total_count = 0

    for i in range(len(tids)):
        tid = tids[i]
        count = counts[i]

        # The inner loop, now running at C-speed
        for _ in range(count):
            if total_count < max_kmers:
                if tid == target_tid:
                    eng_count += 1
                else:
                    non_eng_count += 1
                total_count += 1
            else:
                if tid == target_tid:
                    # Logic directly translated, but fast
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

            # Check threshold
            if eng_count >= threshold:
                return True

    return False


class Workflow:

    def __init__(self, fasta_file: str, output_report_path: str, kraken_db: str, threads: int, kraken_raw_output: str,
                 window_size: int = 200, engineered_kmer_threshold: int = 15):
        self.fasta_file: Path = Path(fasta_file)
        self.report_output_path = Path(output_report_path)
        self.kraken_db = Path(kraken_db)
        self.kraken_output_path = kraken_raw_output
        self.threshold = engineered_kmer_threshold
        self.window_size = window_size
        self.max_threads = threads

    def run_kraken(self, ):
        try:
            logging.info(f"Running Kraken on {self.fasta_file}")
            report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
            logging.info(f"Running the command: kraken2 --db {self.kraken_db} --report-minimizer-data "
                         f"--report {report_file_path} --output {self.kraken_output_path} --use-names {self.fasta_file} "
                         f"--threads {self.max_threads}")
            subprocess.run(
                ["kraken2", "--db", self.kraken_db, "--report-minimizer-data", "--report", report_file_path,
                 "--output", self.kraken_output_path, "--use-names", self.fasta_file,
                 "--threads", str(self.max_threads)], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            # Now you can log the exact error message Kraken2 spit out
            logging.error(f"Kraken2 failed (Code {e.returncode}). Error details: {e.stderr}")
            sys.exit(1)
        logging.info("All files processed.")

    @staticmethod
    def parse_and_run(kmer_pos_info: str, window_size: int, threshold: int):
        # Quick string parsing to get integer arrays
        # This replaces the slow item.split(':') inside the loop
        raw_data = kmer_pos_info.replace('|:|', '').split()

        # Pre-allocate arrays for speed
        n = len(raw_data)
        tids = np.zeros(n, dtype=np.int32)
        counts = np.zeros(n, dtype=np.int32)

        for i, item in enumerate(raw_data):
            t, c = item.split(':')
            if t == 'A':
                tids[i] = -1
            else:
                tids[i] = int(t)
            counts[i] = int(c)

        return fast_window_logic(tids, counts, window_size, threshold)

    def scan_engineered_blocks_kraken(self):
        report = open(self.report_output_path, "w")
        with open(self.kraken_output_path, 'r') as kraken_file:
            entries = kraken_file.readlines()
            for entry in tqdm(entries, total=len(entries)):
                categories = entry.split("\t")
                synthetic_boolean: bool = self.parse_and_run(categories[-1], self.window_size, self.threshold)
                logging.debug(f"Entry: {categories[1]} found to be synthetic from windowing logic: {synthetic_boolean}")
                if synthetic_boolean:
                    report.write(f"Synthetic\t{categories[1]}\n")
                else:
                    report.write(f"Natural\t{categories[1]}\n")

    def run(self):
        self.run_kraken()
        self.scan_engineered_blocks_kraken()
