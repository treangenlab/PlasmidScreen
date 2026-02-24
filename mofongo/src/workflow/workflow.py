"""

"""
import abc
import logging
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm
from collections import deque
from mofongo.lib.const import KrakenConfig, MetaMDBGConfig, MegahitConfig
from collections import deque
import numpy as np
from numba import jit


@jit(nopython=True)
def fast_window_logic(tids, counts):
    window_size = 500
    overlap_max = 25

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
            if eng_count >= overlap_max:
                return True

    return False


class Workflow(abc.ABC):

    def __init__(self, fasta_file: str, output_report_path: str, kraken_db: str, kraken_raw_output: str,
                 kraken_config: KrakenConfig,
                 assembly: bool, assembly_config):
        self.fasta_file: Path = Path(fasta_file)
        self.report_output_path = Path(output_report_path)
        self.kraken_db = Path(kraken_db)
        self.kraken_config = kraken_config
        self.assembly_config = assembly_config
        self.kraken_output_path = kraken_raw_output
        if assembly:
            self.run_assembly()

    @abc.abstractmethod
    def run_assembly(self):
        raise NotImplementedError

    def run_kraken(self, ):
        #try:
        #    logging.info(f"Running Kraken on {self.fasta_file}")
                            #output_file_path: Path = self.kraken_db.parent.joinpath("output.txt")
        #    report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
        #    logging.info(f"Running the command: kraken2 --db {self.kraken_db} --report-minimizer-data "
        #                 f"--report {report_file_path} --output {self.kraken_output_path} --use-names {self.fasta_file} "
        #                 f"--threads {self.kraken_config.max_threads}")
        #    subprocess.check_call(
        #        ["kraken2", "--db", self.kraken_db, "--report-minimizer-data", "--report", report_file_path,
        #         "--output", self.kraken_output_path, "--use-names", self.fasta_file,
        #         "--threads", str(self.kraken_config.max_threads)])
        #except subprocess.CalledProcessError:
        #    logging.info(f"CRITICAL ERROR: Failed to run kraken2 {self.fasta_file}.")
        #    sys.exit(1)
        logging.info("All files processed.")
        self._post_process_kraken()

    @staticmethod
    def _check_windows_old(kmer_pos_info: str) -> bool:
        window_size = 45
        overlap_max = 25
        max_kmers_in_window = window_size - 21 + 1  # 280

        # Pre-parse segments into a list of tuples (is_engineered, count)
        # Using '32630' comparison during parsing saves time later
        segments = []
        #buffer_max = 10
        #window_buffer = 10
        kmers = kmer_pos_info.split(' ')
        for kmer_index, item in enumerate(kmers):
            if item == "|:|": continue
            try:
                tid, count_str = item.split(':')
                if tid == "32630":
                    #           preceding_hits = 0
                    #           bad_hits_prev = 0
                    #           while preceding_hits < window_buffer:
                    #               preceeding_tid, preceding_count = kmers[kmer_index - 1].split(':')
                    #               preceding_hits+=int(preceding_count)
                    #               if preceeding_tid == "1" or preceeding_tid=="2":
                    #                   bad_hits_prev+=int(preceding_count)
                    #           post_hits = 0
                    #           bad_hits_post = 0
                    #           while post_hits < window_buffer:
                    #               following_tid, following_count = kmers[kmer_index + 1].split(':')
                    #               post_hits += int(following_count)
                    #               if following_tid == "1" or following_tid == "2":
                    #                   bad_hits_post += int(following_count)

                    #preceeding_condition = ((
                    #                               preceeding_tid == "1" or preceeding_tid == "0")
                    #                     and int(preceding_count) > buffer_max)
                    #following_condition = (
                    #                                  following_tid == "1" or following_tid == "0") and int(following_count) > buffer_max
                    segments.append((True, int(count_str)))

                #            segments.append((bad_hits_prev < buffer_max and bad_hits_post < buffer_max, int(count_str)))
                else:
                    segments.append((False, int(count_str)))
            except ValueError:
                continue
            except IndexError as ex:
                segments.append((tid == "32630", int(count_str)))

        eng_count = 0
        non_eng_count = 0
        total_in_window = 0

        for is_eng, count in segments:
            # If the segment is long enough to hit the threshold alone, return early
            if is_eng and count >= overlap_max:
                return True

            # We simulate the sliding window behavior without the loop.
            # If we add 'count' kmers of type X, they replace 'count' kmers of type Y
            # (or X) that were previously in the window.

            if is_eng:
                if total_in_window < max_kmers_in_window:
                    added = min(count, max_kmers_in_window - total_in_window)
                    eng_count += added
                    total_in_window += added
                    # If there's still count left after filling the window
                    remaining = count - added
                    if remaining > 0:
                        # They replace old kmers. In the worst case for speed,
                        # they replace non-engineered ones.
                        replaced_non_eng = min(remaining, non_eng_count)
                        eng_count += replaced_non_eng
                        non_eng_count -= replaced_non_eng
                else:
                    # Window is full, these kmers replace whatever was there.
                    # To be safe/conservative (matching your logic):
                    replaced_non_eng = min(count, non_eng_count)
                    eng_count = min(max_kmers_in_window, eng_count + replaced_non_eng)
                    non_eng_count -= replaced_non_eng
            else:
                if total_in_window < max_kmers_in_window:
                    added = min(count, max_kmers_in_window - total_in_window)
                    non_eng_count += added
                    total_in_window += added
                    remaining = count - added
                    if remaining > 0:
                        replaced_eng = min(remaining, eng_count)
                        non_eng_count += replaced_eng
                        eng_count -= replaced_eng
                else:
                    replaced_eng = min(count, eng_count)
                    non_eng_count = min(max_kmers_in_window, non_eng_count + replaced_eng)
                    eng_count -= replaced_eng

            if eng_count >= overlap_max:
                return True

        return False

    @staticmethod
    def _check_windows(kmer_pos_info: str) -> bool:
        """

        :return:
        """
        # 500 at 15 gives 0.7475
        # 450 at 25 gives 0.7 F1
        window_size = 200

        segments = []
        for item in kmer_pos_info.split(' '):
            if item == "|:|": continue  # Skip paired-end delimiters
            #try:
            tid, count = item.split(':')
            segments.append((str(tid), int(count)))
            #except ValueError:
            #    continue

        # Reset for the optimized approach
        current_kmer_idx = 0

        # Dictionary to hold counts: window_index -> count
        #window_counts = {}

        # 0:1
        overlap_max = 20
        #fixed_size_buffer = deque(maxlen=100)
        max_kmers_in_window = window_size - 21 + 1
        current_kmers_in_window = {"count": 0, "engineered": 0, "non-engineered": 0}
        segs = iter(segments)
        for tid, count in segs:
            seg_start = current_kmer_idx + 21 - 1
            seg_end = current_kmer_idx + count + 21 - 1
            for window in range(seg_start, seg_end, 1):
                if current_kmers_in_window["count"] < max_kmers_in_window:
                    if tid == "32630":
                        current_kmers_in_window["engineered"] += 1
                    else:
                        current_kmers_in_window["non-engineered"] += 1
                    current_kmers_in_window["count"] += 1
                else:
                    if tid == "32630":
                        current_kmers_in_window["engineered"] = min(max_kmers_in_window,
                                                                    current_kmers_in_window["engineered"] + 1)
                        current_kmers_in_window["non-engineered"] = max(current_kmers_in_window["non-engineered"] - 1,
                                                                        0)
                    else:
                        current_kmers_in_window["non-engineered"] = min(max_kmers_in_window,
                                                                        current_kmers_in_window["non-engineered"] + 1)
                        current_kmers_in_window["engineered"] = max(current_kmers_in_window["engineered"] - 1, 0)
                        #fixed_size_buffer.append(tid)
                if current_kmers_in_window["engineered"] >= overlap_max:
                    #try:
                    #tid,count = next(segs)
                    #print(f"{tid}:{count}")
                    #running_count=0
                    #valid_hits = 0
                    # while running_count < 100:
                    #     running_count += count
                    #     if tid=="0" or tid=="1":
                    #         valid_hits += count
                    #     tid, count = next(segs)

                    # if valid_hits>=100:
                    #     return False

                    #  if (tid == "0" or tid == "1") and count >= 15 :
                    #      return False
                    #  tid_count=0
                    #  for item in fixed_size_buffer:
                    #      if item =="1" or item =="0":
                    #          tid_count+=1
                    #  if tid_count>=100:
                    #      return False

                    # print("hit")
                    # break
                    # except StopIteration as ex:
                    #     print("oops")
                    return True
        return False

    #if len(current_kmers_in_window) < max_kmers_in_window:
    #    current_kmers_in_window.append(tid)
    #else:
    #    current_kmers_in_window.popleft()
    #    current_kmers_in_window.append(tid)
    #if

    # If this segment is our target, add counts to appropriate windows
    #if tid == "32630":  #"1001":
    # Find which windows this segment touches
    #first_window_idx = seg_start // window_size

    # last_window_idx = seg_end  // window_size

    #last_window_idx = (seg_end - 1) // window_size

    #    for w_idx in range(first_window_idx, last_window_idx + 1):
    # Calculate overlap with this specific window
    #        w_start_limit = w_idx * window_size
    #        w_end_limit = (w_idx + 1) * window_size

    #        overlap = max(0, min(seg_end, w_end_limit) - max(seg_start, w_start_limit))
    #        if overlap >= 15:
    #            return True

    # current_kmer_idx += count
    #return False

    # 1. Parse data FIRST (Python side)

    # 2. The Logic Loop (Compiled by Numba)
    # nopython=True means "If you can't compile this to C, fail." (Ensures speed)

    @staticmethod
    def parse_and_run(kmer_pos_info: str):
        # Quick string parsing to get integer arrays
        # This replaces the slow item.split(':') inside the loop
        raw_data = kmer_pos_info.replace('|:|', '').split()

        # Pre-allocate arrays for speed
        n = len(raw_data)
        tids = np.zeros(n, dtype=np.int32)
        counts = np.zeros(n, dtype=np.int32)

        for i, item in enumerate(raw_data):
            t, c = item.split(':')
            if t=='A':
                tids[i] = -1
            else:
                tids[i] = int(t)
            counts[i] = int(c)

        return fast_window_logic(tids, counts)

    @abc.abstractmethod
    def _rescue_logic(self) -> bool:
        raise NotImplementedError

    def _post_process_kraken(self):
        report = open(self.report_output_path, "w")
        with open(self.kraken_output_path, 'r') as kraken_file:
            entries = kraken_file.readlines()
            for entry in tqdm(entries, total=len(entries)):
                categories = entry.split("\t")
                synthetic_boolean: bool = self.parse_and_run(categories[-1])
                #synthetic_boolean: bool = self._check_windows(categories[-1])
                logging.debug(f"Entry: {categories[1]} found to be synthetic from windowing logic: {synthetic_boolean}")
                if synthetic_boolean:
                    #contig_path = Path(self.kraken_output_path).parent.joinpath("final.contigs.fa")
                    #out = subprocess.run(["grep", categories[1], contig_path], capture_output=True).stdout
                    #cov = str(out).split(" ")[-1].split("=")[-1]
                    #cov = float(cov.strip("\\n'"))
                    #if cov >= 300:
                    report.write(f"Synthetic\t{categories[1]}\n")
                    #else:
                    #    report.write(f"Natural\t{categories[1]}\n")
                else:
                    synthetic_boolean = self._rescue_logic()
                    logging.debug(
                        f"Entry: {categories[1]} found to be synthetic from rescue logic: {synthetic_boolean}")
                    if synthetic_boolean:
                        report.write(f"Synthetic\t{categories[1]}\n")
                    else:
                        report.write(f"Natural\t{categories[1]}\n")

    def run(self):
        self.run_kraken()
        self._post_process_kraken()


class LRWorkflow(Workflow):

    def _rescue_logic(self) -> bool:
        return False

    def run_assembly(self):
        try:
            logging.info(f"Running metaMDBG on {self.fasta_file}")
            output_file_path: Path = self.kraken_db.parent.joinpath("output.txt")
            report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
            cmd = ["metaMDBG", "asm", "--report-minimizer-data", "--report",
                   report_file_path, "--output", output_file_path, "--use-names",
                   self.fasta_file, "--threads", self.kraken_config.max_threads]
            logging.info(f"Running the command: metaMDBG asm  --report-minimizer-data "
                         f"--report {report_file_path} --output {output_file_path} --use-names {self.fasta_file} "
                         f"--threads {self.kraken_config.max_threads}")
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            logging.info(f"CRITICAL ERROR: Failed to run kraken2 {self.fasta_file}.")
            sys.exit(1)
        logging.info("All files processed.")


class SRWorkflow(Workflow):

    def __init__(self, fasta_file: str, paired_end: str, output_report_path: str, kraken_db: str,
                 kraken_raw_output: str,
                 kraken_config: KrakenConfig,
                 assembly: bool, assembly_config: MegahitConfig):
        super().__init__(fasta_file, output_report_path, kraken_db, kraken_raw_output, kraken_config, assembly,
                         assembly_config)
        self.paired_end = paired_end

    def _rescue_logic(self) -> bool:
        return False

    def run_assembly(self):
        try:
            logging.info(f"Running megahit on {self.fasta_file}")
            if self.assembly_config.paired:
                fastq_files_cmd = ["-1", self.fasta_file, "-2", self.paired_end]
            else:
                fastq_files_cmd = ["-1", self.fasta_file]

                #output_file_path: Path = self.kraken_db.parent.joinpath("output.txt")
            report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
            cmd = ["megahit", "--num-cpu-threads", self.kraken_config.max_threads,
                   "--output_dir", self.report_output_path]
            cmd.extend(fastq_files_cmd)
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            logging.info(f"CRITICAL ERROR: Failed to run kraken2 {self.fasta_file}.")
            sys.exit(1)
        logging.info("All files processed.")
