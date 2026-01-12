import abc
import logging
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm

from mofongo.lib.const import KrakenConfig


class Workflow(abc.ABC):

    def __init__(self, fasta_file: str, kraken_db: str, kraken_config: KrakenConfig):
        self.fasta_file: Path = Path(fasta_file)
        self.kraken_db = Path(kraken_db)
        self.kraken_config = kraken_config

    def run_kraken(self, ):
        try:
            logging.info(f"Running Kraken on {self.fasta_file}")
            output_file_path: Path = self.kraken_db.parent.joinpath("output.txt")
            report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
            logging.debug(f"Running the command: kraken2 --db {self.kraken_db} --report-minimizer-data "
                          f"--report {report_file_path} --output {output_file_path} --use-names {self.fasta_file}")
            subprocess.check_call(
                ["kraken2", "--db", self.kraken_db, "--report-minimzer-data", "--report", report_file_path,
                 "--output", output_file_path, "--use-names", self.fasta_file])
        except subprocess.CalledProcessError:
            logging.info(f"CRITICAL ERROR: Failed to run kraken2 {self.fasta_file}.")
            sys.exit(1)
        logging.info("All files processed.")

    @staticmethod
    def _check_windows(kmer_pos_info: str) -> bool:
        """

        :return:
        """
        window_size_to_check = 0
        current_bp = 0
        pairs = kmer_pos_info.split(' ')
        number_of_synthetic_kmers_in_window = 0
        window_size = 200

        segments = []
        for item in kmer_pos_info.split(' '):
            if item == "|:|": continue  # Skip paired-end delimiters
            try:
                tid, count = item.split(':')
                segments.append((tid, int(count)))
            except ValueError:
                continue

        # Iterate through fixed-size windows for this sequence
        # Total k-mers = seq_len - k_len + 1
        #total_kmers = sum(count for _, count in segments)

        # We process windows based on K-MER INDICES, not raw bases,
        # because Kraken classifies k-mers.
        # Window 0 = k-mers 0 to 1000
        # for win_start in range(0, total_kmers, window_size):
        #     win_end = win_start + window_size
        #     target_hits_in_window = 0

        # We need to track where we are in the k-mer stream
        #     current_kmer_idx = 0

        #     for tid, count in segments:
        # Define the span of this specific taxID segment
        #         seg_start = current_kmer_idx
        #         seg_end = current_kmer_idx + count

        # Check for overlap between this segment and the current window
        #         # Overlap logic: max(0, min(end1, end2) - max(start1, start2))
        #         overlap = max(0, min(seg_end, win_end) - max(seg_start, win_start))

        #         if overlap > 0 and tid == "Synthetic (taxid 1001)":
        #             target_hits_in_window += overlap

        # Advance the stream
        #         current_kmer_idx += count

        # Optimization: If we passed the window, break early for this window
        #         if current_kmer_idx > win_end:
        #             break  # But we must not break the outer loop (windows), just this segment loop?
        # Actually, re-iterating segments for every window is slow O(N*M).
        # Let's optimize the loop below.

        # --- OPTIMIZED SINGLE PASS APPROACH ---
        # Instead of nested loops, we walk through segments and fill windows.

        # Reset for the optimized approach
        current_kmer_idx = 0

        # Dictionary to hold counts: window_index -> count
        # window_counts = {}

        for tid, count in segments:
            seg_start = current_kmer_idx
            seg_end = current_kmer_idx + count

            # If this segment is our target, add counts to appropriate windows
            if tid == "Synthetic (taxid 1001)":
                # Find which windows this segment touches
                first_window_idx = seg_start // window_size
                last_window_idx = (seg_end - 1) // window_size

                for w_idx in range(first_window_idx, last_window_idx + 1):
                    # Calculate overlap with this specific window
                    w_start_limit = w_idx * window_size
                    w_end_limit = (w_idx + 1) * window_size

                    overlap = max(0, min(seg_end, w_end_limit) - max(seg_start, w_start_limit))
                    if overlap >= 25:
                        return True
                    #if w_idx not in window_counts: window_counts[w_idx] = 0
                    #window_counts[w_idx] += overlap

            current_kmer_idx += count

        # Report results for this sequence
        # Calculate total windows based on total kmers
        #num_windows = (total_kmers + window_size - 1) // window_size
        return False

    #   for w_idx in range(num_windows):
    #       hits = window_counts.get(w_idx, 0)
    #       density = hits / window_size#

    # Check user threshold (e.g. only show if > 0 hits)
    #       if density > threshold:
    # Convert k-mer window to approximate genomic BP for display
    # Start BP = kmer_index
    # End BP = kmer_index + window_size + k_len
    #           bp_start = w_idx * window_size
    #           bp_end = min(seq_len, bp_start + window_size + k_len - 1)

    #           print(f"{seq_id}\t{bp_start}\t{bp_end}\t{hits}\t{density:.2%}")

    #for pair in pairs:
    # Handle special cases like "|:|" for paired-end delimiter
    #   if pair == "|:|": continue

    #    try:
    #        taxid, count = pair.split(':')
    #        count = int(count)
    #    except ValueError:
    #        continue
    #    start_pos = current_bp
    #    end_pos = current_bp + count + self.kraken_config.kmer_size - 1
    #    if taxid == "Synthetic (taxid 1001)":

    #                number_of_synthetic_kmers_in_window = count
    #                if number_of_synthetic_kmers_in_window >= 25:
    #                    return True
    #                else:
    #                    remaining_length_in_window = 200 - end_pos-start_pos % 200

    #            else:
    #                current_bp += count

    #remaining_kmers_in_window =
    # Calculate start and end of this block
    # A block of 'count' k-mers covers 'count' + k - 1 bases
    #  if end_pos-start_pos > window_size_to_check:
    #                if taxid == "Synthetic (taxid 1001)":
    #                    return True
    #else:
    #    continue
    #            else:

    #if end_pos - start_pos > window_size_to_check:

    #            if end_pos-start_pos<=window_size_to_check:
    #                pass
    #            else:

    # If this block matches our target unique TaxID, print BED line
    #if taxid == target_taxid:
    # Print BED format: chrom start end
    #    print(f"{seq_id}\t{start_pos}\t{end_pos}")

    # Advance position
    # We advance by 'count' because the next k-mer starts 1 base after the previous
    #           current_bp += count
    @abc.abstractmethod
    def _rescue_logic(self) -> bool:
        raise NotImplementedError

    def _post_process_kraken(self):
        report = open(self.kraken_db.parent.joinpath("post_processed_kraken_report.txt"), "w")
        with open(self.kraken_db.parent.joinpath("output.txt"), 'r') as kraken_file:
            entries = kraken_file.readlines()
            for entry in tqdm(entries, tot=len(entries)):
                categories = entry.split("\t")
                synthetic_boolean: bool = self._check_windows(categories[-1])
                logging.debug(f"Entry: {categories[1]} found to be synthetic from windowing logic: {synthetic_boolean}")
                if synthetic_boolean:
                    report.write(f"Synthetic\t{categories[1]}")
                else:
                    synthetic_boolean = self._rescue_logic()
                    logging.debug(
                        f"Entry: {categories[1]} found to be synthetic from rescue logic: {synthetic_boolean}")

    def run(self):
        self.run_kraken()
        self._post_process_kraken()


class LRWorkflow(Workflow):

    def __init__(self, fasta_file: str, kraken_db: str, kraken_config: KrakenConfig):
        super().__init__(fasta_file, kraken_db, kraken_config)

    def _rescue_logic(self) -> bool:
        return False

