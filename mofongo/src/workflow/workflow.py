"""

"""
import abc
import logging
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm

from mofongo.lib.const import KrakenConfig, MetaMDBGConfig, MegahitConfig


class Workflow(abc.ABC):

    def __init__(self, fasta_file: str, output_report_path: str, kraken_db: str, kraken_raw_output:str, kraken_config: KrakenConfig,
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
        try:
            logging.info(f"Running Kraken on {self.fasta_file}")
            #output_file_path: Path = self.kraken_db.parent.joinpath("output.txt")
            report_file_path: Path = self.kraken_db.parent.joinpath("report.txt")
            logging.info(f"Running the command: kraken2 --db {self.kraken_db} --report-minimizer-data "
                         f"--report {report_file_path} --output {self.kraken_output_path} --use-names {self.fasta_file} "
                         f"--threads {self.kraken_config.max_threads}")
            subprocess.check_call(
                ["kraken2", "--db", self.kraken_db, "--report-minimizer-data", "--report", report_file_path,
                 "--output", self.kraken_output_path, "--use-names", self.fasta_file,
                 "--threads", str(self.kraken_config.max_threads)])
        except subprocess.CalledProcessError:
            logging.info(f"CRITICAL ERROR: Failed to run kraken2 {self.fasta_file}.")
            sys.exit(1)
        logging.info("All files processed.")
        self._post_process_kraken()

    @staticmethod
    def _check_windows(kmer_pos_info: str) -> bool:
        """

        :return:
        """
        window_size = 200

        segments = []
        for item in kmer_pos_info.split(' '):
            if item == "|:|": continue  # Skip paired-end delimiters
            try:
                tid, count = item.split(':')
                segments.append((tid, int(count)))
            except ValueError:
                continue

        # Reset for the optimized approach
        current_kmer_idx = 0

        # Dictionary to hold counts: window_index -> count
        # window_counts = {}

        for tid, count in segments:
            seg_start = current_kmer_idx
            seg_end = current_kmer_idx + count

            # If this segment is our target, add counts to appropriate windows
            if tid == "32630":  #"1001":
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

            current_kmer_idx += count
        return False

    @abc.abstractmethod
    def _rescue_logic(self) -> bool:
        raise NotImplementedError

    def _post_process_kraken(self):
        report = open(self.report_output_path, "w")
        with open(self.kraken_output_path , 'r') as kraken_file:
            entries = kraken_file.readlines()
            for entry in tqdm(entries, total=len(entries)):
                categories = entry.split("\t")
                synthetic_boolean: bool = self._check_windows(categories[-1])
                logging.debug(f"Entry: {categories[1]} found to be synthetic from windowing logic: {synthetic_boolean}")
                if synthetic_boolean:
                    report.write(f"Synthetic\t{categories[1]}\n")
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
            cmd = ["metaMDBG", "asm", "--report-minimizer-data","--report",
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

    def __init__(self, fasta_file: str, paired_end: str, output_report_path: str, kraken_db: str,kraken_raw_output: str,
                 kraken_config: KrakenConfig,
                 assembly: bool, assembly_config: MegahitConfig):
        super().__init__(fasta_file, output_report_path, kraken_db, kraken_raw_output,kraken_config, assembly, assembly_config)
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
