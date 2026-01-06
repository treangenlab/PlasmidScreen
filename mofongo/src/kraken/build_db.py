import concurrent.futures
import os
from mofongo.lib.const import custom_taxonomy
from pathlib import Path
from typing import List
import multiprocessing as mp

valid_extensions = ('.fasta', '.fa', '.fna', '.faa')
import logging
import threading


def find_fastas(fasta_files_directory: Path) -> List[Path]:
    logging.info(f"Finding fastas from path {fasta_files_directory}")
    fasta_paths = []
    for root, dirs, files in os.walk(fasta_files_directory):
        for file in files:
            if file.lower().endswith(valid_extensions):
                logging.debug(f"Added path: {Path(root).joinpath(file).as_posix()}")
                fasta_paths.append(Path(root).joinpath(file))
            else:
                logging.debug(f"Ignoring path: {Path(root).joinpath(file)}")
    return fasta_paths


class BuildDB:
    """
    This class is responsible for building the kraken DB and
    """

    def __init__(self, fasta_files_directory: Path, threads_allowed: int):
        self._fasta_paths = find_fastas(fasta_files_directory)
        self.thread_budget = threads_allowed
        self._preprocess_sequences()

    def _preprocess_sequences(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.thread_budget) as executor:
            executor.map()
        pass

    @staticmethod
    def _add_kraken_taxa_to_fasta(fasta_file, output_fasta_location):
        logging.info("Reformating headers for kraken compliance ")
        with open(fasta_file, "r") as f_in, open(output_fasta_location, "w") as f_out:
            for line in f_in:
                if line.startswith(">"):
                    # Parse the ID (everything after > up to the first space)
                    original_header = line.strip()
                    seq_id = original_header[1:].split()[0]

                    # if seq_id in mapping_dict:
                    #    tax_id = mapping_dict[seq_id]
                    # Create new header: >OriginalID|kraken:taxid|TAXID OriginalDescription
                    # We strip the '>' first, add the tag, then print
                    if "IMGPR" in original_header:
                        new_header = f">{seq_id}|kraken:taxid|1012 {original_header[len(seq_id) + 1:]}\n"
                        logging.debug(f"Added kraken:taxid|1012 to {original_header}")
                    elif "Escherichia coli" in original_header:
                        new_header = f">{seq_id}|kraken:taxid|562 {original_header[len(seq_id) + 1:]}\n"
                        logging.debug(f"Added kraken:taxid|562 to {original_header}")
                    else:
                        new_header = f">{seq_id}|kraken:taxid|1001 {original_header[len(seq_id) + 1:]}\n"
                        logging.debug(f"Added kraken:taxid|1001 to {original_header}")
                    f_out.write(new_header)
                else:
                    if line != "\n":
                        # Write sequence lines unchanged
                        f_out.write(line)
