import concurrent.futures
import os
import shutil
import tempfile

from mofongo.lib.const import custom_taxonomy, KrakenConfig
from pathlib import Path
from typing import List
import multiprocessing as mp
import subprocess
import sys
import glob
import json
from tqdm import tqdm
valid_extensions = ('.fasta', '.fa', '.fna', '.faa')
import logging
import threading


def find_fastas(fasta_files_directory: Path) -> List[Path]:
    logging.info(f"Finding fastas from path {fasta_files_directory}")
    fasta_paths = []
    for root, dirs, files in os.walk(fasta_files_directory):
        for file in files:
            if file.lower().endswith(valid_extensions):
                logging.debug(f"Added path: {Path(root).joinpath(file)}")
                fasta_paths.append(Path(root).joinpath(file))
            else:
                logging.debug(f"Ignoring path: {Path(root).joinpath(file)}")
    return fasta_paths


def write_nodes_names(taxonomy_list, output_dir="."):
    nodes_path = os.path.join(output_dir, "nodes.dmp")
    names_path = os.path.join(output_dir, "names.dmp")

    with open(nodes_path, "w") as f_nodes, open(names_path, "w") as f_names:
        for tax_id, parent_id, rank, name in taxonomy_list:
            # Write to nodes.dmp
            # Columns: tax_id | parent_tax_id | rank | embl_code | division_id | ...
            # We pad with empty fields to match standard NCBI format breadth approximately
            f_nodes.write(
                f"{tax_id}\t|\t{parent_id}\t|\t{rank}\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")

            # Write to names.dmp
            # Columns: tax_id | name_txt | unique_name | name_class
            f_names.write(f"{tax_id}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    logging.info(f"Generated {nodes_path} and {names_path}")


# --- Configuration ---
output_file = "combined_all.fasta"


# --- 1. The Worker Function ---
class BuildDB:
    """
    This class is responsible for building the kraken DB and
    """

    def __init__(self, working_directory: Path, kraken_db_name: str, natural_sequences: str, engineered_sequences: str,
                 additional_sequences: str, kraken_config: KrakenConfig):
        self.working_directory = working_directory
        self.kraken_db = self.working_directory.joinpath(kraken_db_name)
        self.kraken_db.mkdir(exist_ok=True)
        self.kraken_db.joinpath("taxonomy").mkdir(exist_ok=True)
        natural_sequences_with_kraken = self.add_kraken_taxa_to_fasta(natural_sequences)
        engineered_sequences_with_kraken = self.add_kraken_taxa_to_fasta(engineered_sequences)
        additional_sequences_with_kraken = self.add_kraken_taxa_to_fasta(additional_sequences)
        self.add_sequence_to_kraken_db(natural_sequences_with_kraken,self.kraken_db)
        self.add_sequence_to_kraken_db(engineered_sequences_with_kraken, self.kraken_db)
        self.add_sequence_to_kraken_db(additional_sequences_with_kraken, self.kraken_db)
        self.kraken_config = kraken_config
        write_nodes_names(custom_taxonomy, str(self.kraken_db.joinpath("taxonomy")))

    @staticmethod
    def add_kraken_taxa_to_fasta(fasta_file):
        logging.info("Reformating headers for kraken compliance ")
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_out, open(fasta_file, "r") as f_in:
            logging.debug(f"Modifying file {f_in}")
            lines = f_in.readlines()
            for line in tqdm(lines):
                if line.startswith(">"):
                    # Parse the ID (everything after > up to the first space)
                    original_header = line.strip()
                    # Remove '>' for parsing, get the first word as ID
                    seq_id = original_header[1:].split()[0]

                    # Extract the description part (everything after the ID)
                    # We add +1 to account for the space after the ID
                    description = original_header[len(seq_id) + 1:].strip()

                    # Logic to Determine new header
                    if "IMGPR" in original_header:
                        tax_id = "1012"
                    elif "Escherichia coli" in original_header:
                        tax_id = "562"
                    else:
                        tax_id = "1001"

                    # Construct new header
                    new_header = f">{seq_id}|kraken:taxid|{tax_id} {description}\n"
                    logging.debug(f"Modified header: {new_header.strip()}")

                    temp_out.write(new_header)

                else:
                    # Handle sequence lines
                    if line.strip():  # Ensure we don't write purely empty lines if undesired
                        temp_out.write(line)

        # 2. At this point, the temp file is closed and saved on disk.
        #    Now, move the temp file OVER the original file.
        new_file_path=fasta_file.split('.')[0]+"_modified_taxa.fasta"

        shutil.move(temp_out.name, new_file_path)
        logging.info(f"Successfully wrote {new_file_path}")
        return new_file_path
    #def _preprocess_sequences(self, fasta_file):
    #    #with concurrent.futures.ThreadPoolExecutor(max_workers=self.thread_budget) as executor:
    #        executor.map(BuildDB.add_kraken_taxa_to_fasta, fasta_file)
    #    with open(self.working_directory.parent.joinpath("concatenated_sequences.fasta"), "w") as outfile:
    #        for filename in fasta_paths:
    #            logging.info(f"Processing {filename}...")
    #        try:
    #            with open(filename, "r") as infile:
    #                # Read the file line by line (memory efficient)
    #                for line in infile:
    #                    outfile.write(line)
    #        except FileNotFoundError:
    #            logging.debug(f"Warning: Could not find {filename}, skipping.")
    #    logging.info(f"Done! Created {output_file}")

    @staticmethod
    def add_sequence_to_kraken_db(sequence, db_name):
        """
        Wrapper to add sequences and build a Kraken2 DB.

        :param db_name: Path/Name of the Kraken2 database directory
        :param input_pattern: Wildcard string for input files (e.g., "genomes/*.fasta")
        :param threads: Number of CPU threads to use for the build step
        """

        # Get list of files based on pattern
        # files = glob.glob(input_pattern)

        try:
            logging.info(f"Adding {sequence} to kraken db {db_name}")
            # subprocess.check_call raises an error if the command fails
            subprocess.check_call([
                "kraken2-build",
                "--add-to-library", sequence,
                "--db", db_name
            ])
        except subprocess.CalledProcessError:
            logging.error(f"CRITICAL ERROR: Failed to add {sequence}. Stopping.")
            sys.exit(1)

        logging.info("File added.")

    def build_kraken_db(self):
        # 2. Build database
        try:
            subprocess.check_call([
                "kraken2-build",
                "--build",
                "--db", self.kraken_db,
                "--threads", str(self.kraken_config.max_threads)
            ])
            logging.info(f"Success! Database '{self.kraken_db}' is ready.")
        except subprocess.CalledProcessError:
            logging.info("Error: The build step failed.")
            sys.exit(1)
