import os, sys
import subprocess
from pathlib import Path
import logging

from Bio import SeqIO
import tempfile

# work_folder='/mnt/c/Users/miken/Rice/Research/dtra_ww/addgene/moclo_tool_kit_genbank_files_2_1'
work_folder = '/mnt/c/Users/miken/Rice/Research/dtra_ww/addgene/crispr_bacterial'
gb_folder = os.path.join(work_folder, 'addgene_genbank_files')
blast_out_path = os.path.join(work_folder, 'crispr_bacterial_blastout.tsv')
plot_output_folder = os.path.join(work_folder, 'pileup_plots')


def run_build_blast_db(fasta_file: Path):
    try:
        subprocess.check_call([
            "makeblastdb",
            "-in",
            fasta_file.as_posix(),
            "-title",
            fasta_file.with_suffix(".mini"),
            "-out",
            fasta_file.with_suffix(".mini"),
            "-dbtype",
            "nucl"
        ])
        logging.info(f"Success!")
    except subprocess.CalledProcessError:
        print("Error: The build step failed.")
        sys.exit(1)


def run_blastn(self, fasta_file: Path, blastdb_path):
    outfile = fasta_file.parent.joinpath("blastout")
    try:
        subprocess.check_call([
            "blastn",
            "-query",
            fasta_file.as_posix(),
            "-db",
            blastdb_path,
            "-out",
            outfile,
            "-outfmt",
            "6",
            "-num_threads",
            self.threads
        ])
        logging.info(f"Success!")
    except subprocess.CalledProcessError:
        print("Error: The build step failed.")
        sys.exit(1)



class BlastVisualization:

    def __init__(self, hits_file_to_parse, fasta_file_of_interest, blast_db, threads):
        """

        :param hits_file_to_parse: File generated from PlasmidScreen after parsing kraken2 taxa hits.
        """
        self.query = self.parse_file(hits_file_to_parse, fasta_file_of_interest)
        self.blast_db = blast_db
        self.threads = threads

    def run(self):
        blast_output_file = self.run_blastn()
        blast_out_parsed_path = blast_output_file.parent.joinpath("blastout_parsed.tsv")
        self.parse_blast(blast_output_file, blast_out_parsed_path)

    def run_blastn(self) -> Path:
        outfile = self.query.parent.joinpath("blastout.tsv")
        try:
            subprocess.check_call([
                "blastn",
                "-query",
                self.query.name,
                "-db",
                self.blast_db,
                "-out",
                outfile,
                "-outfmt",
                "6",
                "-num_threads",
                self.threads
            ])
            logging.info(f"Successfully Ran Blast on {self.query.name} against {self.blast_db}!")
        except subprocess.CalledProcessError:
            logging.error("Error: The build step failed.")
            sys.exit(1)
        return outfile


    @staticmethod
    def blast_hits_visualization(parsed_blast_hits, library_of_genbank_files:Path):
        parsed_blast_hits_file = open(parsed_blast_hits, "r")
        for entry in parsed_blast_hits_file:
            entry = entry.strip().split('\t')
            plasmid_id = entry[0]
            start_pos = entry[4]
            end_pos = entry[5]
            record = SeqIO.read(library_of_genbank_files.joinpath(plasmid_id), "genbank")
            for feature in record.features:
                # 10 start # 20 end
                # 8          25
                overlap_start = max(int(start_pos), feature.locations.start)
                overlap_end = min(int(end_pos), feature.locations.end)

                overlap_len = overlap_end - overlap_start

                # If overlap is positive, we have a hit
                if overlap_len > 0:
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                    product = feature.qualifiers.get("product", ["Unknown"])[0]

                    # Calculate % coverage of the gene
                    gene_len = feature.locations.end - feature.locations.start
                    coverage = (overlap_len / gene_len) * 100

                    print(f"Match found in gene: {gene_name}")
                    print(f"  Product: {product}")
                    print(f"  Overlap: {overlap_len} bp ({coverage:.1f}% of gene covered)")
                    print("-" * 30)
                    #if feature.locations.start >= start_pos and feature.locations.start <=end_pos
                #if feature.type == "CDS":
                #    gene_name = feature.qualifiers.get("gene", ["N/A"])[0]
                #    product = feature.qualifiers.get("product", ["N/A"])[0]
                #    start = feature.location.start
                #    end = feature.location.end            #plasmid_id_genbank_file = open(library_of_genbank_files.joinpath(plasmid_id), 'r')


    @staticmethod
    def parse_blast(input_file, output_file, min_cov=0.8, min_bitscore=100):
        """
        Parses a BLAST tab-separated file (-outfmt "6 std qlen slen")
        and filters based on coverage and/or bitscore.
        """

        # Column indices based on -outfmt "6 std qlen slen"
        # Python uses 0-based indexing
        COL_ALIGN_LEN = 3  # length (4th column)
        COL_BITSCORE = 11  # bitscore (12th column)
        COL_QLEN = 12  # qlen (13th column) - Custom added column

        count_kept = 0
        count_total = 0

        try:
            with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    # Skip comment lines if any
                    if line.startswith("#"):
                        fout.write(line)
                        continue
                    parts = line.strip().split('\t')

                    # Ensure the line has enough columns
                    if len(parts) < 13:
                        sys.stderr.write(f"Warning: Line {count_total + 1} malformed or missing 'qlen'. Skipping.\n")
                        continue

                    count_total += 1

                    # Parse numerical values
                    align_len = float(parts[COL_ALIGN_LEN])
                    bitscore = float(parts[COL_BITSCORE])
                    qlen = float(parts[COL_QLEN])

                    # Calculate Coverage
                    # Note: align_len includes gaps.
                    # For strict coverage (matches only), use (qend - qstart + 1) / qlen
                    coverage = (align_len / qlen) * 100.0

                    # Filtering Logic
                    pass_cov = True
                    pass_bits = True

                    if min_cov is not None:
                        if coverage < min_cov:
                            pass_cov = False

                    if min_bitscore is not None:
                        if bitscore < min_bitscore:
                            pass_bits = False

                        # Write if it passes all checks
                    if pass_cov and pass_bits:
                        fout.write(line)
                        count_kept += 1

            logging.info(f"Done. Processed {count_total} hits.")
            logging.info(f"Retained {count_kept} hits matching criteria.")
            logging.info(f"Output saved to: {output_file}")

        except FileNotFoundError:
            print(f"Error: The file '{input_file}' was not found.")
        except Exception as e:
            print(f"An error occurred: {e}")

    @staticmethod
    def parse_file(hits_file, fasta_file):
        records = []

        # Use a context manager for automatic cleanup

        tmp_file = tempfile.NamedTemporaryFile(mode='w+t', delete=False, suffix='.txt')

        with open(hits_file, "r") as f_in:
            synthetic_hits = [line.split(":") for line in f_in if "Synthetic" in line]
        for record in SeqIO.parse(fasta_file, format="fasta"):
            if record.id in synthetic_hits:
                records.append(record)
        SeqIO.write(records, tmp_file, format="fasta")
        return tmp_file
