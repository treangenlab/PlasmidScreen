import os, sys
import subprocess
from pathlib import Path
import logging

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


def parse_blast(input_file, output_file, min_cov, min_bitscore):
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


def blast_hits_visualization(parsed_blast_hits):
    parsed_blast_hits_file = open(parsed_blast_hits, "r")
    for entry in parsed_blast_hits_file:
        entry = entry.strip().split('\t')
        plasmid_id = entry[0]
        plasmid_id_genbank_file = open(plasmid_id.join(".gbk"), 'r')
