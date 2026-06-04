from plasmidScreen.api import run_screen, ScreenResult
from pathlib import Path

# Have the path ready of the fastq file of interest

fastq = Path("data/example_univec_sequence.fasta")

if __name__ == '__main__':
    my_screen_result: ScreenResult = run_screen(fasta_file=fastq, kraken_db="/dodo/dbs/PlasmidScreenMinimizers",
                                                diamond_db="/dodo/dbs/uniref_march_2025_with_tax.dmnd",threads=60)
    # Here is the found engineered reads
    print(my_screen_result.engineered_read_ids)
    # Here is the engineered reads based on kmer scanning
    print(my_screen_result.engineered_scan.engineered_read_ids)
    # Here is the engineered reads from codon optimization detection
    print(my_screen_result.codon_adaptation.engineered_read_ids)