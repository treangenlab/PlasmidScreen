from pathlib import Path

from plasmidScreen.api import ScreenResult, build_codon_database, run_screen

# Have the path ready of the fastq file of interest

fastq = Path("data/example_univec_sequence.fasta")
example_2 = Path("/dodo/fmq1/pKT0174.fasta")
if __name__ == "__main__":

    # build_codon_database(output_dir="/home/Users/fmq1/tmp/")
    my_screen_result: ScreenResult = run_screen(
        fasta_file=example_2,
        seq_type=1,
        engineered_kmer_threshold=20,
        kraken_db="/dodo/fmq1/PlasmidScreenDBNoAddgene",
        diamond_db="/dodo/dbs/uniref_march_2025_with_tax.dmnd",
        threads=60,
        quiet_mode=True,
        run_codon_usage=True,
        memory_mode=0,
        background_rate=0.0003,
    )
    print(my_screen_result.engineered_read_ids)  # overall
    print(my_screen_result.engineered_scan.engineered_read_ids)  # k-mer
    print(my_screen_result.codon_cai_engineered_read_ids)  # CAI
    print(my_screen_result.natural_read_ids_overall)
