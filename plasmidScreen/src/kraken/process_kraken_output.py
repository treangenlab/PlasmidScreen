from pathlib import Path


class FilteringKrakenOutputReport:

    def check_windows(self, entry):
        kmer_list = entry.split(" ")
        for kmer_index in range(0):
            pass

    def process_file(self, kraken_output_file: Path):
        report = open(kraken_output_file.parent.joinpath("post_processed_kraken_report.txt"), "w")
        with open(kraken_output_file, 'r') as kraken_file:
            entries = kraken_file.readlines()
            for entry in entries:
                categories = entry.split("\t")
                self.check_windows(entry[-1])
                #if categories[2] == "Synthetic (taxid 1001)":

                report.write(f"Synthetic\t{categories[1]}")

#"C       0b1a2182-7409-4231-a012-c9ca3a8defee    Ecoli (taxid 562)       4544    0:298 1012:1 0:13 1012:1 0:199 1012:1 0:981 1012:1 0:486 1012:1 0:186 562:1 0:137 562:31 0:32 562:30 0:27 562:9 0:390 1012:5 0:698 1012:1 0:995"
