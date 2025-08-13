# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
from workflow.scripts import generate_bakta_commands as gbc

GTDB_TK_SUMMARY_PATH = "workflow/scripts/test/data/gtdbtk.bac120.summary.tsv"
BAKTA_DB = "fake/path/to/bakta_db"
BINS_DIR = "data/genomes"
OUTPUT_DIR_ANNOTATIONS = "data/bakta_annotations"
# it is a folder that already contains Bakta results (in JSON)
OUTPUT_DIR_ANNOTATIONS_WITH_RESULTS = "workflow/scripts/test/data/bakta/results"
EXTENSION = ".fna"
THREADS = 4


class TestGenerateBaktaCommands(unittest.TestCase):
    def test_generate_bakta_commands_annot(self):

        expected_commands = [
            "bakta --genus Streptococcus --species pneumoniae --skip-plot --output data/bakta_annotations/bin_2 --prefix bin_2 --threads 4 --verbose --db fake/path/to/bakta_db data/genomes/bin_2.fna",
            "bakta --genus Escherichia --species coli --skip-plot --output data/bakta_annotations/bin_34 --prefix bin_34 --threads 4 --verbose --db fake/path/to/bakta_db data/genomes/bin_34.fna",
            "bakta --genus Klebsiella --species pneumoniae --skip-plot --output data/bakta_annotations/bin_36 --prefix bin_36 --threads 4 --verbose --db fake/path/to/bakta_db data/genomes/bin_36.fna",
            "bakta --genus Staphylococcus --species aureus --skip-plot --output data/bakta_annotations/bin_4 --prefix bin_4 --threads 4 --verbose --db fake/path/to/bakta_db data/genomes/bin_4.fna",
            "bakta --genus Enterococcus_B --species faecium --skip-plot --output data/bakta_annotations/bin_6 --prefix bin_6 --threads 4 --verbose --db fake/path/to/bakta_db data/genomes/bin_6.fna",
            "bakta --genus Cutibacterium --species acnes --skip-plot --output data/bakta_annotations/bin_8_1 --prefix bin_8_1 --threads 4 --verbose --db fake/path/to/bakta_db data/genomes/bin_8_1.fna",
        ]

        bakta_commands = gbc.generate_bakta_commands_annot(
            gtdb_summary_path=GTDB_TK_SUMMARY_PATH,
            genomes_dir=BINS_DIR,
            out_dir=OUTPUT_DIR_ANNOTATIONS,
            bakta_database=BAKTA_DB,
            threads=THREADS,
            extension=EXTENSION,
        )

        self.assertEqual(len(bakta_commands), len(expected_commands))
        for cmd in expected_commands:
            self.assertIn(cmd, bakta_commands)

    def test_generate_bakta_commands_plot(self):

        expected_commands = [
            "bakta_plot --output data/bakta_annotations/bin_2 --type cog --verbose workflow/scripts/test/data/bakta/results/bin_2/bin_2.json",
            "bakta_plot --output data/bakta_annotations/bin_34 --type cog --verbose workflow/scripts/test/data/bakta/results/bin_34/bin_34.json",
            "bakta_plot --output data/bakta_annotations/bin_36 --type cog --verbose workflow/scripts/test/data/bakta/results/bin_36/bin_36.json",
            "bakta_plot --output data/bakta_annotations/bin_4 --type cog --verbose workflow/scripts/test/data/bakta/results/bin_4/bin_4.json",
            "bakta_plot --output data/bakta_annotations/bin_6 --type cog --verbose workflow/scripts/test/data/bakta/results/bin_6/bin_6.json",
            "bakta_plot --output data/bakta_annotations/bin_8_1 --type cog --verbose workflow/scripts/test/data/bakta/results/bin_8_1/bin_8_1.json",
        ]

        bakta_commands = gbc.generate_bakta_commands_plot(
            gtdb_summary_path=GTDB_TK_SUMMARY_PATH,
            bakta_annot_dir=OUTPUT_DIR_ANNOTATIONS_WITH_RESULTS,
            out_dir=OUTPUT_DIR_ANNOTATIONS,
        )

        self.assertEqual(len(bakta_commands), len(expected_commands))
        for cmd in expected_commands:
            self.assertIn(cmd, bakta_commands)

        print(bakta_commands)


if __name__ == "__main__":
    unittest.main()
