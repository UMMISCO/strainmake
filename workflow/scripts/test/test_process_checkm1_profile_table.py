# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import os
import pandas as pd
from workflow.scripts import process_checkm1_profile_table as pcpt

BEFORE = "workflow/scripts/test/data/profile.tsv"
AFTER = "workflow/scripts/test/data/profile.processed.tsv"

AFTER_DF = pd.read_csv(AFTER, sep="\t")

class TestProcessCheckm1ProfileTable(unittest.TestCase):
    def test_rename_columns(self):

        output_file = "workflow/scripts/test/data/profile.processed.test.tsv"
        pcpt.rename_columns(BEFORE, output_file)

        # read the output file and compare it with the expected DataFrame
        output_df = pd.read_csv(output_file, sep="\t")
        
        # check if DataFrame are the same
        pd.testing.assert_frame_equal(output_df, AFTER_DF)

        os.remove(output_file)

    def test_cli(self):
        """
        Here we specifically test the CLI interface of the script.
        """

        output_file = "workflow/scripts/test/data/profile.processed.cli.test.tsv"
        os.system(f"python3 workflow/scripts/process_checkm1_profile_table.py --input_table {BEFORE} {output_file}")

        # read the output file and compare it with the expected DataFrame
        output_df = pd.read_csv(output_file, sep="\t")

        # check if DataFrame are the same
        pd.testing.assert_frame_equal(output_df, AFTER_DF)

        os.remove(output_file)

if __name__ == "__main__":
    unittest.main()