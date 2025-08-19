# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
from workflow.scripts import filter_dereplicated_bins_by_quality as fdbq
import os

CHECKM2_REPORT = "workflow/scripts/test/data/quality_report.tsv"
BINS_DIRECTORY = "workflow/scripts/test/data/bins"
OUTPUT_DIR = "workflow/scripts/test/data/filtered_bins"


class TestFilterDereplicatedBinsByQuality(unittest.TestCase):
    def test_filter_bins(self):

        os.makedirs(OUTPUT_DIR, exist_ok=True)

        # defining test parameters
        min_completeness = 50.0
        max_contamination = 10.0

        # calling the function to test
        fdbq.filter_bins(
            CHECKM2_REPORT,
            BINS_DIRECTORY,
            min_completeness,
            max_contamination,
            OUTPUT_DIR,
        )

        # asserting the expected outcome: only bin_2.fa and bin_34.fa should be present in OUTPUT_DIR
        expected_files = {"bin_2.fa", "bin_34.fa"}
        bin_files = {
            f
            for f in os.listdir(OUTPUT_DIR)
            if os.path.isfile(os.path.join(OUTPUT_DIR, f))
        }
        self.assertSetEqual(
            bin_files,
            expected_files,
            f"Expected only {expected_files} in the output directory",
        )

        # Remove all files in OUTPUT_DIR
        for f in os.listdir(OUTPUT_DIR):
            file_path = os.path.join(OUTPUT_DIR, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
        os.rmdir(OUTPUT_DIR)

    def test_filter_bins_no_bins(self):
        """
        Test the case where no bins pass the filtration criteria.
        """

        os.makedirs(OUTPUT_DIR, exist_ok=True)

        # defining test parameters
        # given CHECKM2_REPORT, with such parameters, no bins should pass the filtration
        min_completeness = 99.0
        max_contamination = 1.0

        # calling the function to test
        fdbq.filter_bins(
            CHECKM2_REPORT,
            BINS_DIRECTORY,
            min_completeness,
            max_contamination,
            OUTPUT_DIR,
        )

        # asserting the expected outcome: no bins should be present in OUTPUT_DIR
        expected_files = set()
        bin_files = {
            f
            for f in os.listdir(OUTPUT_DIR)
            if os.path.isfile(os.path.join(OUTPUT_DIR, f))
        }
        self.assertSetEqual(
            bin_files,
            expected_files,
            f"Expected only {expected_files} in the output directory",
        )

        # Remove all files in OUTPUT_DIR
        for f in os.listdir(OUTPUT_DIR):
            file_path = os.path.join(OUTPUT_DIR, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
        os.rmdir(OUTPUT_DIR)


if __name__ == "__main__":
    unittest.main()
