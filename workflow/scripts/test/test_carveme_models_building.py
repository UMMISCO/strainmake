# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import io
import sys
from workflow.scripts import carveme_models_building as cmb

BINS_DIR = "workflow/scripts/test/data/bins"


class TestCarvemeModelsBuilding(unittest.TestCase):
    def test_list_mag(self):

        expected_files = [
            "workflow/scripts/test/data/bins/bin_34.fa",
            "workflow/scripts/test/data/bins/bin_2.fa",
        ]

        list_bins = cmb.list_mag(BINS_DIR, ".fa", verbose=False)

        self.assertEqual(sorted(list_bins), sorted(expected_files))

    def test_list_mag_no_results(self):
        """
        In this test, we check that the function returns an empty list when no files match the criteria.
        """

        expected_files = []

        list_bins = cmb.list_mag(BINS_DIR, ".fna", verbose=False)

        self.assertEqual(list_bins, expected_files)

    def test_carveme_carve(self):

        expected_output = (
            "carve --solver gurobi --dna --output workflow/scripts/test/data/carveme_models/bin_34.fa.xml workflow/scripts/test/data/bins/bin_34.fa\n"
            "carve --solver gurobi --dna --output workflow/scripts/test/data/carveme_models/bin_2.fa.xml workflow/scripts/test/data/bins/bin_2.fa\n"
        )

        list_bins = cmb.list_mag(BINS_DIR, ".fa", verbose=False)
        captured_output = io.StringIO()
        sys_stdout = sys.stdout
        sys.stdout = captured_output

        try:
            cmb.carveme_carve(
                mag_list=list_bins,
                out_dir="workflow/scripts/test/data/carveme_models",
                cpu=1,
                verbose=False,
                dryrun=True,
            )
        finally:
            sys.stdout = sys_stdout

        # asserting that the captured output matches the expected output, i.e., the commands that would be run
        self.assertEqual(captured_output.getvalue(), expected_output)

    def test_carveme_merge_community(self):

        expected_command = "merge_community workflow/scripts/test/data/carveme_models/* -o workflow/scripts/test/data/carveme_models/community.xml\n"

        community_dir = "workflow/scripts/test/data/carveme_models"
        captured_output = io.StringIO()
        sys_stdout = sys.stdout
        sys.stdout = captured_output

        try:
            cmb.carveme_merge_community(
                communities=community_dir,
                out_dir="workflow/scripts/test/data/carveme_models",
                verbose=False,
                dryrun=True,
            )
        finally:
            sys.stdout = sys_stdout

        # asserting that the captured output matches the expected output, i.e., the command that would be run
        self.assertEqual(captured_output.getvalue(), expected_command)


if __name__ == "__main__":
    unittest.main()
