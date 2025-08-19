# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import sys
import os
from io import StringIO
from workflow.scripts import genes_prediction as gp

INPUT_GENOME = "workflow/scripts/test/data/prodigal/contigs_example.fa"
OUTPUT_GENOME = "workflow/scripts/test/data/prodigal/contigs_example.genes.fna"


class TestGenesPrediction(unittest.TestCase):
    def test_run_prodigal_dry_run(self):

        expected_output = "Dry run:  prodigal -i workflow/scripts/test/data/prodigal/contigs_example.fa -d workflow/scripts/test/data/prodigal/contigs_example.genes.fna\n"

        captured_output = StringIO()
        sys_stdout = sys.stdout
        sys.stdout = captured_output

        try:
            gp.run_prodigal(INPUT_GENOME, OUTPUT_GENOME, dry_run=True)
        finally:
            sys.stdout = sys_stdout

        self.assertEqual(expected_output, captured_output.getvalue())

    def test_run_prodigal(self):
        """
        This test will actually run Prodigal on a small example file.
        Make sure Prodigal is installed and available in the PATH.
        """
        gp.run_prodigal(INPUT_GENOME, OUTPUT_GENOME)

        # Check if the output file was created
        self.assertTrue(os.path.exists(OUTPUT_GENOME))

        # Clean up the output file after the test
        os.remove(OUTPUT_GENOME)


if __name__ == "__main__":
    unittest.main()
