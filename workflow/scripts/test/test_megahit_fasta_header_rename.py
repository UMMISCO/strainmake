# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import hashlib
import os
import workflow.scripts.megahit_fasta_header_rename as mfh

TEST_ASSEMBLY = (
    "workflow/scripts/test/data/assembly/before_deduplicate_contigs_name.fna"
)
TEST_ASSEMBLY_md5 = "529373287f35bac3d36a4eecb964d6e0"


class TestMegahitFastaHeaderRename(unittest.TestCase):

    def test_replace_spaces_in_fasta(self):

        expected_sequence_names = ["contig.A", "contig.B", "contig.C"]

        # ensuring the test FASTA has not been altered
        with open(TEST_ASSEMBLY, "rb") as f:
            file_md5 = hashlib.md5(f.read()).hexdigest()

        self.assertEqual(
            file_md5,
            TEST_ASSEMBLY_md5,
            f"Test assembly {TEST_ASSEMBLY} has been altered (md5: {file_md5})",
        )

        output_file = (
            "workflow/scripts/test/data/assembly/after_deduplicate_contigs_name.fna"
        )

        # running the function to replace spaces in FASTA headers
        mfh.replace_spaces_in_fasta(TEST_ASSEMBLY, output_file)

        # checking the output file
        # first, extracting the sequence names
        with open(output_file, "r") as f:
            sequence_names = [line.strip()[1:] for line in f if line.startswith(">")]
            # sorting the sequence names to ensure order does not affect the test
            sequence_names.sort()

        # then, checking if the sequence names match the expected ones
        self.assertEqual(
            sequence_names,
            expected_sequence_names,
            f"Expected sequence names {expected_sequence_names}, but got {sequence_names}",
        )

        # clean up the output file after the test
        os.remove(output_file)


if __name__ == "__main__":
    unittest.main()
