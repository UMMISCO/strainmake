# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import os
from workflow.scripts import deduplicate_contigs_name as dcn

FASTA_DIR = "workflow/scripts/test/data/assembly"

class TestDeduplicateContigsName(unittest.TestCase):
    def test_list_fasta(self):

        # copying all *.fna files in FASTA_DIR to the same directory with .fa extension
        for filename in os.listdir(FASTA_DIR):
            if filename.endswith(".fna"):
                src = os.path.join(FASTA_DIR, filename)
                dst = os.path.join(FASTA_DIR, os.path.splitext(filename)[0] + ".fa")

                with open(src, "rb") as fsrc, open(dst, "wb") as fdst:
                    fdst.write(fsrc.read())

        expected_files = {
            os.path.join(FASTA_DIR, "before_deduplicate_contigs_name.fa"),
            os.path.join(FASTA_DIR, "before_deduplicate_contigs_name2.fa")
        }
        list_fasta_files = dcn.list_fasta(FASTA_DIR)

        # asserting that the function idenitified the correct FASTA files
        self.assertSetEqual(set(list_fasta_files), expected_files,
                            f"Expected {expected_files} but got {set(list_fasta_files)}")
        
    def test_rename_fasta_headers(self):

        list_fasta_files = dcn.list_fasta(FASTA_DIR)
        dcn.rename_fasta_headers(list_fasta_files)

        # checking that the headers in the output file are as expected after renaming
        output_file_1 = os.path.join(FASTA_DIR, "before_deduplicate_contigs_name.fa")
        expected_headers_fasta_1 = [
            ">before_deduplicate_contigs_name_contig A",
            ">before_deduplicate_contigs_name_contig B",
            ">before_deduplicate_contigs_name_contig C"
        ]

        with open(output_file_1, "r") as f:
            headers = [line.strip() for line in f if line.startswith(">")]

        self.assertListEqual(headers, expected_headers_fasta_1, f"Expected headers {expected_headers_fasta_1} but got {headers}")

        # checking the second fasta file
        output_file_2 = os.path.join(FASTA_DIR, "before_deduplicate_contigs_name2.fa")
        expected_headers_fasta_2 = [
            ">before_deduplicate_contigs_name2_contig A",
            ">before_deduplicate_contigs_name2_contig B"
        ]

        with open(output_file_2, "r") as f:
            headers = [line.strip() for line in f if line.startswith(">")]

        self.assertListEqual(headers, expected_headers_fasta_2, f"Expected headers {expected_headers_fasta_2} but got {headers}")

        # cleaning up the test files
        for filename in os.listdir(FASTA_DIR):
            if filename.endswith(".fa"):
                os.remove(os.path.join(FASTA_DIR, filename))

if __name__ == "__main__":
    unittest.main()