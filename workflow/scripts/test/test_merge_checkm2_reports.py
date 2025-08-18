# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import hashlib
import pandas as pd
import workflow.scripts.merge_checkm2_reports as mcr

DUMMY_CHECKM2_REPORTS = [
    "workflow/scripts/test/data/checkm2_reports/hybridspades_semibin2_dummy.tsv",
    "workflow/scripts/test/data/checkm2_reports/hybridspades_vamb_dummy.tsv",
    "workflow/scripts/test/data/checkm2_reports/megahit_metabat2_dummy.tsv"
]
EXPECTED_MERGING_RESULT = "workflow/scripts/test/data/checkm2_reports/expected_results.tsv"
EXPECTED_MERGING_RESULT_md5 = "aeeeced9522ae3bcfa436f267e0d7d77" 

class TestMergeCheckM2Reports(unittest.TestCase):

    def test_identify_assembly_program(self):

        # test cases for different assembly programs
        test_cases = {
            "path/to/metaspades_report.txt": "metaspades",
            "path/to/megahit_report.txt": "megahit",
            "path/to/metaflye_report.txt": "metaflye",
            "path/to/hybridspades_report.txt": "hybridspades",
            "path/to/hylight_report.txt": "hylight",
            "path/to/hylight.txt": "hylight",
            "path/to/unknown_program_report.txt": None
        }

        # iterating through test cases, asserting the expected function output
        for report_path, expected_program in test_cases.items():
            with self.subTest(report_path=report_path):
                result = mcr.identify_assembly_program(report_path)
                self.assertEqual(result, expected_program,
                                 f"Expected {expected_program} for {report_path}, but got {result}")
                
    def test_identify_binning_program(self):

        # test cases for different binning programs
        test_cases = {
            "path/to/semibin2_report.txt": "sembin2",
            "path/to/metabat2_report.txt": "metabat2",
            "path/to/vamb_report.txt": "vamb",
            "path/to/vamb/report.txt": "vamb",
            "path/to/unknown_binning_program_report.txt": None
        }

        # iterating through test cases, asserting the expected function output
        for report_path, expected_program in test_cases.items():
            with self.subTest(report_path=report_path):
                result = mcr.identify_binning_program(report_path)
                self.assertEqual(result, expected_program,
                                 f"Expected {expected_program} for {report_path}, but got {result}")
                
    def test_merge_reports(self):

        expected_dataframe = pd.read_csv(EXPECTED_MERGING_RESULT, sep="\t")

        # ensuring the reference file has not been altered
        with open(EXPECTED_MERGING_RESULT, "rb") as f:
            file_md5 = hashlib.md5(f.read()).hexdigest()

        self.assertEqual(file_md5, EXPECTED_MERGING_RESULT_md5, 
                 f"Reference file {EXPECTED_MERGING_RESULT} has been altered (md5: {file_md5})")

        # merging dummy reports
        merged_df = mcr.merge_reports(DUMMY_CHECKM2_REPORTS)

        # checking if the merged dataframe matches the expected dataframe
        pd.testing.assert_frame_equal(merged_df, expected_dataframe)    

if __name__ == "__main__":
    unittest.main()