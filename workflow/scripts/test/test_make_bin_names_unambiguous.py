# run from root of the repository
# python3 -m unittest discover -s workflow/scripts/test/

import unittest
import hashlib
import pandas as pd
import os
import workflow.scripts.make_bin_names_unambiguous as mbnu

LIST_BINS = "workflow/scripts/test/data/bins_list.txt"
LIST_BINS_WITH_ONE_DUPLICATE_ONE_GZIPPED = (
    "workflow/scripts/test/data/bins_list_with_duplicate_names.txt"
)
UNDUPLICATED_BINS_TABLE_EXPECTED = (
    "workflow/scripts/test/data/bins_list_unduplicated_expected.tsv"
)
UNDUPLICATED_BINS_TABLE_EXPECTED_md5 = "57bac8d9d81220ff9cedf1a851175d78"


class TestMakeBinNamesUnambiguous(unittest.TestCase):

    def test_get_bin_filename(self):

        expected_filenames = ["bin_2.fa", "bin_34.fa"]

        bins_df = pd.read_csv(LIST_BINS, sep="\t", names=["path"])

        bins_filename = mbnu.get_bin_filename(bins_df)
        bins_filename_list = bins_filename["filename"].tolist()

        # asserting that the filenames match the expected ones
        self.assertEqual(
            bins_filename_list,
            expected_filenames,
            f"Expected filenames {expected_filenames}, but got {bins_filename_list}",
        )

    def test_make_unduplicated_filenames(self):

        # ensuring the reference file has not been altered
        with open(UNDUPLICATED_BINS_TABLE_EXPECTED, "rb") as f:
            file_md5 = hashlib.md5(f.read()).hexdigest()

        self.assertEqual(
            file_md5,
            UNDUPLICATED_BINS_TABLE_EXPECTED_md5,
            f"Expected md5 {UNDUPLICATED_BINS_TABLE_EXPECTED_md5}, but got {file_md5}",
        )

        expected_table = pd.read_csv(UNDUPLICATED_BINS_TABLE_EXPECTED, sep="\t")
        bins_df = pd.read_csv(LIST_BINS, sep="\t", names=["path"])

        # adding a duplicate for testing
        bins_df = pd.concat(
            [
                bins_df,
                pd.DataFrame(
                    [{"path": "workflow/scripts/test/DUMMY_PATH/data/bin_2.fa"}]
                ),
            ],
            ignore_index=True,
        )

        bins_df = mbnu.get_bin_filename(bins_df)
        bins_df = mbnu.make_unduplicated_filenames(bins_df)
        bins_df = bins_df.reset_index(drop=True)

        # asserting that `make_unduplicated_filenames` correctly renamed the duplicated filenames
        pd.testing.assert_frame_equal(bins_df, expected_table)

    def test_copy_rename(self):

        # preparing input data
        bins_df = pd.read_csv(
            LIST_BINS_WITH_ONE_DUPLICATE_ONE_GZIPPED, sep="\t", names=["path"]
        )
        bins_df = mbnu.get_bin_filename(bins_df)
        bins_df = mbnu.make_unduplicated_filenames(bins_df)

        print(bins_df)

        # running the copy and rename operation
        destination_folder = "workflow/scripts/test/data/copied_and_renamed_bins"
        mbnu.copy_rename(bins_df, destination_folder)

        # asserting only expected files exist in destination_folder
        expected_files = ["a.fa", "bin_2_1.fa", "bin_2.fa", "bin_34.fa"]
        actual_files = os.listdir(destination_folder)

        self.assertEqual(
            sorted(actual_files),
            sorted(expected_files),
            f"Expected files {expected_files}, but found {actual_files} in {destination_folder}",
        )

        # cleaning up by removing all files from destination_folder after test
        for f in os.listdir(destination_folder):
            os.remove(os.path.join(destination_folder, f))

        # removing the destination folder itself
        os.rmdir(destination_folder)


if __name__ == "__main__":
    unittest.main()
