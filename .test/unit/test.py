import unittest
import subprocess
import os
import re
import shutil
import pyfastx
import pandas as pd
from pathlib import Path

class TestPipeline(unittest.TestCase):
    """
    A set of unit tests for the Snakemake pipeline
    """

    @classmethod
    def setUpClass(cls):
        """
        Set up the environment before any test cases are run
        """
        # one folder for pipeline configs for each testing case
        cls.case1 = "1._SR_only"
        os.makedirs(cls.case1, exist_ok=True)

        cls.case2 = "2._LR_only"
        os.makedirs(cls.case2, exist_ok=True)

        cls.case3 = "3._ALL_seq"
        os.makedirs(cls.case3, exist_ok=True)

        cls.case4 = "4._SR_fail_if_LR_only"
        os.makedirs(cls.case4, exist_ok=True)

        cls.case5 = "5._LR_fail_if_SR_only"
        os.makedirs(cls.case5, exist_ok=True)

        cls.case6 = "6._SR_not_paired"
        os.makedirs(cls.case6, exist_ok=True)

        cls.case7 = "7._ALL_already_preprocessed"
        os.makedirs(cls.case6, exist_ok=True)

        cls.case8 = "8._ALL_LR_FASTA"
        os.makedirs(cls.case8, exist_ok=True)
        
    @classmethod
    def tearDownClass(cls):
        """
        """
        pass

    def parse_snakemake_dryrun_output(self, output: str) -> pd.DataFrame:
        """
        Parse the job table from Snakemake dry-run output and return it as a DataFrame
        """
        # fing the lines corresponding to the job table
        job_table_lines = []
        in_table = False

        for line in output.splitlines():
            # identify where is the start of the table
            if line.startswith("job") and "count" in line:
                in_table = True
                continue
            # detect the end of the table
            if in_table and line.strip().startswith("total "):
                job_table_lines.append(line.strip())
                break
            # collect lines within the table
            if in_table:
                job_table_lines.append(line.strip())

        # then, create a DataFrame from the job table lines
        job_data = []
        # to exclude header and total row
        for job_line in job_table_lines[1:-1]: 
            # to split on multiple spaces
            parts = re.split(r'\s{2,}', job_line)  
            if len(parts) == 2:
                job_data.append(parts)

        df = pd.DataFrame(job_data, columns=["job", "count"])
        df["count"] = pd.to_numeric(df["count"])

        return df

    def test_only_short_reads_dryrun(self):
        """
        Test everything would run well based on dry-run
        """
        # run Snakemake with the specified configuration file for short reads
        config_path = os.path.join(self.case1, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was successful
        self.assertEqual(result.returncode, 0, f"Snakemake failed: {result.stderr} {result.stdout}")

        # checking number of tasks for a rule is consistent
        snakemake_jobs = self.parse_snakemake_dryrun_output(result.stdout)
        reads_mapping_tasks = snakemake_jobs[snakemake_jobs['job'] == 'reads_mapping']
        self.assertEqual(int(reads_mapping_tasks['count'].values[0]), 2, "Job 'reads_mapping' count is not 2")

    def test_only_long_reads_dryrun(self):
        """
        Test everything would run well based on dry-run
        """
        # run Snakemake with the specified configuration file for long reads
        config_path = os.path.join(self.case2, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was successful
        self.assertEqual(result.returncode, 0, f"Snakemake failed: {result.stderr} {result.stdout}")

    def test_all_reads_dryrun(self):
        """
        Test everything would run well based on dry-run
        """
        # run Snakemake with the specified configuration file for short and long reads
        config_path = os.path.join(self.case3, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was successful
        self.assertEqual(result.returncode, 0, f"Snakemake failed: {result.stderr} {result.stdout}")

        # checking number of tasks for rules is consistent
        snakemake_jobs = self.parse_snakemake_dryrun_output(result.stdout)
        reads_mapping_tasks = snakemake_jobs[snakemake_jobs['job'] == 'reads_mapping']
        self.assertEqual(int(reads_mapping_tasks['count'].values[0]), 3, "Job 'reads_mapping' count is not 3")

        reads_mapping_LR_tasks = snakemake_jobs[snakemake_jobs['job'] == 'reads_mapping_LR']
        self.assertEqual(int(reads_mapping_LR_tasks['count'].values[0]), 1, "Job 'reads_mapping_LR' count is not 1")

    def test_given_config_dryrun(self):
        """
        Ensure the default config would run
        """
        # run Snakemake with the default config
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was successful
        self.assertEqual(result.returncode, 0, f"Snakemake failed: {result.stderr} {result.stdout}")

        # checking number of tasks for a rule is consistent
        snakemake_jobs = self.parse_snakemake_dryrun_output(result.stdout)
        reads_mapping_tasks = snakemake_jobs[snakemake_jobs['job'] == 'reads_mapping']
        self.assertEqual(int(reads_mapping_tasks['count'].values[0]), 3, "Job 'reads_mapping' count is not 3")

    def test_short_reads_fail_if_long_reads_only_dryrun(self): 
        """
        Pipeline should fail if the user wants to use short read assemblers but only have
        long reads sequences
        """
        # run Snakemake with the specified configuration file
        config_path = os.path.join(self.case4, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was not successful
        self.assertEqual(result.returncode, 1, f"Test failed: {result.stderr} {result.stdout}")

        # check if a ValueError was raised for this problem
        self.assertIn("ValueError", result.stdout)
        
    def test_long_reads_fail_if_short_reads_only_dryrun(self): 
        """
        Pipeline should fail if the user wants to use long read assemblers but only have
        short reads sequences
        """
        # run Snakemake with the specified configuration file
        config_path = os.path.join(self.case5, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was not successful
        self.assertEqual(result.returncode, 1, f"Test failed: {result.stderr} {result.stdout}")

        # check if a ValueError was raised for this problem
        self.assertIn("ValueError", result.stdout)

    def test_if_fail_if_no_pair_dryrun(self): 
        """
        Pipeline should fail if only one (eg. R1) of SR FASTQ is given since we should have
        R1 and R2 FASTQ in this case
        """
        # run Snakemake with the specified configuration file
        config_path = os.path.join(self.case6, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )
        
        # check if Snakemake command was not successful
        self.assertEqual(result.returncode, 1, f"Test failed: {result.stderr} {result.stdout}")

        # check if a ValueError was raised for this problem
        self.assertIn("ValueError", result.stdout)

    def test_reads_already_preprocessed_implementation(self):
        """
        Method for testing if the method implemented for taking into account already preprocessed 
        samples works well
        """

        prepare_results_folder_script = "../../workflow/scripts/prepare/already_preprocessed_seq.py"
        metadata_table = os.path.join(self.case7, "metadata.tsv")

        # testing if if the results directory folder is well created, with good symlinks
        result = subprocess.run(
            [
                "python3",
                prepare_results_folder_script,
                "--results_dir",
                "../../results",
                metadata_table
            ],
            capture_output=True, text=True
        )
        
        self.assertEqual(result.returncode, 0, f"The script did not run successfully: {result.stderr} {result.stdout}")

        # we should have a `results` folder and  three symlinks:
        # results/02_preprocess/bowtie2/SAMPLE1_1.clean.fastq.gz -> /path/to/pipeline/.test/unit/data/fake_illumina_R1.SAMPLE1.fastq.gz
        # results/02_preprocess/bowtie2/SAMPLE1_2.clean.fastq.gz -> /path/to/pipeline/.test/unit/data/fake_illumina_R2.SAMPLE1.fastq.gz
        # results/02_preprocess/fastp_long_read/SAMPLE1_1.fastq.gz -> /path/to/pipeline/.test/unit/data/fake_nanopore_sample0_aligned_reads.SAMPLE1.fastq.gz

        # test if "results" exist
        self.assertTrue(os.path.exists("../../results"), "The 'results' folder has not been created")

        # checking if results/02_preprocess/bowtie2/SAMPLE1_1.clean.fastq.gz, results/02_preprocess/bowtie2/SAMPLE1_2.clean.fastq.gz
        # and results/02_preprocess/fastp_long_read/SAMPLE1_1.fastq.gz are symlinks
        expected_symlinks = [
            "../../results/02_preprocess/bowtie2/SAMPLE1_1.clean.fastq.gz",
            "../../results/02_preprocess/bowtie2/SAMPLE1_2.clean.fastq.gz",
            "../../results/02_preprocess/fastp_long_read/SAMPLE1_1.fastq.gz",
        ]

        for symlink in expected_symlinks:
            symlink_path = Path(symlink)

            # checking if the symlink exists
            self.assertTrue(symlink_path.exists(), f"Symlink {symlink} does not exist")

            # checking if it's actually a symlink
            self.assertTrue(symlink_path.is_symlink(), f"{symlink} is not a symlink")

        # we then run the pipeline in dry mode and test that there is no rule of preprocessing
        # (since here we take these samples as already preprocessed)
        config_path = os.path.join(self.case7, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p",
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile.already_preprocessed_seq.smk",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )   

        # check if Snakemake command was successful
        self.assertEqual(result.returncode, 0, f"Snakemake failed: {result.stderr} {result.stdout}")
        
        snakemake_jobs = self.parse_snakemake_dryrun_output(result.stdout)
        # list of preprocessing jobs that shouldn't appear in the dry-run
        preprocessing_jobs = [
            "fastqc_before_preprocessing", 
            "fastp", 
            "host_decontamination", 
            "fastqc_after_preprocessing", 
            "fastp_long_read"
        ]

        # asserting that none of the preprocessing jobs are scheduled
        for index, row in snakemake_jobs.iterrows():
            self.assertNotIn(row["job"], preprocessing_jobs, 
                             f"Preprocessing job {row['job']} should not be scheduled for already preprocessed samples\n{snakemake_jobs}")

        # cleaning by removing symlinks and the "results" directory
        result = subprocess.run(
            [
                "python3",
                prepare_results_folder_script,
                "--clean",
                "--results_dir",
                "../../results",
                metadata_table
            ])  

        shutil.rmtree("../../results")

    def test_lr_reads_fasta(self):
        """
        Method for testing if the implementation of FASTA managing for LR samples works well
        """

        # preparing the results directory such as reads were already preprocessed
        prepare_results_folder_script = "../../workflow/scripts/prepare/already_preprocessed_seq.py"
        metadata_table = os.path.join(self.case8, "metadata.tsv")

        # testing if if the results directory folder is well created, with good symlinks
        result = subprocess.run(
            [
                "python3",
                prepare_results_folder_script,
                "--results_dir",
                "../../results",
                metadata_table
            ],
            capture_output=True, text=True
        )

        subprocess.run(
                [
                    "python3",
                    prepare_results_folder_script,
                    "--results_dir",
                    "../../results",
                    "--long-reads-seq-format",
                    "fasta",
                    metadata_table
                ],
                capture_output=True, text=True
            )
        
        # run Snakemake with the specified configuration file for short and long reads
        config_path = os.path.join(self.case8, "config.yaml")
        result = subprocess.run(
            [
                "snakemake", 
                "-c", "4", 
                "-p", 
                "--conda-frontend", "conda", 
                "--use-conda", 
                "--configfile", config_path,
                "-s", "../../workflow/Snakefile.already_preprocessed_seq.smk",
                "--directory", "../../",
                "-n"
            ],
            capture_output=True, text=True
        )

        # cleaning by removing symlinks and the "results" directory
        subprocess.run(
            [
                "python3",
                prepare_results_folder_script,
                "--clean",
                "--results_dir",
                "../../results",
                metadata_table
            ])  

        shutil.rmtree("../../results")
        
        # check if Snakemake command was successful
        self.assertEqual(result.returncode, 0, f"Snakemake failed: {result.stderr} {result.stdout}")

        # # checking number of tasks for rules is consistent
        # snakemake_jobs = self.parse_snakemake_dryrun_output(result.stdout)
        # reads_mapping_tasks = snakemake_jobs[snakemake_jobs['job'] == 'reads_mapping']
        # self.assertEqual(int(reads_mapping_tasks['count'].values[0]), 3, "Job 'reads_mapping' count is not 3")

        # reads_mapping_LR_tasks = snakemake_jobs[snakemake_jobs['job'] == 'reads_mapping_LR']
        # self.assertEqual(int(reads_mapping_LR_tasks['count'].values[0]), 1, "Job 'reads_mapping_LR' count is not 1")
    
    ################ Utils

    def test_sequences_renaming_in_fasta(self):
        """ 
        Method to ensure the method to rename headers in FASTA works well
        """
        # test FASTA we will use
        os.makedirs("fasta/tmp")
        fasta_test = "fasta/tmp/test_tmp_for_test.fa"
        shutil.copy2("fasta/test.fa", fasta_test)

        fa = pyfastx.Fasta(fasta_test)

        headers_before_renaming = ['contig_1', 'contig_2', 'contig_3']

        self.assertEqual(headers_before_renaming, list(fa.keys()))

        # cleaning index
        os.remove(f"{fasta_test}.fxi")

        # renaming step using the script
        script = "../../workflow/scripts/deduplicate_contigs_name.py"
        result = subprocess.run(
            [
                "python3",
                script,
                "fasta/tmp"
            ],
            capture_output=True, text=True
        )

        # check if the command was successful
        self.assertEqual(result.returncode, 0, f"FASTA headers renaming failed")

        # checking how the FASTA headers were renamed
        fa_new = pyfastx.Fasta(fasta_test)
        headers_after_renaming = ['test_tmp_for_test_contig_1', 'test_tmp_for_test_contig_2', 'test_tmp_for_test_contig_3']
        headers_after_renaming_test = list(fa_new.keys())

        self.assertEqual(headers_after_renaming, headers_after_renaming_test, f"Renamed sequences do not match what is expected (what we have: {headers_after_renaming_test})")

        # cleaning
        shutil.rmtree("fasta/tmp")
        

if __name__ == "__main__":
    unittest.main()