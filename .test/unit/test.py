import unittest
import subprocess
import os
import re
import pandas as pd

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

if __name__ == "__main__":
    unittest.main()