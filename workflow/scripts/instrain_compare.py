""" 
This script will check which inStrain results it can use (those which inStrain did not fail) and will run `inStrain compare`
on them
"""

import os
import subprocess
from snakemake.script import snakemake

def check_results_folder(results_folder: list) -> dict:
    """
    Iterates over directories listed in `results_folder` and checks the content of their `results.txt` file.
    If it contains "SUCCESS" on the first line, adds this folder path to the dictionary under the key "run".
    If it contains "FAILED" on the first line, adds this folder path to the dictionary under the key "dont_run"
    """
    
    results = {"run": [], "dont_run": []}
    
    for folder in results_folder:
        results_file = os.path.join(folder, "results.txt")
        
        # checking if results.txt exists
        if os.path.isfile(results_file):
            # reading the first line of results.txt
            with open(results_file, "r") as file:
                first_line = file.readline().strip()

            # classifying the folder based on the first line of the file
            if first_line == "SUCCESS":
                results["run"].append(folder)
            elif first_line == "FAILED":
                results["dont_run"].append(folder)
    
    return results


def run_instrain_compare(run_or_not: dict, stb_file: str, output: str,
                         stdout: str, stderr: str):
    """ 
    Run `inStrain compare` using folders listed at key "run" of `run_or_not` dictionary
    """
    
    instrain_results_to_process = " ".join(run_or_not["run"])
    cmd = f"inStrain compare -i {instrain_results_to_process} --output {output} --database_mode -s {stb_file} > {stdout} 2> {stderr}"

    subprocess.run(cmd, shell=True)

run_or_not = check_results_folder(snakemake.input[0])
run_instrain_compare(run_or_not, snakemake.input[1], snakemake.output, snakemake.log[0], snakemake.log[1])