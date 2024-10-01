#!/usr/bin/env python3
"""
Script for gathering binning results from multiple runs into a single target directory,
after processing with tools that split data across multiple directories. For example, after
having used `workflow/scripts/prepare/divide_metadata_table.py`.

The script checks that the target directory is not among the result directories.
It collects results from the "05_binning" and "07_bins_refinement" subdirectories of each result directory
and moves them to the corresponding folders in the target directory, maintaining the structure.
"""

import argparse
import os
import shutil
import logging

# configuring logging for the script
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def move_files(src_dir, dest_dir):
    """
    Moves all contents from src_dir to dest_dir. dest_dir should exist.
    """
    logging.info(f"Moving files from {src_dir} to {dest_dir}")

    for item in os.listdir(src_dir):
        src_path = os.path.join(src_dir, item)
        dest_path = os.path.join(dest_dir, item)
        if os.path.isdir(src_path):
            shutil.move(src_path, dest_path)
        else:
            shutil.move(src_path, dest_dir)

def gather_binning_results(target_dir, res_dirs, folder_name):
    """
    Gathers results from the specified subdirectory of each result directory into the corresponding target directory.
    
    Parameters:
    - `target_dir`: The root directory where the results will be gathered
    - `res_dirs`: A list of directories with individual binning results
    - `folder_name`: The subdirectory name to process (e.g., '05_binning' or '07_bins_refinement')
    """
    for res_dir in res_dirs:
        binning_dir = os.path.join(res_dir, folder_name)
        if not os.path.exists(binning_dir):
            logging.warning(f"No '{folder_name}' directory found in {res_dir}. Skipping...")
            continue
        
        logging.info(f"Processing directory {binning_dir}")
        
        dest_binning_dir = os.path.join(target_dir, folder_name)
        if not os.path.exists(dest_binning_dir):
            logging.info(f"Creating directory {dest_binning_dir}")
            os.makedirs(dest_binning_dir)
        
        move_files(binning_dir, dest_binning_dir)

def main():
    parser = argparse.ArgumentParser(description="Gather binning results from multiple runs into a single target directory")
    parser.add_argument('--target_dir', required=True, help="Directory where the results will be gathered")
    parser.add_argument('res_dirs', nargs='+', help="List of directories with individual binning results")

    args = parser.parse_args()

    target_dir = os.path.abspath(args.target_dir)
    res_dirs = [os.path.abspath(res_dir) for res_dir in args.res_dirs]

    # checking if target_dir is in res_dirs
    if target_dir in res_dirs:
        logging.error("Target directory cannot be one of the result directories.")
        exit(1)

    logging.info(f"Gathering results into {target_dir} from {len(res_dirs)} result directories.")

    # check for target binning directories (05_binning and 07_bins_refinement)
    for folder_name in ['05_binning', '07_bins_refinement']:
        full_target_dir = os.path.join(target_dir, folder_name)
        if not os.path.exists(full_target_dir):
            logging.error(f"{full_target_dir} does not exist.")
            exit(1)

        logging.info(f"Binning directory for target is considered to be {full_target_dir}.")

    # gathering the binning results from both 05_binning and 07_bins_refinement
    gather_binning_results(target_dir, res_dirs, '05_binning')
    gather_binning_results(target_dir, res_dirs, '07_bins_refinement')

if __name__ == "__main__":
    main()
