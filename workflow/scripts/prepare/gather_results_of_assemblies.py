#!/usr/bin/env python3
"""
Script for gathering results from multiple assembly runs into a single target directory, for 
example after having split your data using `workflow/scripts/prepare/divide_metadata_table.py`.
    
The script checks that the target directory is not among the result directories.
It collects results from the "03_assembly" subdirectory of each result directory and moves them to
corresponding folders in the target directory, maintaining structure like "hybridspades", "megahit", 
and "LR/metaflye".
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

def gather_results(target_dir, res_dirs):
    """
    Gathers results from multiple assembly directories into the target directory.
    Ensures that the content of `hybridspades`, `megahit`, and `metaflye` (inside `LR`) 
    are moved into their respective locations in the target directory.
    """
    for res_dir in res_dirs:
        assembly_dir = os.path.join(res_dir, '03_assembly')
        if not os.path.exists(assembly_dir):
            logging.warning(f"No '03_assembly' directory found in {res_dir}. Skipping...")
            continue
        
        logging.info(f"Processing directory {assembly_dir}")
        
        # handling hybridspades/metaspades and megahit folders
        for subfolder in ['hybridspades', 'megahit', 'metaspades']:
            src_subdir = os.path.join(assembly_dir, subfolder)
            if os.path.exists(src_subdir):
                dest_subdir = os.path.join(target_dir, subfolder)
                move_files(src_subdir, dest_subdir)
        
        # handling LR/metaflye
        lr_dir = os.path.join(assembly_dir, 'LR', 'metaflye')
        if os.path.exists(lr_dir):
            dest_lr_metaflye = os.path.join(target_dir, 'LR', 'metaflye')
            move_files(lr_dir, dest_lr_metaflye)

def main():
    parser = argparse.ArgumentParser(description="Gather results from multiple assembly runs into a single target directory")
    parser.add_argument('--target_dir', required=True, help="Directory where the results will be gathered")
    parser.add_argument('res_dirs', nargs='+', help="List of directories with individual assembly results")

    args = parser.parse_args()

    target_dir = os.path.abspath(args.target_dir)
    res_dirs = [os.path.abspath(res_dir) for res_dir in args.res_dirs]

    # checking if target_dir is in res_dirs
    if target_dir in res_dirs:
        logging.error("Target directory cannot be one of the result directories.")
        exit(1)

    logging.info(f"Gathering results into {target_dir} from {len(res_dirs)} result directories.")

    target_dir = os.path.join(target_dir, "03_assembly")
    if not os.path.exists(target_dir):
        logging.error(f"{target_dir} does not exist.")
        exit(1)

    logging.info(f"Assemblies directory of target is considered to be {target_dir}.")
    
    # gathering the results
    gather_results(target_dir, res_dirs)

if __name__ == "__main__":
    main()