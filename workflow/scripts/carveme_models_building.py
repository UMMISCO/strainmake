""" 
Using a script to call carveme to build metabolic model of several organisms, since the 
carveme's `-r` command is buggy (https://github.com/cdanielmachado/carveme/issues/99)
"""

#!/usr/bin/env python3

import argparse
import os
import logging
from concurrent.futures import ThreadPoolExecutor

# configuring logging for the script
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def list_mag(dir: str, ext: str, verbose: bool) -> list:
    """
    Function to list all the genomes (MAG) files in the directory
    with the extension `ext`
    """
    mag_list = []
    # walking through the directory tree
    for root, dirs, files in os.walk(dir):
        for file in files:
            # checking if the file ends with the specified extension
            if file.endswith(ext):
                mag_path = os.path.join(root, file)
                mag_list.append(mag_path)
                if verbose:
                    logging.info(f"Found MAG file: {mag_path}")
    return mag_list

def carveme_carve(mag_list: list, out_dir: str, cpu: int, verbose: bool, dryrun: bool = False) -> None:
    """ 
    Function to infer metabolic models using CarveMe on 
    the genomes (MAG) given in the `mag_list`
    Models are saved in the `out_dir` directory
    """
    def carve(mag):
        # getting the base name of the MAG file
        mag_basename = os.path.basename(mag)
        # degining the output model file path
        output_file = os.path.join(out_dir, f"{mag_basename}.xml")

        if verbose:
            logging.info(f"Processing MAG file: {mag}")
            logging.info(f"Output will be saved to: {output_file}")

        # running the CarveMe command to infer the metabolic model
        command = f"carve --solver gurobi --dna --output {output_file} {mag}"

        if verbose:
            logging.info(f"Running command: {command}")

        if dryrun:
            print(command)
        else:
            os.system(command)

    # using ThreadPoolExecutor to parallelize the carving process
    with ThreadPoolExecutor(max_workers=cpu) as executor:
        executor.map(carve, mag_list)

def carveme_merge_community(communities: str, out_dir: str, verbose: bool, dryrun: bool = False) -> None:
    """ 
    Function to merge the individual metabolic models of the organisms that were inferred using `carve`
    """
    output_dir = os.path.join(out_dir, "community.xml")
    command = f"merge_community {communities}/* -o {output_dir}"

    if verbose:
        logging.info(f"Running command: {command}")

    if dryrun:
        print(command)
    else:
        os.system(command)


def main():
    parser = argparse.ArgumentParser(description="Infer metabolic models using CarveMe")
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # Subcommand for carveme_carve
    parser_carve = subparsers.add_parser("carve", help="Infer metabolic models using CarveMe")
    parser_carve.add_argument("-i", "--input_dir", required=True, help="Input directory containing MAG files")
    parser_carve.add_argument("-o", "--output_dir", required=True, help="Output directory to save the models")
    parser_carve.add_argument("-e", "--extension", default=".fa", help="Extension of MAG files (default: .fa)")
    parser_carve.add_argument("-c", "--cpu", type=int, default=1, help="Number of CPU cores to use (default: 1)")
    parser_carve.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    # Subcommand for carveme_merge_community
    parser_merge = subparsers.add_parser("merge", help="Merge individual metabolic models into a community model")
    parser_merge.add_argument("-i", "--input_dir", required=True, help="Input directory containing individual models")
    parser_merge.add_argument("-o", "--output_dir", required=True, help="Output directory to save the community model")
    parser_merge.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.command == "carve":
        mag_list = list_mag(args.input_dir, args.extension, args.verbose)
        carveme_carve(mag_list, args.output_dir, args.cpu, args.verbose)
    elif args.command == "merge":
        carveme_merge_community(args.input_dir, args.output_dir, args.verbose)

if __name__ == "__main__":
    main()