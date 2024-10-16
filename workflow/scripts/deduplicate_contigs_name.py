#!/usr/bin/env python3
""" 
Ensure there is not duplicated headers in FASTA stored in a given folder, by adding in every header
the filename
Script to identify if two or more bins, stored in a given folder, share the same contig names. If so, 
it deduplicate them in order to have unique contigs name
"""

import os
import logging
import argparse

# configuring logging for the script
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def list_fasta(dir: str):
    """ 
    Returns a list of FASTA files stored in `dir`
    """
    fasta = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            # only selecting FASTA files
            if file.endswith(".fa"):
                fasta.append(os.path.join(root, file))
    return fasta

def rename_fasta_headers(fasta: list):
    """ 
    Renames the headers of all FASTA files in the list `fasta` by 
    adding the filename as a prefix to each header.
    """
    for file_path in fasta:
        # Extract the base filename without extension
        base_name = os.path.basename(file_path).split('.')[0]

        logging.info(f"Treating file: {base_name}")
        
        # Create a new file to write the updated fasta
        new_file_path = f"{file_path}_renamed.fa"
        
        with open(file_path, 'r') as infile, open(new_file_path, 'w') as outfile:
            for line in infile:
                # If the line is a header, modify it
                if line.startswith(">"):
                    # Strip the newline and prefix with the filename
                    new_header = f">{base_name}_{line[1:].strip()}\n"
                    outfile.write(new_header)
                else:
                    # Write the sequence line as-is
                    outfile.write(line)

        # we delete the original file
        os.remove(file_path)

        # we rename the FASTA with the new contigs name as the original FASTA
        os.rename(new_file_path, file_path)
        logging.info(f"New file written in {os.path.basename(file_path)}")

def main():
    """ 
    Main program logic
    """
    parser = argparse.ArgumentParser(description="Rename sequences of all FASTA stored in a given folder.")
    parser.add_argument('fasta_dir', help="Folder where the FASTA files are stored.")

    args = parser.parse_args()

    fastas = list_fasta(args.fasta_dir)
    rename_fasta_headers(fastas)

if __name__ == "__main__":
    main()