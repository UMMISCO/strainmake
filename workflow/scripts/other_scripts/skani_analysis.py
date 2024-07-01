#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform Skani analysis on bins.')
    parser.add_argument('--bins', required=True, choices=['refined', 'dereplicated'], 
                        help='Type of bins to analyze')
    parser.add_argument('--tmp', required=True, help='Temporary directory for intermediate files')
    parser.add_argument('--output_file', required=True, help='File to save the output results')
    parser.add_argument('--cpu', type=int, required=True, help='Number of CPU cores to use')
    parser.add_argument('--tsv_output', required=True, help='File to save the Skani matrix in TSV format')
    return parser.parse_args()

def copy_and_rename_bins_refined(src_dir, tmp_dir):
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    
    list_bins = []

    print("Copying and renaming bins files")
    for assembly in os.listdir(src_dir):
        assembly_dir = os.path.join(src_dir, assembly)
        print(f"Copying and renaming bins from {assembly} assembly (progress bar displays samples)")
        if os.path.isdir(assembly_dir):
            for sample in tqdm(os.listdir(assembly_dir)):
                sample_dir = os.path.join(assembly_dir, sample, 'final_bins')
                if os.path.exists(sample_dir):
                    for bin_file in os.listdir(sample_dir):
                        if bin_file.endswith('.fa'):
                            src_bin = os.path.join(sample_dir, bin_file)
                            new_bin_name = f"{assembly}.{sample}.{bin_file}"
                            dst_bin = os.path.join(tmp_dir, new_bin_name)
                            shutil.copy2(src_bin, dst_bin)
                            list_bins.append(dst_bin)
    
    print("Creating the list_bins.txt file")
    list_bins_path = os.path.join(tmp_dir, 'list_bins.txt')
    with open(list_bins_path, 'w') as f:
        for bin_path in list_bins:
            f.write(bin_path + '\n')
    
    return list_bins_path

def copy_and_rename_bins_dereplicated(src_dir, tmp_dir):
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    
    list_bins = []

    print("Copying and renaming dereplicated bins files")
    for assembly in os.listdir(src_dir):
        assembly_dir = os.path.join(src_dir, assembly, 'bins')
        if os.path.exists(assembly_dir):
            print(f"Copying and renaming bins from {assembly} assembly")
            for bin_file in tqdm(os.listdir(assembly_dir)):
                if bin_file.endswith('.fa'):
                    src_bin = os.path.join(assembly_dir, bin_file)
                    new_bin_name = f"{assembly}.{bin_file}"
                    dst_bin = os.path.join(tmp_dir, new_bin_name)
                    shutil.copy2(src_bin, dst_bin)
                    list_bins.append(dst_bin)
    
    print("Creating the list_bins.txt file")
    list_bins_path = os.path.join(tmp_dir, 'list_bins.txt')
    with open(list_bins_path, 'w') as f:
        for bin_path in list_bins:
            f.write(bin_path + '\n')

    return list_bins_path

def run_skani(list_bins_path, output_file, cpu):
    print("Running Skani analysis")
    cmd = ['skani', 'triangle', '--medium', '-t', str(cpu), '-l', list_bins_path, '-o', output_file, '--full-matrix']
    subprocess.run(cmd, check=True)

def read_phylip_lower_triangular(filepath):
    """
    A function to read the Skani results (Phylip lower triangular matrix) into a 
    numpy matrix and returns it as a pandas dataframe for having the dimensions
    names
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # get the number of elements
    num_elements = int(lines[0].strip())

    # initialize an empty numpy array
    matrix = np.zeros((num_elements, num_elements))

    bin_names = []

    # fill the numpy array with values
    for i in range(1, len(lines)):
        elements = lines[i].strip().split()
        values = list(map(float, elements[1:]))

        bin_name = os.path.basename(elements[0])
        bin_names.append(bin_name)
        
        matrix[i-1, :i] = values
        matrix[:i, i-1] = values

    skani_results = pd.DataFrame(matrix, index=bin_names, columns=bin_names)

    return skani_results

def main():
    args = parse_arguments()

    if args.bins == 'refined':
        src_dir = 'results/07_bins_refinement/binette'
        list_bins_path = copy_and_rename_bins_refined(src_dir, args.tmp)
    elif args.bins == 'dereplicated':
        src_dir = 'results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality'
        list_bins_path = copy_and_rename_bins_dereplicated(src_dir, args.tmp)

    run_skani(list_bins_path, args.output_file, args.cpu)

    # read the skani result and save as tsv
    skani_matrix = read_phylip_lower_triangular(args.output_file)
    skani_matrix.to_csv(args.tsv_output, sep='\t', index=True, header=True)

    print(f"Skani matrix saved to {args.tsv_output}")

if __name__ == "__main__":
    main()