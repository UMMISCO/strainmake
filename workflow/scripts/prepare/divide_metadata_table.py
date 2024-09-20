"""
A script to split a metadata TSV table into a given number of subtables, keeping same samples in the
same table
"""

import argparse
import pandas as pd
import os
import math

def split_tsv(input_file, n, output_dir):
    """
    Splits the input TSV file into multiple smaller TSV files, ensuring that all rows
    with the same `sample_id` remain in the same output file
    """
    df = pd.read_csv(input_file, sep='\t')

    sample_ids = df['sample_id'].unique()

    num_samples = len(sample_ids)
    if n > num_samples:
        n = num_samples

    # splitting sample_ids into n groups, as evenly as possible
    group_size = math.ceil(num_samples / n)
    sample_id_groups = [sample_ids[i:i + group_size] for i in range(0, num_samples, group_size)]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, sample_group in enumerate(sample_id_groups):
        subset_df = df[df['sample_id'].isin(sample_group)]
        output_file = os.path.join(output_dir, f'output_part_{i+1}.tsv')
        subset_df.to_csv(output_file, sep='\t', index=False)
        print(f"File {output_file} generated with {len(subset_df)} rows")

def main():
    parser = argparse.ArgumentParser(description="Splits a TSV file into multiple smaller files while keeping the same sample IDs in the same file")
    parser.add_argument('input_file', help="Path to the input TSV file")
    parser.add_argument('n', type=int, help="Number of output files to generate")
    parser.add_argument('output_dir', help="Directory to save the output TSV files")

    args = parser.parse_args()

    split_tsv(args.input_file, args.n, args.output_dir)

if __name__ == "__main__":
    main()