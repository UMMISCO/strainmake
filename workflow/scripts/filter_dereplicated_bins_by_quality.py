"""
This command line tool will read a CheckM2 quality report,
keep bins that pass users given parameters and copy the selected bins
in the output folder
"""

import argparse
import pandas as pd
import os
import shutil

def filter_bins(checkm_report, bins_directory, min_comp, max_cont, outdir):
    # load the TSV file
    df = pd.read_csv(checkm_report, sep='\t')

    print("Opened CheckM2 quality report")
    
    # filter the dataframe based on completeness and contamination criteria
    filtered_df = df[(df['Completeness'] >= min_comp) & (df['Contamination'] <= max_cont)]
    print("After filtration:")
    print(filtered_df)

    # if no bins pass the filtration, print a message and return
    if filtered_df.empty:
        print("No bins passed the filtration criteria.")
        return
    
    # make sure the output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # copy the selected bins to the output directory
    for bin_name in filtered_df['Name']:
        src_path = os.path.join(bins_directory, bin_name)
        dest_path = os.path.join(outdir, bin_name)
        if os.path.exists(src_path):
            shutil.copy(src_path, dest_path)
        else:
            print(f"Warning: {src_path} does not exist and will be skipped.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and copy dereplicated bins by quality criteria.")
    parser.add_argument('--checkm-report', required=True, help='Path to the quality report TSV file.')
    parser.add_argument('--bins-directory', required=True, help='Directory containing the bin files.')
    parser.add_argument('--min-comp', type=float, required=True, help='Minimum completeness threshold.')
    parser.add_argument('--max-cont', type=float, required=True, help='Maximum contamination threshold.')
    parser.add_argument('--outdir', required=True, help='Output directory for filtered bins.')
    
    args = parser.parse_args()
    
    filter_bins(args.checkm_report, args.bins_directory, args.min_comp, args.max_cont, args.outdir)
