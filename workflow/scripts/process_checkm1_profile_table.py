""" 
Renaming columns of table produced by `checkm profile` command
"""

import pandas as pd
import argparse

def rename_columns(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')

    # detecting column names
    column_mapping = {}
    for col in df.columns:
        if ": mapped reads" in col:
            column_mapping[col] = "mapped reads"
        elif ": % mapped reads" in col:
            column_mapping[col] = "% mapped reads"
        elif ": % binned populations" in col:
            column_mapping[col] = "% binned populations"
        elif ": % community" in col:
            column_mapping[col] = "% community"

    # renaming columns
    df.rename(columns=column_mapping, inplace=True)

    # exporting DataFrame into a TSV
    df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="""Rename "checkm profile"'s table columns.""")
    parser.add_argument('--input_table', required=True, help='Path to the input TSV file.')
    parser.add_argument('output', help='Path to save the output TSV file.')

    args = parser.parse_args()

    rename_columns(args.input_table, args.output)

if __name__ == "__main__":
    main()
