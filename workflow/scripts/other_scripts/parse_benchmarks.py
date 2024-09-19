"""
Script to produce a unique TSV made of the results of the Snakemake's
benchmark directive
"""

import argparse
import os
import re
import pandas as pd

def get_sample_name_from_path(path: str):
    """
    Extracts the sample name from the benchmark file path.
    Assumes that the sample name is between a "/" and a "." in the filename.
    Example: benchmarks/03_assembly/megahit/ADDE_TTK-09.rename.benchmark.txt -> "ADDE_TTK-09"
    """
    global_rules = ['get_bowtie_index', 'bins_refinement', 'checkm2.db', 'list_unduplicated_filenames', 
                    'list.benchmark', 'copy_and_rename', 'assessment_and_filtration', '.concatenate',
                    'concatenate_ref_genes', 'stb.', 'merge_reports', 'prodigal', 'plot_reports',
                    'dereplication.']
    
    # checking if any of the global rules are in the path
    if any(rule in path for rule in global_rules):
        # this rule applies globally, not per sample
        return 'global'

    basename = os.path.basename(path)  # Get file name from path
    sample_name = basename.split('.')[0]  # Extract part before the first dot
    return sample_name

def list_tables(path: str):
    """
    Finds all benchmark tables in the given folder.
    It returns a list of paths with .txt extension (benchmark files).
    """
    benchmark_files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".txt"):  # Only consider .txt files (assuming benchmarks)
                benchmark_files.append(os.path.join(root, file))
    return benchmark_files

def extract_part_and_tool(path: str):
    """
    Extracts 'part' and 'tool' from the benchmark file path.
    'part' is the directory that starts with a number (e.g., 03_assembly).
    'tool' is the directory that follows 'part' (e.g., megahit).
    """
    # splitting the path into its components
    parts = path.split(os.sep)
    
    # finding the 'part' (the first directory that starts with a number)
    part = None
    tool = None
    
    for i, part_candidate in enumerate(parts):
        if re.match(r"^\d+_", part_candidate):
            part = part_candidate
            # the tool is the next directory after 'part'
            if i + 1 < len(parts):
                tool = parts[i + 1]
            break
    
    if not part or not tool:
        raise ValueError(f"Could not extract part and tool from path: {path}")
    
    return part, tool

def concatenate_benchmarks(list_tables: list):
    """
    Produce a single DataFrame by concatenating all benchmark tables.
    Adds columns for sample name, directory, part, and tool.
    """
    all_dataframes = []

    for table_path in list_tables:
        sample_name = get_sample_name_from_path(table_path)
        part, tool = extract_part_and_tool(table_path)
        
        # loading the benchmark file into a DataFrame
        df = pd.read_csv(table_path, sep='\t')  # Assuming the benchmark files are TSV

        # adding columns for the sample, path, part, and tool
        df['sample'] = sample_name
        df['path'] = table_path
        df['part'] = part
        df['tool'] = tool

        all_dataframes.append(df)

    # concatenating all the DataFrames into one
    final_df = pd.concat(all_dataframes, ignore_index=True)
    
    return final_df

def main():
    """ 
    CLI logic
    """
    parser = argparse.ArgumentParser(description="Concatenate Snakemake benchmark results into a single TSV")
    parser.add_argument("benchmark_dir", help="Directory containing benchmark .txt files")
    parser.add_argument("output_tsv", help="Output file to save the concatenated TSV")
    args = parser.parse_args()

    # getting list of all benchmark files
    tables = list_tables(args.benchmark_dir)

    # concatenating the benchmarks into a single DataFrame
    final_df = concatenate_benchmarks(tables)

    # saving the final DataFrame as a TSV file
    final_df.to_csv(args.output_tsv, sep='\t', index=False)

    print(f"Concatenated benchmark results saved to {args.output_tsv}")

if __name__ == "__main__":
    main()