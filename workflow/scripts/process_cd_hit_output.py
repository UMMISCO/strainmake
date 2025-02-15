"""
This script processes the output from CD-HIT, identifies the best gene in each cluster based on length, and filters the input FASTA file to 
include only these best genes.
"""

import argparse
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import gzip
import time
import os

def process_cd_hit_output(input_file: str) -> pd.DataFrame:
    """
    Process the output of CD-HIT and return a DataFrame with the cluster information.
    """
    
    # open the file and read the lines
    with open(input_file, 'r') as f:
        lines = f.readlines()

    data = []

    # iterate through each line in the file
    for line in tqdm(lines, desc="Processing CD-HIT output"):
        # check if the line starts with '>Cluster' to get the cluster number
        if line.startswith('>Cluster'):
            cluster_number = line.split()[1]
        else:
            parts = line.split()
            sequence_length = parts[1].replace('nt,', '')
            sequence_name = parts[2].replace('...', '').lstrip('>')
            data.append([cluster_number, sequence_name, sequence_length])

    # create a DataFrame with the collected data
    df = pd.DataFrame(data, columns=['cluster_number', 'sequence_name', 'sequence_length'])
    
    # print the number of clusters processed
    print(f"Number of clusters processed: {df['cluster_number'].nunique()}")
    # print the total number of sequences processed
    print(f"Total number of sequences processed: {len(df)}")
    
    return df

def identify_best_gene(df: pd.DataFrame, minimal_len: int = 0) -> pd.DataFrame:
    """
    Identify the best gene in each cluster based on the length of the sequence.
    If minimal_len > 0, discard clusters whose best gene has a length < minimal_len.
    """
    # convert sequence_length to numeric type
    df['sequence_length'] = pd.to_numeric(df['sequence_length'])
    
    # group by cluster_number and get the sequence with the maximum length in each group
    best_genes = df.loc[df.groupby('cluster_number')['sequence_length'].idxmax()]
    
    # if minimal_len > 0, filter out clusters with best gene length < minimal_len
    if minimal_len > 0:
        best_genes = best_genes[best_genes['sequence_length'] >= minimal_len]
    
    return best_genes

def filter_genes_list(best_genes: pd.DataFrame, input_fasta: str, output_fasta: str) -> None:
    """
    Filters a list of genes based on a DataFrame of best genes and writes the filtered genes to a FASTA file.

    Args:
        best_genes (pd.DataFrame): A DataFrame containing a column 'sequence_name' with the names of the sequences to be filtered.
        input_fasta (str): Path to the input FASTA file containing all sequences.
        output_fasta (str): Path to the output FASTA file where the filtered sequences will be written.

    Returns:
        None
    """

    # Export the sequence_name column to a temporary TXT file without header
    temp_txt_file = "temp_sequence_names.txt"
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

    # Export the sequence_name column to a temporary TXT file without header
    temp_txt_file = os.path.join(os.path.dirname(output_fasta), "temp_sequence_names.txt")
    best_genes['sequence_name'].to_csv(temp_txt_file, index=False, header=False)

    # seqkit command to extract these sequences from the input FASTA file
    seqkit_cmd = f"seqkit grep -f {temp_txt_file} {input_fasta} -o {output_fasta}"
    print(f"Running seqkit command: {seqkit_cmd}")
    os.system(seqkit_cmd)

    # Remove the temporary TXT file
    os.remove(temp_txt_file)

def main():

    parser = argparse.ArgumentParser(description="Process CD-HIT output and identify the best gene in each cluster.")
    parser.add_argument("input_file", type=str, help="Path to the CD-HIT output file")
    parser.add_argument("input_fasta", type=str, help="Path to the input fasta.gz file")
    parser.add_argument("output_file", type=str, help="Path to the output file where the filtered table will be saved")
    parser.add_argument("output_fasta", type=str, help="Path to the output FASTA file containing only the best genes")
    parser.add_argument("--minimal_len", type=int, default=0, help="Minimal length of the best gene to consider (default: 0, no filtering if set to 0)")
    
    args = parser.parse_args()
    
    df = process_cd_hit_output(args.input_file)
    print(df)
    
    best_genes = identify_best_gene(df, args.minimal_len)
    print("Best genes in each cluster:")
    print(best_genes)
    
    # Save the filtered table to the specified output file
    best_genes.to_csv(args.output_file, index=False)
        
    # Filter the genes to get only the best genes
    filter_genes_list(best_genes, args.input_fasta, args.output_fasta)

if __name__ == "__main__":
    main()