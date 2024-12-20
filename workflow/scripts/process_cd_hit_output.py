"""
This script processes the output from CD-HIT, identifies the best gene in each cluster based on length, and filters the input FASTA file to 
include only these best genes.
The script performs the following steps:
1. Parses command-line arguments to get the paths for the input CD-HIT output file, input FASTA file, output filtered table file, and output FASTA file.
2. Processes the CD-HIT output to create a DataFrame.
3. Identifies the best gene in each cluster based on the specified minimal length.
4. Saves the filtered table of best genes to the specified output file.
5. Opens the input FASTA file and filters it to include only the best genes, saving the result to the specified output FASTA file.
"""

import argparse
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import gzip

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

def open_genes_list(fasta: str) -> list:
    with gzip.open(fasta, "rt") as handle:
        return list(SeqIO.parse(handle, "fasta"))

def filter_genes_list(genes_list: list, best_genes: pd.DataFrame, output_fasta: str) -> None:
    """
    Filters a list of genes based on a DataFrame of best genes and writes the filtered genes to a FASTA file.

    Args:
        genes_list (list): A list of gene objects to be filtered.
        best_genes (pd.DataFrame): A DataFrame containing the best genes with a column 'sequence_name' to match against gene IDs.
        output_fasta (str): The file path where the filtered genes will be written in FASTA format.
    """

    filtered_genes = [gene for gene in genes_list if gene.id in best_genes['sequence_name'].values]
    print(len(filtered_genes))
    # Write the filtered genes to the specified output FASTA file
    SeqIO.write(filtered_genes, output_fasta, "fasta")

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
    
    # Open the genes list from the input FASTA file
    genes_list = open_genes_list(args.input_fasta)
    
    # Filter the genes list to get only the best genes
    filter_genes_list(genes_list, best_genes, args.output_fasta)

if __name__ == "__main__":
    main()