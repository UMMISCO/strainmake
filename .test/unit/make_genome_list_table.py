import os
import numpy as np
import pandas as pd

np.random.seed(999)

genome_dir = "genomes"
genome_files = os.listdir(genome_dir)
genome_files = [f for f in genome_files if f.endswith('.fna')]
genomes_path = [os.path.abspath(os.path.join(genome_dir, genome)) for genome in genome_files]

genomes_list_for_nanosim = pd.DataFrame({
    "Genome": genome_files,
    "Path": genomes_path
})

# for the DNA_TYPE_LIST table
chromosome_name = [1] * 5
dna_type = ["circular"] * 5

dna_type_list_for_nanosim = pd.DataFrame({
    "Genome": genome_files,
    "Chromosome": chromosome_name,
    "DNA_type": dna_type
})

genomes_list_for_nanosim.to_csv("genomes_list_for_nanosim", sep='\t', index=False, header=False)
dna_type_list_for_nanosim.to_csv("dna_type_list_for_nanosim.tsv", sep='\t', index=False, header=False)
