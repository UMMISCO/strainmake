import os
import numpy as np
import pandas as pd

np.random.seed(999)

genome_dir = "genomes"
genome_files = os.listdir(genome_dir)
genome_files = [f for f in genome_files if f.endswith('.fna')]

assert len(genome_files) == 5, f"5 genomes are needed, what we have {genome_files}"

mean = 2
std_dev = 1
abundances = np.random.normal(loc=mean, scale=std_dev, size=5)

abundances = np.abs(abundances)
abundances = abundances / abundances.sum()
abundances_percent = abundances * 100

abundance_table_for_nanosim = pd.DataFrame({
    "Size": genome_files,
    "150000": abundances_percent # will generate 150 000 reads using the "abundances_percent" distribution of abundance
})

abundance_table_for_nanosim.to_csv("abundance_table_for_nanosim.tsv", sep='\t', index=False, header=True)