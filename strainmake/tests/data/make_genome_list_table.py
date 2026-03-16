"""Build genome list and DNA type tables for NanoSim test data generation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


def generate_genome_tables(genome_dir: str = "genomes", seed: int = 999) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate genome list and DNA type tables for NanoSim based on the genome files in `genome_dir`. 
    The genome list table contains columns "Genome" and "Path", while the DNA type table contains columns "Genome", "Chromosome", and "DNA_type". 
    Both tables are returned as DataFrames.
    """
    np.random.seed(seed)
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith(".fna")]
    genomes_path = [os.path.abspath(os.path.join(genome_dir, genome)) for genome in genome_files]

    genomes_list_for_nanosim = pd.DataFrame({
        "Genome": genome_files,
        "Path": genomes_path,
    })

    dna_type_list_for_nanosim = pd.DataFrame({
        "Genome": genome_files,
        "Chromosome": [1] * len(genome_files),
        "DNA_type": ["circular"] * len(genome_files),
    })

    return genomes_list_for_nanosim, dna_type_list_for_nanosim


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate genome list and DNA type tables for NanoSim.")
    parser.add_argument("--genome-dir", default="genomes")
    parser.add_argument("--genomes-output", default="genomes_list_for_nanosim")
    parser.add_argument("--dna-type-output", default="dna_type_list_for_nanosim.tsv")
    args = parser.parse_args(argv)

    genomes_df, dna_df = generate_genome_tables(genome_dir=args.genome_dir)

    genomes_output = Path(args.genomes_output)
    genomes_output.parent.mkdir(parents=True, exist_ok=True)
    dna_output = Path(args.dna_type_output)
    dna_output.parent.mkdir(parents=True, exist_ok=True)

    genomes_df.to_csv(genomes_output, sep="\t", index=False, header=False)
    dna_df.to_csv(dna_output, sep="\t", index=False, header=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
