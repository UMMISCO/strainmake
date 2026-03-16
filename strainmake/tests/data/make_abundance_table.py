"""Build abundance table for NanoSim test data generation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


def generate_abundance_table(genome_dir: str = "genomes", seed: int = 999) -> pd.DataFrame:
    """
    Generate an abundance table for NanoSim based on the genome files in `genome_dir`. The abundances are randomly generated from a 
    normal distribution, normalized to sum to 100%, and returned as a DataFrame with columns "Size" and "150000".
    """
    np.random.seed(seed)
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith(".fna")]
    assert len(genome_files) == 5, f"5 genomes are needed, got: {genome_files}"

    abundances = np.random.normal(loc=2, scale=1, size=5)
    abundances = np.abs(abundances)
    abundances = abundances / abundances.sum()
    abundances_percent = abundances * 100

    return pd.DataFrame({
        "Size": genome_files,
        "150000": abundances_percent,
    })


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate abundance table for NanoSim.")
    parser.add_argument("--genome-dir", default="genomes")
    parser.add_argument("--output", default="abundance_table_for_nanosim.tsv")
    args = parser.parse_args(argv)

    df = generate_abundance_table(genome_dir=args.genome_dir)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False, header=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
