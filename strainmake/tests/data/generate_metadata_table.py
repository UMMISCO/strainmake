"""
Automatically generate TSV table with paths to FASTQ/FASTA files for pipeline tests.

Migrated and slightly refactored from `.test/unit/generate_metadata_table.py`.
"""

from __future__ import annotations

import argparse
import os
import re
from enum import Enum
from pathlib import Path

import pandas as pd


class MetadataType(str, Enum):
    SR = "SR"
    LR = "LR"
    ALL = "all"


def extract_sample_id(filename: str, pattern: str) -> str | None:
    """Extract sample ID from filename using the provided regex pattern."""
    match = re.search(pattern, filename)
    if match:
        return match.group(1)
    return None


def generate_metadata(directory: str, metadata_type: MetadataType | str, lr_format: str = "fastq") -> pd.DataFrame:
    """Return metadata DataFrame for test sequencing files found in `directory`."""
    metadata_type_value = metadata_type.value if isinstance(metadata_type, MetadataType) else str(metadata_type)

    sr_pattern = r"fake_illumina_R[12]\.(.*)\.fastq\.gz"
    lr_pattern = fr"fake_nanopore_.*\.(.*)\.{lr_format}\.gz"

    files = [f for f in os.listdir(directory) if f.endswith((".fastq.gz", ".fasta.gz"))]

    metadata: list[list[str]] = []
    for seq_file in files:
        file_path = os.path.abspath(os.path.join(directory, seq_file))

        if metadata_type_value == "SR" and ("R1" in seq_file or "R2" in seq_file):
            sample_id = extract_sample_id(seq_file, sr_pattern)
            file_type = "R1" if "R1" in seq_file else "R2"
            if sample_id:
                metadata.append([sample_id, file_type, file_path])

        elif metadata_type_value == "LR" and "nanopore" in seq_file and f".{lr_format}" in seq_file:
            sample_id = extract_sample_id(seq_file, lr_pattern)
            if sample_id:
                metadata.append([sample_id, "long", file_path])

        elif metadata_type_value == "all":
            sample_id = None
            file_type = None
            if "R1" in seq_file or "R2" in seq_file:
                sample_id = extract_sample_id(seq_file, sr_pattern)
                file_type = "R1" if "R1" in seq_file else "R2"
            elif "nanopore" in seq_file and f".{lr_format}" in seq_file:
                sample_id = extract_sample_id(seq_file, lr_pattern)
                file_type = "long"
            if sample_id and file_type:
                metadata.append([sample_id, file_type, file_path])

    return pd.DataFrame(metadata, columns=["sample_id", "type", "sample"])


def write_metadata_tsv(df: pd.DataFrame, output_file: Path | str) -> None:
    """Write metadata DataFrame in expected TSV format (`sample`, `sample_id`, `type`)."""
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    out = df[["sample", "sample_id", "type"]].copy()
    out.to_csv(output_file, sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate metadata TSV from simulated reads directory.")
    parser.add_argument("directory", type=str)
    parser.add_argument("metadata_type", choices=["SR", "LR", "all"])
    parser.add_argument("--lr-format", default="fastq", choices=["fastq", "fasta"])
    parser.add_argument("--output", default="metadata.tsv")
    args = parser.parse_args(argv)

    df = generate_metadata(args.directory, args.metadata_type, lr_format=args.lr_format)
    write_metadata_tsv(df, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
