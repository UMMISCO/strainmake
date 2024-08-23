"""
Automatically generates TSV table with path to FASTQ for use in pipeline

./generate_metadata_table.py --help
"""

import os
import typer
import re
import pandas as pd
from enum import Enum

app = typer.Typer()

def extract_sample_id(filename: str, pattern: str) -> str:
    match = re.search(pattern, filename)
    if match:
        return match.group(1)
    return None

class MetadataType(str, Enum):
    SR = "SR"
    LR = "LR"
    ALL = "all"

@app.command()
def generate_metadata(directory: str, metadata_type: MetadataType):
    """
    Generate a TSV file listing .fastq.gz files in the directory with metadata
    
    Args:
    - directory: The path to the directory containing .fastq.gz files
    - metadata_type: The type of metadata to filter by (SR, LR, or all)
    """
    # prepare patterns for sample ID extraction
    sr_pattern = r"fake_illumina_R[12]\.(.*)\.fastq\.gz"
    lr_pattern = r"fake_nanopore_.*\.(.*)\.fastq\.gz"

    # list all .fastq.gz files in the directory
    fastq_files = [f for f in os.listdir(directory) if f.endswith(".fastq.gz")]

    metadata = []

    # process files based on metadata_type
    for fastq_file in fastq_files:
        file_path = os.path.abspath(os.path.join(directory, fastq_file))
        if metadata_type == "SR" and ("R1" in fastq_file or "R2" in fastq_file):
            sample_id = extract_sample_id(fastq_file, sr_pattern)
            file_type = "R1" if "R1" in fastq_file else "R2"
            metadata.append([sample_id, file_type, file_path])

        elif metadata_type == "LR" and "nanopore" in fastq_file:
            sample_id = extract_sample_id(fastq_file, lr_pattern)
            metadata.append([sample_id, "long", file_path])

        elif metadata_type == "all":
            if "R1" in fastq_file or "R2" in fastq_file:
                sample_id = extract_sample_id(fastq_file, sr_pattern)
                file_type = "R1" if "R1" in fastq_file else "R2"
            elif "nanopore" in fastq_file:
                sample_id = extract_sample_id(fastq_file, lr_pattern)
                file_type = "long"
            else:
                continue
            metadata.append([sample_id, file_type, file_path])

    # create a DataFrame and write to a TSV file
    metadata_df = pd.DataFrame(metadata, columns=["sample_id", "type", "sample"])
    output_file = os.path.join("metadata.tsv")
    metadata_df.to_csv(output_file, sep="\t", index=False)

    typer.echo(f"Metadata file generated: {output_file}")

if __name__ == "__main__":
    app()
