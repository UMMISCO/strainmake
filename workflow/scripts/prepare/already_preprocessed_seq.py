"""
This script sets up a `results` directory for sequencing data that has already been preprocessed. 

It reads a metadata table in the pipeline's accepted format, which contains paths to the FASTQ files. 
The script then creates symbolic links in the `results` folder pointing to these original FASTQ files, 
rather than copying them. This helps save storage space.

After running this script successfully, you can execute the pipeline with a command like:
    `snakemake -c 15 -p --conda-frontend conda --use-conda -s workflow/Snakefile.already_preprocessed_seq.smk`

This will start the pipeline from the assembly steps.

If you want to remove the symbolic links created by this script, use the `--clean` option.
"""

import argparse
import pandas as pd
from pathlib import Path

def get_SR_samples(metadata_table: pd.DataFrame):
    """ 
    Returns a dictionary with one key by sample (column `sample_id`), and 
    for each sample another dictionary with path to the R1 FASTQ and path
    to the R2 FASTQ (info in column `type`).

    If there are no short reads samples in the metadata, it will return
    None.
    """
    sr_samples = {}
    sr_data = metadata_table[metadata_table['type'].isin(['R1', 'R2'])]

    if sr_data.empty:
        return None

    for sample_id, group in sr_data.groupby('sample_id'):
        sr_samples[sample_id] = {
            'R1': group[group['type'] == 'R1']['sample'].values[0],
            'R2': group[group['type'] == 'R2']['sample'].values[0],
        }

    return sr_samples

def get_LR_samples(metadata_table: pd.DataFrame):
    """ 
    Returns a dictionary with one key by sample (column `sample_id`), and 
    for each sample the path to the long read FASTQ.

    If there are no long reads samples in the metadata, it will return
    None.
    """
    lr_samples = {}
    lr_data = metadata_table[metadata_table['type'] == 'long']

    if lr_data.empty:
        return None

    for sample_id, group in lr_data.groupby('sample_id'):
        lr_samples[sample_id] = group['sample'].values[0]

    return lr_samples

def prepare_results_dir(
    SR_samples,  # type "dict" or None if "get_SR_samples" returned None
    LR_samples,  # type "dict" or None if "get_LR_samples" returned None
    results_dir  # path where the results folder will be created
):
    """ 
    Create a `results` folder (fail if it already exists).
    
    For short reads samples, it creates `results/02_preprocess/bowtie2` folder 
    with inside files named `{sample}_1.clean.fastq.gz` symlinked to 
    the real FASTQ and `{sample}_2.clean.fastq.gz` symlinked to 
    the real FASTQ.

    For long reads samples, it creates `results/02_preprocess/fastp_long_read`
    folder with inside files named `{sample}_1.fastq.gz` symlinked to the real FASTQ.
    """

    # create results directory structure
    bowtie2_dir = results_dir / "02_preprocess" / "bowtie2"
    fastp_long_read_dir = results_dir / "02_preprocess" / "fastp_long_read"

    bowtie2_dir.mkdir(parents=True, exist_ok=False)
    fastp_long_read_dir.mkdir(parents=True, exist_ok=False)

    # handle short reads
    if SR_samples:
        for sample, paths in SR_samples.items():
            r1_target = bowtie2_dir / f"{sample}_1.clean.fastq.gz"
            r2_target = bowtie2_dir / f"{sample}_2.clean.fastq.gz"

            # convert paths['R1'] and paths['R2'] to absolute paths
            r1_source = Path(paths['R1']).resolve()
            r2_source = Path(paths['R2']).resolve()

            r1_target.symlink_to(r1_source)
            r2_target.symlink_to(r2_source)

            print(f"Symlinked {r1_source} -> {r1_target}")
            print(f"Symlinked {r2_source} -> {r2_target}")

    # handle long reads
    if LR_samples:
        for sample, path in LR_samples.items():
            lr_target = fastp_long_read_dir / f"{sample}_1.fastq.gz"
            
            # convert path to an absolute path
            lr_source = Path(path).resolve()

            lr_target.symlink_to(lr_source)
            print(f"Symlinked {lr_source} -> {lr_target}")

def clean_results_dir(results_dir):
    """
    Safely remove the symbolic links created in the `results` folder.
    """
    # paths to clean
    sr_dir = results_dir / "02_preprocess" / "bowtie2"
    lr_dir = results_dir / "02_preprocess" / "fastp_long_read"

    # remove symlinks in short reads directory
    if sr_dir.exists():
        for file in sr_dir.iterdir():
            if file.is_symlink():
                file.unlink()
                print(f"Removed symlink: {file}")
        print(f"Cleaned symbolic links in {sr_dir}")

    # remove symlinks in long reads directory
    if lr_dir.exists():
        for file in lr_dir.iterdir():
            if file.is_symlink():
                file.unlink()
                print(f"Removed symlink: {file}")
        print(f"Cleaned symbolic links in {lr_dir}")

    print("Symlinks removed successfully.")

def main():
    parser = argparse.ArgumentParser(description="Prepare results folder with symbolic links for FASTQ files.")
    parser.add_argument("tsv_file", help="path to the input TSV file with sample metadata.")
    parser.add_argument("--clean", action="store_true", help="remove symbolic links from the results folder.")
    parser.add_argument("--results_dir", default="results", help="directory where the 'results' folder will be created (default is './results').")
    args = parser.parse_args()

    # convert results_dir to path object
    results_dir = Path(args.results_dir).resolve()

    if args.clean:
        clean_results_dir(results_dir)
    else:
        # load the tsv file
        metadata_table = pd.read_csv(args.tsv_file, sep='\t')

        # get short reads and long reads samples
        SR_samples = get_SR_samples(metadata_table)
        LR_samples = get_LR_samples(metadata_table)

        # prepare the results directory and create symlinks
        prepare_results_dir(SR_samples, LR_samples, results_dir)

if __name__ == "__main__":
    main()
