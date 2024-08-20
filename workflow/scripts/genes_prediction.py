"""
A CLI to predict genes in genomes stored in FASTA or gzipped FASTA in a folder using
Prodigal
"""

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import argparse


def run_prodigal(input_file, output_file):
    """
    Run Prodigal on an input genome (`input_file`) and write the prediction results 
    into `output_file
    """
    cmd = ['prodigal', '-i', input_file, '-d', output_file]
    subprocess.run(cmd, check=True)

def process_file(file_path, output_dir):
    """
    Calling genes in a genome stored in .fa or .fa.gz
    """
    file_name = os.path.basename(file_path)
    output_file = os.path.join(output_dir, file_name.replace('.fa', '.genes.fna').replace('.gz', ''))

    if file_path.endswith('.gz'):
        # in the case of a gzipped FASTA file, it temporarily un-gzip it to 
        # process it with Prodigal
        temporarily_file = output_file.replace('.genes.fna', '')
        with gzip.open(file_path, 'rb') as f_in:
            with open(temporarily_file, 'wb') as f_out:
                f_out.write(f_in.read())

        # running Prodigal using the un-gzipped file
        run_prodigal(temporarily_file, output_file)
        # removing the un-gzipped file
        os.remove(temporarily_file)
    else:
        run_prodigal(file_path, output_file)

def main(input_dir, output_dir, cpu):
    """
    CLI logic. Allows to run genes calling in parallel
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # listing all .fa and .fa.gz files stored in the input directory
    input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fa') or f.endswith('.fa.gz')]

    # then, using ThreadPoolExecutor, it runs in parallel, genes prediction on the genomes 
    with ThreadPoolExecutor(max_workers=cpu) as executor:
        future_to_file = {executor.submit(process_file, file_path, output_dir): file_path for file_path in input_files}
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                future.result()
                print(f"Traitement terminé pour {file_path}")
            except Exception as exc:
                print(f"Erreur lors du traitement de {file_path}: {exc}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to run Prodigal in parallel on .fa or .fa.gz files")
    parser.add_argument("input", help="Path to folder containing genomes")
    parser.add_argument("output", help="Path to folder that will store output")
    parser.add_argument("--cpu", type=int, default=1, help="CPUs to use")

    args = parser.parse_args()

    main(args.input, args.output, args.cpu)
