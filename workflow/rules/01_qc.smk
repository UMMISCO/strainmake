include: "../Snakefile"

import pandas as pd

SAMPLES_TABLE = config['samples']
SAMPLES_DF = pd.read_csv(SAMPLES_DF, sep="\t")
FASTQ_FILES = SAMPLES_DF['sample'].tolist()

rule fastqc_before_preprocessing:
    input:
        expand("{FASTQ_FILES}", fastq_files=FASTQ_FILES)
    output: directory("results/01_qc/fastqc/")
    conda: 
        "../envs/fastqc.yaml"
    log:
        stdout = "logs/01_qc/fastqc/{sample}_{read}.stdout",
        stderr = "logs/01_qc/fastqc/{sample}_{read}.stderr"
    shell:
        "fastqc {input} -o results/01_qc/fastqc/ > {log.stdout} 2> {log.stderr}"