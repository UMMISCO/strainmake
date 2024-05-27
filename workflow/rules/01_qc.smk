import pandas as pd 

SAMPLES_TABLE = config['samples']
SAMPLES_DF = pd.read_csv(SAMPLES_TABLE, sep="\t")
FASTQ_FILES = SAMPLES_DF['sample'].tolist()
FASTQ_FILES = [f[:-9] for f in FASTQ_FILES]

rule fastqc_before_preprocessing:
    input: expand("{fastq}.fastq.gz", fastq=FASTQ_FILES)
    output: directory("results/01_qc/fastqc/{fastq_files}/")
    conda: 
        "../envs/fastqc.yaml"
    log:
        stdout = "logs/01_qc/fastqc/{fastq_files}.stdout",
        stderr = "logs/01_qc/fastqc/{fastq_files}.stderr"   
    wildcard_constraints:
        fastq_files="|".join(FASTQ_FILES)
    shell:
        """
        mkdir results/01_qc/fastqc/{wildcards.fastq_files} \
        && \
        fastqc {input} -o results/01_qc/fastqc/{wildcards.fastq_files} > {log.stdout} 2> {log.stderr}
        """