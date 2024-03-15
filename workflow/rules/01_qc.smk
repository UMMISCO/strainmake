rule fastqc_before_preprocessing:
    input:
        "data/{sample}_{read}.fastq.gz",
    output:
        html_report="results/01_qc/fastqc/{sample}_{read}_fastqc.html",
        zip_report="results/01_qc/fastqc/{sample}_{read}_fastqc.zip"
    conda: 
        "../envs/fastqc.yaml"
    log:
        stdout = "logs/01_qc/fastqc/{sample}_{read}.stdout",
        stderr = "logs/01_qc/fastqc/{sample}_{read}.stderr"
    shell:
        "fastqc {input} -o results/01_qc/fastqc/ > {log.stdout} 2> {log.stderr}"