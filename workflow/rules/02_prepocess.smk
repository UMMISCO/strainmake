rule fastp:
    input:
       r1 = "data/{sample}_1.fastq.gz",
       r2 = "data/{sample}_2.fastq.gz"
    output:
        r1 = "results/02_preprocess/fastp/{sample}_1.fastq.gz",
        r2 = "results/02_preprocess/fastp/{sample}_2.fastq.gz"
    conda: 
        "../envs/fastp.yaml"
    log:
        stdout = "logs/02_preprocess/fastp/{sample}.stdout",
        stderr = "logs/02_preprocess/fastp/{sample}.stdout"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            --length_required {MIN_READ_LENGTH} \
            --qualified_quality_phred {MIN_PHRED} \
            > {log.stdout} 2> {log.stderr}
        """

rule fastqc_after_preprocessing:
    input:
        r1 = "results/02_preprocess/fastp/{sample}_{read}.fastq.gz"
    output:
        html_report="results/02_preprocess/fastqc/{sample}_{read}_fastqc.html",
        zip_report="results/02_preprocess/fastqc/{sample}_{read}_fastqc.zip"
    conda: 
        "../envs/fastqc.yaml"
    log:
        stdout = "logs/02_preprocess/fastqc/{sample}_{read}.stdout",
        stderr = "logs/02_preprocess/fastqc/{sample}_{read}.stderr"
    shell:
        "fastqc {input} -o results/02_preprocess/fastqc/ > {log.stdout} 2> {log.stderr}"
