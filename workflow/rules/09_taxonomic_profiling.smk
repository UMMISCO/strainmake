import os

rule metaphlan_profiling:
    input: "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz" # on short reads only
    output: "results/09_taxonomic_profiling/metaphlan/{sample}.profile.txt"
    conda:
        "../envs/metaphlan.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/metaphlan/{sample}.profile.stdout",
        stderr = "logs/09_taxonomic_profiling/metaphlan/{sample}.profile.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/metaphlan/{sample}.profile.benchmark.txt"
    threads: config['taxonomic_profiling']['metaphlan']['threads']
    shell:
        """
        metaphlan --input_type fastq --nproc {threads} \
            {input} {output} \
        > {log.stdout} 2> {log.stderr}
        """

# we create a folder with the short reads METEOR will work on (we use symbolic links)
rule meteor_prepare_fastq:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        # same fastq files, but in another folder
        r1 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_1.fastq.gz",
        r2 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_2.fastq.gz"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_prepare_fastq_folder.stdout",
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.prepare_fastq_folder.benchmark.txt"
    shell:
        """
        ln -sv $(realpath {input.r1}) {output.r1} > {log} 2>&1
        ln -sv $(realpath {input.r2}) {output.r2} >> {log} 2>&1
        """

rule meteor_fastq_indexing:
    input:
        r1 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_1.fastq.gz",
        r2 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_2.fastq.gz"
    output:
        # same fastq files, but in another folder
        directory("results/09_taxonomic_profiling/meteor/{sample}/fastq_index")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_fastq_indexing.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_fastq_indexing.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.fastq_indexing.benchmark.txt"
    params:
        fastq_data = lambda wildcards: os.path.dirname("results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_1.fastq.gz".format(sample=wildcards.sample))
    shell:
        """
        meteor fastq -i $(realpath {params.fastq_data}) -p -o {output} > {log.stdout} 2> {log.stderr}
        """

rule meteor_mapping:
    input:
       "results/09_taxonomic_profiling/meteor/{sample}/fastq_index"
    output:
        directory("results/09_taxonomic_profiling/meteor/{sample}/mapping")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_mapping.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_mapping.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.fastq_mapping.benchmark.txt"
    threads:
        config['taxonomic_profiling']['meteor']['threads']
    params:
        indexed_fastq_file_with_sample = lambda wildcards: os.path.join(f"results/09_taxonomic_profiling/meteor/{wildcards.sample}/fastq_index", f"{wildcards.sample}")
    shell:
        """
        meteor mapping -i {params.indexed_fastq_file_with_sample} -o {output} -r $REFERENCE \
            --trim 0 --ka --kf -t {threads} \
            > {log.stdout} 2> {log.stderr} 
        """

rule meteor_profiling:
    input:
        "results/09_taxonomic_profiling/meteor/{sample}/mapping"
    output:
        directory("results/09_taxonomic_profiling/meteor/{sample}/profiling")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_profiling.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_profiling.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.profile.benchmark.txt"
    shell:
        """ 
        meteor profile -i {input} -o {output} -r $REFERENCE \
            -n coverage \
            > {log.stdout} 2> {log.stderr}
        """