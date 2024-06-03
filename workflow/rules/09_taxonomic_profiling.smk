rule metaphlan_profiling:
    input: "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz"
    output: "results/09_taxonomic_profiling/metaphlan/{sample}.profile.txt"
    conda:
        "../envs/metaphlan.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/metaphlan/{sample}.profile.stdout",
        stderr = "logs/09_taxonomic_profiling/metaphlan/{sample}.profile.stdout"
    threads: config['taxonomic_profiling']['metaphlan']['threads']
    shell:
        """
        metaphlan --input_type fastq --nproc {threads} \
            {input} {output} \
        > {log.stdout} 2> {log.stderr}
        """