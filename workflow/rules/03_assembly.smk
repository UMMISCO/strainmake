rule megahit_assembly:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        touch("results/03_assembly/{sample}_assembly")
    conda:
        "../envs/megahit.yaml"
    log:
        stdout = "logs/03_assembly/megahit/{sample}.stdout",
        stderr = "logs/03_assembly/megahit/{sample}.stdout"
    shell:
        """
        mkdir -p results/03_assembly/megahit/tmp \
        && \
        megahit -1 {input.r1} -2 {input.r2} \
            --num-cpu-threads {MEGAHIT_THREADS} \
            --tmp-dir results/03_assembly/megahit/tmp
            --out-dir results/03_assembly/megahit > {log.stdout} 2> {log.stderr}
        """
