# no matter the assembler used next is the assembly name to obtain
assembly_name = "{sample}.final_assembly.fa"

rule megahit_assembly:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        touch("results/03_assembly/{sample}_assembly"),
        assembly = f'results/03_assembly/assembly/{assembly_name}'
    conda:
        "../envs/megahit.yaml"
    log:
        stdout = "logs/03_assembly/megahit/{sample}.stdout",
        stderr = "logs/03_assembly/megahit/{sample}.stdout"
    params:
        threads = config['megahit']['threads'],
        tmp_dir = "tmp/"
    shell:
        """
        mkdir -p {params.tmp_dir} \
        && \
        megahit -1 {input.r1} -2 {input.r2} \
            --num-cpu-threads {params.threads} \
            --tmp-dir {params.tmp_dir} \
            --out-dir results/03_assembly/megahit > {log.stdout} 2> {log.stderr} \
        && \
        mv results/03_assembly/megahit/final.contigs.fa {output.assembly}
        """
