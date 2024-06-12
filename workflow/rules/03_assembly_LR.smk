# metaFlye for long read
rule metaflye_assembly:
    input:
        # files produced by fastp (should have already been decontaminated before used in the pipeline)
        long_read = "results/02_preprocess/fastp_long_read/{sample}_1.fastq.gz",
    output:
        assembly = "results/03_assembly/LR/metaflye/{sample}/assembly.fa.gz"
    conda:
        "../envs/flye.yaml"
    log:
        stdout = "logs/03_assembly/metaflye/{sample}.stdout",
        stderr = "logs/03_assembly/metaflye/{sample}.stderr"
    params:
        out_dir = "results/03_assembly/LR/metaflye/{sample}",
        method_flag = "--nano-hq" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio-hifi"
    threads: config['assembly']['metaflye']['threads']
    shell:
        """
        flye {params.method_flag} {input.long_read} --out-dir {params.out_dir} \
            --meta --threads {threads} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.out_dir}/assembly.fasta {params.out_dir}/assembly.fa \
        && \
        pigz {params.out_dir}/assembly.fa
        """

# hybrid assembly using hybridSPADes
rule hybridspades_assembly:
    input:
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz",
        long_read = "results/02_preprocess/fastp_long_read/{sample}_1.fastq.gz"
    output:
        assembly = "results/03_assembly/LR/hybridspades/{sample}/assembly.fa.gz"
    conda:
        "../envs/spades.yaml"
    log:
        stdout = "logs/03_assembly/hybridspades/{sample}.stdout",
        stderr = "logs/03_assembly/hybridspades/{sample}.stdout"
    params:
        out_dir = "results/03_assembly/LR/hybridspades/{sample}",
        memory_limit = config['assembly']['hybridspades']['memory_limit'],
        method_flag = "--nanopore" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio"
    threads: config['assembly']['hybridspades']['threads']
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} \
            {params.method_flag} {input.long_read} \
            --threads {threads} \
            -o {params.out_dir} \
            -m {params.memory_limit} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.out_dir}/scaffolds.fasta {params.out_dir}/assembly.fa \
        && \
        pigz {params.out_dir}/assembly.fa
        """
