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