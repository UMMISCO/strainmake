rule quast_qc_long_read:
    input:
        # long read assemblies produced in step 03
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz",
        long_reads = "results/02_preprocess/fastp_long_read/{sample_lr}_1.fastq.gz"
    output:
        "results/04_assembly_qc/quast_long_read/{assembler_lr}/{sample_lr}/report.html"
    conda:
        "../envs/quast.yaml"
    log:
        stdout = "logs/04_assembly_qc/quast/{assembler_lr}/{sample_lr}.stdout",
        stderr = "logs/04_assembly_qc/quast/{assembler_lr}/{sample_lr}.stdout"
    params:
        out_dir = "results/04_assembly_qc/quast_long_read/{assembler_lr}/{sample_lr}",
        method = "--nanopore" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio"
    threads: config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            --circos \
            {params.method} {input.long_reads} \
            {input.assembly} \
            > {log.stdout} 2> {log.stderr} 
        """