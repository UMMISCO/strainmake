rule quast_qc:
    input:
        # assemblies produced in step 03
        "results/03_assembly/{assembler}/{sample}/assembly.fa.gz"
    output:
        "results/04_assembly_qc/quast/{assembler}/{sample}/report.html"
    conda:
        "../envs/quast.yaml"
    log:
        stdout = "logs/04_assembly_qc/quast/{assembler}/{sample}.stdout",
        stderr = "logs/04_assembly_qc/quast/{assembler}/{sample}.stdout"
    params:
        out_dir = "results/04_assembly_qc/quast/{assembler}/{sample}"
    threads: config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            --circos \
            {input} \
            > {log.stdout} 2> {log.stderr} 
        """

rule quast_qc_long_read:
    input:
        # long read assemblies produced in step 03
        assembly = "results/03_assembly/{long_read_assembler}/{sample}/assembly.fa.gz",
        long_reads = "results/02_preprocess/fastp_long_read/{sample}_1.fastq.gz"
    output:
        "results/04_assembly_qc/quast_long_read/{long_read_assembler}/{sample}/report.html"
    conda:
        "../envs/quast.yaml"
    log:
        stdout = "logs/04_assembly_qc/quast/{long_read_assembler}/{sample}.stdout",
        stderr = "logs/04_assembly_qc/quast/{long_read_assembler}/{sample}.stdout"
    params:
        out_dir = "results/04_assembly_qc/quast_long_read/{long_read_assembler}/{sample}",
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