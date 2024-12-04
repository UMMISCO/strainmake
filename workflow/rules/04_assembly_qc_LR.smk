rule quast_qc_long_read:
    input:
        # long read assemblies produced in step 03
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz"
    output:
        "results/04_assembly_qc/quast_long_read/{assembler_lr}/{sample_lr}/report.html"
    conda:
        "../envs/quast.yaml"
    log:
        stdout = "logs/04_assembly_qc/quast/{assembler_lr}/{sample_lr}.stdout",
        stderr = "logs/04_assembly_qc/quast/{assembler_lr}/{sample_lr}.stdout"
    benchmark:
        "benchmarks/04_assembly_qc/quast/{assembler_lr}/{sample_lr}.benchmark.txt"
    params:
        out_dir = "results/04_assembly_qc/quast_long_read/{assembler_lr}/{sample_lr}",
        method = "--nanopore" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio"
    threads: config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            --circos \
            {input.assembly} \
            > {log.stdout} 2> {log.stderr} 
        """