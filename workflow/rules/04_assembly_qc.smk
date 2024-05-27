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
    threads: threads = config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            {input} \
            > {log.stdout} 2> {log.stderr} 
        """