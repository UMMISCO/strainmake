# no matter the assembler used next is the assembly name to obtain
assembly_name = "{sample}.final_assembly.fasta"


rule quast_qc:
    input:
        # assembly produced in step 03
        f'results/03_assembly/assembly/{assembly_name}'
    output:
        touch("results/04_assembly_qc/{sample}_assembly_qc")
    conda:
        "../envs/quast.yaml"
    log:
        stdout = "logs/04_assembly_qc/quast/{sample}.stdout",
        stderr = "logs/04_assembly_qc/quast/{sample}.stdout"
    params:
        threads = config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o results/04_assembly_qc/quast \
            --max-ref-number 0 \
            {input} \
            > {log.stdout} 2> {log.stderr} 
        """