ASSEMBLER = config['assembly']['assembler']
HYBRID_ASSEMBLER = config['assembly']['hybrid_assembler'] 
ASSEMBLER_LR = config['assembly']['long_read_assembler']

# taking into account the case where we don't have SR assembly
if ASSEMBLER == None:
       ASSEMBLER = []

# taking into account the case where we don't have hybrid assembly
if HYBRID_ASSEMBLER == None:
       HYBRID_ASSEMBLER = []

# taking into account the case where we don't have LR
if ASSEMBLER_LR == None:
       ASSEMBLER_LR = []

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
        stderr = "logs/04_assembly_qc/quast/{assembler}/{sample}.stderr"
    benchmark:
        "benchmarks/04_assembly_qc/quast/{assembler}/{sample}.benchmark.txt"
    params:
        out_dir = "results/04_assembly_qc/quast/{assembler}/{sample}"
    threads: config['quast']['threads']
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
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
        assembly = "results/03_assembly/{assembler_lr}/{sample_lr}/assembly.fa.gz"
    output:
        "results/04_assembly_qc/quast/{assembler_lr}/{sample_lr}/report.html"
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
    wildcard_constraints:
        assembler_lr = "|".join(ASSEMBLER_LR)
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            --circos \
            {input.assembly} \
            > {log.stdout} 2> {log.stderr} 
        """