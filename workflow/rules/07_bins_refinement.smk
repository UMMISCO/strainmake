from utils import *

ASSEMBLER = config['assembly']['assembler'] 

HYBRID_ASSEMBLER = config['assembly']['hybrid_assembler'] 
ASSEMBLER_LR = config['assembly']['long_read_assembler'] 

# taking into account the case where we don't have SR
if ASSEMBLER == None:
       ASSEMBLER = []

# taking into account the case where we don't have LR
if HYBRID_ASSEMBLER == None:
       HYBRID_ASSEMBLER = []

# taking into account the case where we don't have LR
if ASSEMBLER_LR == None:
       ASSEMBLER_LR = []

SHORT_READ_BINNER = config['binning']['binner']
LONG_READ_BINNER = config['binning']['long_read_binner']

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)

rule binette_refinement:
    input: 
        bins_dirs = expand("results/05_binning/{binner}/bins/{{assembler}}/{{sample}}", binner=SHORT_READ_BINNER),
        assembly = "results/03_assembly/{assembler}/{sample}/assembly.fa.gz",
        checkm2_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output: directory("results/07_bins_refinement/binette/{assembler}/{sample}") # bins will be in a "final_bins" folder
    conda:
        "../envs/binette.yaml"
    log:
        stdout = "logs/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.stdout",
        stderr = "logs/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.stderr",
        stdout_check = "logs/07_bins_refinement/binette/{assembler}/{sample}/check.stdout",
        stderr_check = "logs/07_bins_refinement/binette/{assembler}/{sample}/check.stderr"
    benchmark:
        "benchmarks/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.benchmark.txt"
    params:
        # constructing the precise folder path with bins from the input
        bins_folder = lambda wildcards, input: [f"{dir}/bins" for dir in input.bins_dirs],
        low_mem = config["bins_refinement"]["binette"]["low_mem"]
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    threads: config["bins_refinement"]["binette"]["threads"]
    shell:
        """
        binette --bin_dirs {params.bins_folder} --contigs {input.assembly} \
        --checkm2_db {input.checkm2_database} --outdir {output} \
        --threads {threads} {params.low_mem} --verbose \
        > {log.stdout} 2> {log.stderr} \
        && \
        bash workflow/scripts/check_binette_produced_files.sh {output}/final_bins {input.bins_dirs} \
        > {log.stdout_check} 2> {log.stderr_check}
        """

rule binette_refinement_LR:
    input: 
        bins_dirs = expand("results/05_binning/{binner_lr}/bins/{{assembler_lr}}/{{sample_lr}}", binner_lr = LONG_READ_BINNER),
        assembly = "results/03_assembly/{assembler_lr}/{sample_lr}/assembly.fa.gz",
        checkm2_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output: directory("results/07_bins_refinement/binette/{assembler_lr}/{sample_lr}") # bins will be in a "final_bins" folder
    conda:
        "../envs/binette.yaml"
    log:
        stdout = "logs/07_bins_refinement/binette/{assembler_lr}/{sample_lr}/bins_refinement.stdout",
        stderr = "logs/07_bins_refinement/binette/{assembler_lr}/{sample_lr}/bins_refinement.stderr",
        stdout_check = "logs/07_bins_refinement/binette/{assembler_lr}/{sample_lr}/check.stdout",
        stderr_check = "logs/07_bins_refinement/binette/{assembler_lr}/{sample_lr}/check.stderr"
    benchmark:
        "benchmarks/07_bins_refinement/binette/{assembler_lr}/{sample_lr}/bins_refinement.benchmark.txt"
    params:
        # constructing the precise folder path with bins from the input
        bins_folder = lambda wildcards, input: [f"{dir}/bins" for dir in input.bins_dirs],
        low_mem = config["bins_refinement"]["binette"]["low_mem"]
    wildcard_constraints:
        sample_lr = "|".join(SAMPLES_LR),
        assembler_lr = "|".join(ASSEMBLER_LR)
    threads: config["bins_refinement"]["binette"]["threads"]
    shell:
        """
        binette --bin_dirs {params.bins_folder} --contigs {input.assembly} \
        --checkm2_db {input.checkm2_database} --outdir {output} \
        --threads {threads} {params.low_mem} --verbose \
        > {log.stdout} 2> {log.stderr} \
        && \
        bash workflow/scripts/check_binette_produced_files.sh {output}/final_bins {input.bins_dirs} \
        > {log.stdout_check} 2> {log.stderr_check}
        """