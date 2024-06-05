ASSEMBLER_LR = config['assembly']['long_read_assembler'] 

# taking into account the case where we don't have LR
if ASSEMBLER_LR == None:
       ASSEMBLER_LR = []

BINNER = config['binning']['long_read_binner'] 

SAMPLES_TABLE = config['samples']
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)

rule binette_refinement_LR:
    input: 
        bins_dirs = expand("results/05_binning/LR/{binner}/bins/{{assembler}}/{{sample}}", binner=BINNER),
        assembly = "results/03_assembly/LR/{assembler}/{sample}/assembly.fa.gz",
        checkm2_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output: directory("results/07_bins_refinement/binette/{assembler}/{sample}") # bins will be in a "final_bins" folder
    conda:
        "../envs/binette.yaml"
    log:
        stdout = "logs/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.stdout",
        stderr = "logs/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.stderr",
        stdout_check = "logs/07_bins_refinement/binette/{assembler}/{sample}/check.stdout",
        stderr_check = "logs/07_bins_refinement/binette/{assembler}/{sample}/check.stderr"
    params:
        # constructing the precise folder path with bins from the input
        bins_folder = lambda wildcards, input: [f"{dir}/bins" for dir in input.bins_dirs]
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER_LR)
    shell:
        """
        binette --bin_dirs {params.bins_folder} --contigs {input.assembly} \
        --checkm2_db {input.checkm2_database} --outdir {output} \
        > {log.stdout} 2> {log.stderr} \
        && \
        bash workflow/scripts/check_binette_produced_files.sh {output}/final_bins {input.bins_dirs} \
        > {log.stdout_check} 2> {log.stderr_check}
        """