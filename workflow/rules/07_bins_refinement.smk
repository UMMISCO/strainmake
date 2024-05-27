from utils import *

ASSEMBLER = config['assembly']['assembler'] 
BINNER = config['binning']['binner'] 

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)

rule binette_refinement:
    input: 
        bins_dirs = expand("results/05_binning/{binner}/bins/{{assembler}}/{{sample}}", binner=BINNER),
        assembly = "results/03_assembly/{assembler}/{sample}/assembly.fa.gz",
        checkm2_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output: directory("results/07_bins_refinement/binette/{assembler}/{sample}") # bins will be in a "final_bins" folder
    conda:
        "../envs/binette.yaml"
    log:
        stdout = "logs/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.stdout",
        stderr = "logs/07_bins_refinement/binette/{assembler}/{sample}/bins_refinement.stderr"
    params:
        # constructing the precise folder path with bins from the input
        bins_folder = lambda wildcards, input: [f"{dir}/bins" for dir in input.bins_dirs]
    shell:
        """
        binette --bin_dirs {params.bins_folder} --contigs {input.assembly} \
        --checkm2_db {input.checkm2_database} --outdir {output} \
        > {log.stdout} 2> {log.stderr}
        """