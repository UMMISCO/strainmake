SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)

rule gtdb_tk_taxonomic_annotation:
    input:
        # folder with refined bins
        refined_bins = "results/07_bins_refinement/binette/{assembler}/{sample}",
        ref_data = "data/gtdb_tk/release220"
    output: directory("results/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}")
    conda:
        "../envs/gtdb_tk.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}/classify.stdout",
        stderr = "logs/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}/classify.stderr"
    threads: config['bins_postprocessing']['gtdbtk']['threads']
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.refined_bins}/final_bins --cpus {threads} --out_dir {output} \
            --extension ".fa" \
            --skip_ani_screen --pplacer_cpus 1 \
            > {log.stdout} 2> {log.stderr}
        """

# this rule produces text files with list of refined genomes for use
# in dRep
rule list_refined_genomes:
    input: expand("results/07_bins_refinement/binette/{{assembler}}/{sample}", sample=SAMPLES)
    output: "results/08_bins_postprocessing/genomes_list/{assembler}/list.txt"
    log:
        stdout = "logs/08_bins_postprocessing/genomes_list/{assembler}/list.stdout",
        stderr = "logs/08_bins_postprocessing/genomes_list/{assembler}/list.stderr"
    shell:
        """
        bash workflow/scripts/list_refined_genomes.sh {output} {input}/final_bins \
            > {log.stdout} 2> {log.stderr}
        """


rule genomes_dereplication:
    # input is formed of every refined produced no matter the sample
    # it is, however, assembler specific
    input: "results/08_bins_postprocessing/genomes_list/{assembler}/list.txt"
    output: directory("results/08_bins_postprocessing/dRep/{assembler}")
    conda:
        "../envs/drep.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/drep/{assembler}/dereplication.stdout",
        stderr = "logs/08_bins_postprocessing/drep/{assembler}/dereplication.stderr"
    params:
        comparison_algorithm = config['bins_postprocessing']['drep']['comparison_algorithm'],
        other_args = config['bins_postprocessing']['drep']['other_args'],
    threads: config['bins_postprocessing']['drep']['threads']
    shell:
        """
        dRep dereplicate --genomes {input} --processors {threads} \
            --S_algorithm {params.comparison_algorithm} \
            {params.other_args} \
            {output} \
        > {log.stdout} 2> {log.stderr}
        """