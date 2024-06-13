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

# this rule produces text files with list of refined genomes
rule list_refined_genomes:
    input: expand("results/07_bins_refinement/binette/{{assembler}}/{sample}", sample=SAMPLES)
    output: "results/08_bins_postprocessing/genomes_list/{assembler}/list.txt"
    log:
        stdout = "logs/08_bins_postprocessing/genomes_list/{assembler}/list.stdout",
        stderr = "logs/08_bins_postprocessing/genomes_list/{assembler}/list.stderr"
    params:
        # constructing the precise folder path with bins from the input
        bins_folder = lambda wildcards, input: [f"{dir}/final_bins" for dir in input]
    shell:
        """
        bash workflow/scripts/list_refined_genomes.sh {output} {params.bins_folder} \
            > {log.stdout} 2> {log.stderr}
        """

# copying bins into another foler and renaming them if needed to avoid duplicated names 
# (what makes dRep fail)
# the copied bins will be removed once dRep is done
rule copy_and_rename_bins:
    input:
        "results/08_bins_postprocessing/genomes_list/{assembler}/list.txt"
    output:
        bins_dest = directory("results/08_bins_postprocessing/genomes_list/genomes/{assembler}"),
        bins_name_link_table = "results/08_bins_postprocessing/genomes_list/{assembler}/unduplicated.tsv"
    log:
        stdout = "logs/08_bins_postprocessing/genomes_list/{assembler}/copy_and_rename.stdout",
        stderr = "logs/08_bins_postprocessing/genomes_list/{assembler}/copy_and_rename.stderr"
    shell:
        """
        python3 workflow/scripts/make_bin_names_unambiguous.py {input} {output.bins_name_link_table} \
            {output.bins_dest} > {log.stdout} 2> {log.stderr}
        """

# this rule produces text files with list of refined genomes that have been copied and renamed
# to avoid bins name duplication for using dRep
rule list_refined_genomes_after_unduplicating_filenames:
    input: "results/08_bins_postprocessing/genomes_list/genomes/{assembler}"
    output: "results/08_bins_postprocessing/genomes_list/{assembler}/list_unduplicated_filenames.txt"
    log:
        stdout = "logs/08_bins_postprocessing/genomes_list/{assembler}/list_unduplicated_filenames.stdout",
        stderr = "logs/08_bins_postprocessing/genomes_list/{assembler}/list_unduplicated_filenames.stderr"
    shell:
        """
        bash workflow/scripts/list_refined_genomes.sh {output} {input} \
            > {log.stdout} 2> {log.stderr}
        """

rule genomes_dereplication:
    # input is formed of every refined produced no matter the sample
    # it is, however, assembler specific
    input: "results/08_bins_postprocessing/genomes_list/{assembler}/list_unduplicated_filenames.txt"
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