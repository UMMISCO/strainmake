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
    params:
        other_args = config['bins_postprocessing']['gtdbtk']['other_args']
    threads: config['bins_postprocessing']['gtdbtk']['threads']
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.refined_bins}/final_bins --cpus {threads} --out_dir {output} \
            --extension ".fa" \
            --skip_ani_screen --pplacer_cpus 1 \
            {params.other_args} \
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

# this rule will perform dereplicated bins quality estimation 
# and filter the bins based on user settings
rule dereplicated_genomes_quality_and_filtering:
    input:
        # folder with dereplicated bins
        bins = "results/08_bins_postprocessing/dRep/{assembler}", 
        diamond_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        out_dir = directory("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/checkm2"),
        selected_bins = directory("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/bins")
    conda:
        "../envs/checkm2.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/checkm2/{assembler}.assessment.stdout",
        stderr = "logs/08_bins_postprocessing/checkm2/{assembler}.assessment.stderr",
        stdout_filtration = "logs/08_bins_postprocessing/genomes_filtration/{assembler}.filtration.stdout",
        stderr_filtration = "logs/08_bins_postprocessing/genomes_filtration/{assembler}.filtration.stderr"
    params:
        minimal_completeness = config['bins_postprocessing']['genomes_quality_filtration']['filtration']['min_completeness'],
        maximal_contamination = config['bins_postprocessing']['genomes_quality_filtration']['filtration']['max_contamination']
    threads: config['bins_postprocessing']['genomes_quality_filtration']['checkm2']['threads']
    shell:
        """
        checkm2 predict --input {input.bins}/dereplicated_genomes --threads {threads} \
            -x .fa \
            --database_path {input.diamond_database} \
            --output-directory {output.out_dir} > {log.stdout} 2> {log.stderr} \
        && \
        python3 workflow/scripts/filter_dereplicated_bins_by_quality.py \
            --checkm-report {output.out_dir}/quality_report.tsv \
            --bins-directory {input.bins}/dereplicated_genomes \
            --min-comp {params.minimal_completeness} \
            --max-cont {params.maximal_contamination} \
            --outdir {output.selected_bins} \
            > {log.stdout_filtration} 2> {log.stderr_filtration}
        """

# predicting genes in dereplicated genomes
rule genes_calling:
    input:
        # it is better to run Prodigal on each genome individually in normal mode, than running it on 
        # a multiple FASTA file 
        # (https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#metagenomes)
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/bins"
    output:
        directory("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/genes")
    conda:
        "../envs/prodigal.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/prodigal/{assembler}.stdout",
        stderr = "logs/08_bins_postprocessing/prodigal/{assembler}.stderr"
    threads: config['bins_postprocessing']['genes_prediction']['prodigal']['threads']
    shell:
        """
        python3 workflow/scripts/genes_prediction.py --cpu {threads} {input} {output} \
            > {log.stdout} 2> {log.stderr}
        """