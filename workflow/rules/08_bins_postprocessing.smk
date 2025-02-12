SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)

# adding to "SAMPLES" samples "SAMPLES_LR" that were not found in "SAMPLES"
for sample in SAMPLES_LR:
    if sample not in SAMPLES:
        SAMPLES.append(sample)

ANI_THRESHOLD = [str(ani) for ani in config['bins_postprocessing']['drep']['ani']]
DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE = str(config['bins_postprocessing']['genes_prediction']['prodigal']['ani'])

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

SHORT_READ_BINNER = config['binning']['binner']
LONG_READ_BINNER = config['binning']['long_read_binner']

# we would have refined bins with Binette only if we used several binning methods
refined = True if len(LONG_READ_BINNER) > 1 or len(SHORT_READ_BINNER) > 1 else False

wildcard_constraints:
    assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
    assembler_lr = "|".join(ASSEMBLER_LR) if ASSEMBLER_LR != [] else "none",
    assembler_sr_hybrid = "|".join(ASSEMBLER + HYBRID_ASSEMBLER) if ASSEMBLER + HYBRID_ASSEMBLER != [] else "none"

rule gtdb_tk_taxonomic_annotation:
    input:
        # folder with dereplicated and filtered bins (= MAG)
        refined_bins = "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins",
        ref_data = "data/gtdb_tk/release220"
    output: directory("results/08_bins_postprocessing/gtdb_tk/{ani}/{assembler}")
    conda:
        "../envs/gtdb_tk.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/gtdb_tk/{ani}/{assembler}/classify.stdout",
        stderr = "logs/08_bins_postprocessing/gtdb_tk/{ani}/{assembler}/classify.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/gtdb_tk/{ani}/{assembler}/classify.benchmark.txt"
    params:
        other_args = config['bins_postprocessing']['gtdbtk']['other_args']
    threads: config['bins_postprocessing']['gtdbtk']['threads']
    wildcard_constraints:
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.refined_bins} --cpus {threads} --out_dir {output} \
            --extension ".fa" \
            --skip_ani_screen --pplacer_cpus 1 \
            {params.other_args} \
            > {log.stdout} 2> {log.stderr}
        """

# this rule produces text files with list of refined genomes
if refined:
    rule list_refined_genomes_binette:
        input: 
            expand("results/07_bins_refinement/binette/{{assembler}}/{sample}", sample=SAMPLES) 
        output: "results/08_bins_postprocessing/genomes_list/{assembler}/list.txt"
        log:
            stdout = "logs/08_bins_postprocessing/genomes_list/{assembler}/list.stdout",
            stderr = "logs/08_bins_postprocessing/genomes_list/{assembler}/list.stderr"
        benchmark:
            "benchmarks/08_bins_postprocessing/genomes_list/{assembler}/list.benchmark.txt"
        params:
            # constructing the precise folder path with bins from the input
            bins_folder = lambda wildcards, input: [f"{dir}/final_bins" for dir in input]
        shell:
            """
            bash workflow/scripts/list_refined_genomes.sh {output} {params.bins_folder} \
                > {log.stdout} 2> {log.stderr}
            """
else:
    rule list_refined_genomes_no_binette_SR:
        input:
            expand("results/05_binning/{binner}/bins/{{assembler_sr_hybrid}}/{sample}", binner=SHORT_READ_BINNER, sample=SAMPLES)
        output: "results/08_bins_postprocessing/genomes_list/{assembler_sr_hybrid}/list.txt"
        log:
            stdout = "logs/08_bins_postprocessing/genomes_list/{assembler_sr_hybrid}/list.stdout",
            stderr = "logs/08_bins_postprocessing/genomes_list/{assembler_sr_hybrid}/list.stderr"
        benchmark:
            "benchmarks/08_bins_postprocessing/genomes_list/{assembler_sr_hybrid}/list.benchmark.txt"
        params:
            # constructing the precise folder path with bins from the input
            bins_folder = lambda wildcards, input: [f"{dir}/bins" for dir in input]
        shell:
            """
            bash workflow/scripts/list_refined_genomes.sh {output} {params.bins_folder} \
                > {log.stdout} 2> {log.stderr}
            """
    rule list_refined_genomes_no_binette_LR:
        input:
            expand("results/05_binning/{binner}/bins/{{assembler_lr}}/{sample}", binner=LONG_READ_BINNER, sample=SAMPLES_LR)
        output: "results/08_bins_postprocessing/genomes_list/{assembler_lr}/list.txt"
        log:
            stdout = "logs/08_bins_postprocessing/genomes_list/{assembler_lr}/list.stdout",
            stderr = "logs/08_bins_postprocessing/genomes_list/{assembler_lr}/list.stderr"
        benchmark:
            "benchmarks/08_bins_postprocessing/genomes_list/{assembler_lr}/list.benchmark.txt"
        params:
            # constructing the precise folder path with bins from the input
            bins_folder = lambda wildcards, input: [f"{dir}/bins" for dir in input]
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
    benchmark:
        "benchmarks/08_bins_postprocessing/genomes_list/{assembler}/copy_and_rename.benchmark.txt"
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
    benchmark:
        "benchmarks/08_bins_postprocessing/genomes_list/{assembler}/list_unduplicated_filenames.benchmark.txt"
    shell:
        """
        bash workflow/scripts/list_refined_genomes.sh {output} {input} \
            > {log.stdout} 2> {log.stderr}
        """

rule genomes_dereplication:
    # input is formed of every refined produced no matter the sample
    # it is, however, assembler specific
    input: "results/08_bins_postprocessing/genomes_list/{assembler}/list_unduplicated_filenames.txt"
    output: directory("results/08_bins_postprocessing/dRep/{ani}/{assembler}")
    conda:
        "../envs/drep.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/drep/{ani}/{assembler}/dereplication.stdout",
        stderr = "logs/08_bins_postprocessing/drep/{ani}/{assembler}/dereplication.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/drep/{ani}/{assembler}/dereplication.benchmark.txt"
    params:
        comparison_algorithm = config['bins_postprocessing']['drep']['comparison_algorithm'],
        other_args = config['bins_postprocessing']['drep']['other_args'],
        ani_dec = lambda wildcards: f"{float(wildcards.ani) / 100}"
    threads: config['bins_postprocessing']['drep']['threads']
    wildcard_constraints:
        ani = "|".join(ANI_THRESHOLD)
    shell:
        """
        dRep dereplicate --genomes {input} --processors {threads} \
            --S_algorithm {params.comparison_algorithm} \
            --S_ani {params.ani_dec} \
            {params.other_args} \
            {output} \
        > {log.stdout} 2> {log.stderr}
        """

# this rule will perform dereplicated bins quality estimation 
# and filter the bins based on user settings
rule dereplicated_genomes_quality_and_filtering:
    input:
        # folder with dereplicated bins
        bins = "results/08_bins_postprocessing/dRep/{ani}/{assembler}", 
        diamond_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        out_dir = directory("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/checkm2"),
        selected_bins = directory("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins")
    conda:
        "../envs/checkm2.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/checkm2/{ani}/{assembler}.assessment.stdout",
        stderr = "logs/08_bins_postprocessing/checkm2/{ani}/{assembler}.assessment.stderr",
        stdout_filtration = "logs/08_bins_postprocessing/genomes_filtration/{ani}/{assembler}.filtration.stdout",
        stderr_filtration = "logs/08_bins_postprocessing/genomes_filtration/{ani}/{assembler}.filtration.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/checkm2/{ani}/{assembler}.assessment_and_filtration.benchmark.txt"
    params:
        minimal_completeness = config['bins_postprocessing']['genomes_quality_filtration']['filtration']['min_completeness'],
        maximal_contamination = config['bins_postprocessing']['genomes_quality_filtration']['filtration']['max_contamination']
    threads: config['bins_postprocessing']['genomes_quality_filtration']['checkm2']['threads']
    wildcard_constraints:
        ani = "|".join(ANI_THRESHOLD)
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
            > {log.stdout_filtration} 2> {log.stderr_filtration} \
        && \
        python3 workflow/scripts/deduplicate_contigs_name.py {output.selected_bins}
        """

# predicting genes in dereplicated genomes
rule genes_calling:
    input:
        # it is better to run Prodigal on each genome individually in normal mode, than running it on 
        # a multiple FASTA file 
        # (https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#metagenomes)
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins"
    output:
        directory("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/genes")
    conda:
        "../envs/prodigal.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/prodigal/{ani}/{assembler}.stdout",
        stderr = "logs/08_bins_postprocessing/prodigal/{ani}/{assembler}.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/prodigal/{ani}/{assembler}.benchmark.txt"
    threads: config['bins_postprocessing']['genes_prediction']['prodigal']['threads']
    wildcard_constraints:
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        python3 workflow/scripts/genes_prediction.py --cpu {threads} {input} {output} \
            > {log.stdout} 2> {log.stderr}
        """

# estimating distribution of dereplicated bins in the metagenomes using CheckM 1
# first we calculate coverage
# then, we profile the metagenome using the coverage
rule coverage_in_mapping:
    input:
        dereplicated_and_filtered_bins = "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins",
        samples_mapped_on_dereplicated_and_filtered_bins = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam",
        # needing the indexed BAM also
        mapping_index = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam.bai"
    output:
        "results/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}/coverage.tsv"
    conda:
        "../envs/checkm1.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.coverage.stdout",
        stderr = "logs/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.coverage.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.coverage.benchmark.txt"
    threads: config['bins_postprocessing']['profiling']['checkm1']['threads']
    wildcard_constraints:
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        checkm coverage -x fa {input.dereplicated_and_filtered_bins} {output} \
            {input.samples_mapped_on_dereplicated_and_filtered_bins} \
            > {log.stdout} 2> {log.stderr}
        """

rule bins_distribution_estimation:
    input: "results/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}/coverage.tsv"
    output: "results/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}/profile.tsv"
    conda:
        "../envs/checkm1.yaml"
    log: 
        stdout = "logs/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.profile.stdout",
        stderr = "logs/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.profile.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.profile.benchmark.txt"
    wildcard_constraints:
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        checkm profile --tab_table -f  {output} {input} > {log.stdout} 2> {log.stderr}
        """

rule process_estimated_bins_distribution:
    input: "results/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}/profile.tsv"
    output: "results/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}/profile.processed.tsv"
    conda:
        "../envs/python.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.profile_processing_table.stdout",
        stderr = "logs/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.profile_processing_table.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/checkm1/{ani}/{assembler}/{sample}.profile_processing_table.benchmark.txt"
    wildcard_constraints:
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        python3 workflow/scripts/process_checkm1_profile_table.py --input_table {input} \
            {output} > {log.stdout} 2> {log.stderr}
        """

# reconstructing metabolic models using CarveMe for MAGs
rule carveme_models_building:
    input:
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins"
    output:
        directory("results/08_bins_postprocessing/carveme/{ani}/{assembler}")
    conda:
        "../envs/carveme.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/carveme/{ani}/{assembler}.stdout",
        stderr = "logs/08_bins_postprocessing/carveme/{ani}/{assembler}.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/carveme/{ani}/{assembler}.benchmark.txt"
    wildcard_constraints:
        ani = "|".join(ANI_THRESHOLD)
    params:
        launch_script = "workflow/scripts/carveme_models_building.py"
    threads:
        config['bins_postprocessing']['carveme']['threads']
    shell:
        """
       python3 {params.launch_script} carve --cpu {threads} -i {input} -o {output} -v > {log.stdout} 2> {log.stderr}
       """

# merging organisms' metabolic models into a community model
# we then compress all individual models
rule carveme_merge_models:
    input:
        "results/08_bins_postprocessing/carveme/{ani}/{assembler}"
    output:
        "results/08_bins_postprocessing/carveme/{ani}/{assembler}/community_model/community.xml"
    conda:
        "../envs/carveme.yaml"
    log:
        stdout = "logs/08_bins_postprocessing/carveme/{ani}/{assembler}/community_model.stdout",
        stderr = "logs/08_bins_postprocessing/carveme/{ani}/{assembler}/community_model.stderr"
    benchmark:
        "benchmarks/08_bins_postprocessing/carveme/{ani}/{assembler}/community_model.benchmark.txt"
    wildcard_constraints:
        ani = "|".join(ANI_THRESHOLD)
    shell:
        """
        python3 {params.launch_script} merge -i {input} -o {output} -v > {log.stdout} 2> {log.stderr} \
        && \
        pigz {input}/*.xml
        """