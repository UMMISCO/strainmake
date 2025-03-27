import os 

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

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)

# adding to "SAMPLES" samples "SAMPLES_LR" that were not found in "SAMPLES"
for sample in SAMPLES_LR:
    if sample not in SAMPLES:
        SAMPLES.append(sample)

# mananing genomes of reference (if any)
reference_genomes_dir = config['quast']['reference_genomes_dir']
if reference_genomes_dir is not None:
    reference_genomes_extension = config['quast']['reference_genomes_extension']
    REFERENCE_GENOMES = ",".join(
        os.path.join(reference_genomes_dir, f) for f in os.listdir(reference_genomes_dir) if f.endswith(reference_genomes_extension)
    )
    # adding the `metaquast` flag (`-r genome1,genome2,...`) to the `REFERENCE_GENOMES` variable
    REFERENCE_GENOMES = f"-r {REFERENCE_GENOMES}"
else:
    REFERENCE_GENOMES = "none"

wildcard_constraints:
       assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER) if ASSEMBLER + HYBRID_ASSEMBLER != [] else "none",
       assembler_lr = "|".join(ASSEMBLER_LR) if ASSEMBLER_LR != [] else "none",
       assembler_all = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR) if ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR != [] else "none", 

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
        out_dir = "results/04_assembly_qc/quast/{assembler}/{sample}",
        ref_genomes = REFERENCE_GENOMES if REFERENCE_GENOMES != "none" else "",
    threads: config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            --circos \
            {params.ref_genomes} \
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
        out_dir = "results/04_assembly_qc/quast/{assembler_lr}/{sample_lr}",
        ref_genomes = REFERENCE_GENOMES if REFERENCE_GENOMES != "none" else "",
        method = "--nanopore" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio"
    threads: config['quast']['threads']
    shell:
        """
        metaquast.py -t {threads} -o {params.out_dir} \
            --max-ref-number 0 \
            --circos \
            {params.ref_genomes} \
            {input.assembly} \
            > {log.stdout} 2> {log.stderr} 
        """

# building a non redundant gene catalog for each assembly approach
rule gene_calling_assembly:
    input:
        # assemblies produced in step 03
        "results/03_assembly/{assembler}/{sample}/assembly.fa.gz"
    output:
        "results/04_assembly_qc/gene_calling/{assembler}/{sample}/genes.fna"
    conda:
        "../envs/prodigal.yaml"
    log:
        stdout = "logs/04_assembly_qc/gene_calling/{assembler}/{sample}.stdout",
        stderr = "logs/04_assembly_qc/gene_calling/{assembler}/{sample}.stderr"
    benchmark:
        "benchmarks/04_assembly_qc/gene_calling/{assembler}/{sample}.benchmark.txt"
    shell:
        """
        gunzip -c {input} | prodigal -i /dev/stdin -d {output} -p meta \
            > {log.stdout} 2> {log.stderr}
        """

rule gene_calling_assembly_long_read:
    input:
        # assemblies produced in step 03
        "results/03_assembly/{assembler_lr}/{sample}/assembly.fa.gz"
    output:
        "results/04_assembly_qc/gene_calling/{assembler_lr}/{sample}/genes.fna"
    conda:
        "../envs/prodigal.yaml"
    log:
        stdout = "logs/04_assembly_qc/gene_calling/{assembler_lr}/{sample}.stdout",
        stderr = "logs/04_assembly_qc/gene_calling/{assembler_lr}/{sample}.stderr"
    benchmark:
        "benchmarks/04_assembly_qc/gene_calling/{assembler_lr}/{sample}.benchmark.txt"
    shell:
        """
        gunzip -c {input} | prodigal -i /dev/stdin -d {output} -p meta \
            > {log.stdout} 2> {log.stderr}
        """

# concatenating all genes from all samples for each assembly approach and making sequence names unique
# we also add the assembly method to the sequence ID and replace the spaces induced by Prodigal by "__" to not make sequences be redundant/duplicated for MMseqs2
rule concatenating_assembly_genes:
    input:
        expand("results/04_assembly_qc/gene_calling/{{assembler}}/{sample}/genes.fna", sample=SAMPLES)
    output:
        "results/04_assembly_qc/gene_calling/{assembler}/genes.fna.gz"
    conda:
        "../envs/seqkit.yaml"
    benchmark:
        "benchmarks/04_assembly_qc/gene_calling/{assembler}.benchmark.txt"
    params:
        uncompressed_output = "results/04_assembly_qc/gene_calling/{assembler}/genes.fna"
    shell:
        """
        cat {input} | seqkit rename | seqkit replace -p " " -r "__" | seqkit replace -p "(.*)" -r "{wildcards.assembler}___\\$1" -o {params.uncompressed_output} \
            && pigz {params.uncompressed_output}
        """

# concatenating all genes from all samples for each assembly approach and making sequence names unique
# we also add the assembly method to the sequence ID and replace the spaces induced by Prodigal by "__" to not make sequences be redundant/duplicated for MMseqs2
rule concatenating_assembly_genes_long_read:
    input:
        expand("results/04_assembly_qc/gene_calling/{{assembler_lr}}/{sample}/genes.fna", sample=SAMPLES_LR)
    output:
        "results/04_assembly_qc/gene_calling/{assembler_lr}/genes.fna.gz"
    conda:
        "../envs/seqkit.yaml"
    benchmark:
        "benchmarks/04_assembly_qc/gene_calling/{assembler_lr}.benchmark.txt"
    params:
        uncompressed_output = "results/04_assembly_qc/gene_calling/{assembler_lr}/genes.fna"
    shell:
        """
        cat {input} | seqkit rename | seqkit replace -p " " -r "__" | seqkit replace -p "(.*)" -r "{wildcards.assembler_lr}___\\$1" -o {params.uncompressed_output} \
            && pigz {params.uncompressed_output}
        """

# clustering genes to obtain a non redundant catalog
# we select best gene for cluster based on their length
rule gene_clustering:
    input:
        "results/04_assembly_qc/gene_calling/{assembler_all}/genes.fna.gz"
    output:
        "results/04_assembly_qc/gene_clustering/{assembler_all}/non_redundant_gene_catalog.fna.gz"
    conda:
        "../envs/mmseqs2.yaml"
    benchmark:
        "benchmarks/04_assembly_qc/gene_clustering/{assembler_all}_gene_clustering.benchmark.txt"
    log:
        stdout = "logs/04_assembly_qc/gene_clustering/{assembler_all}_gene_clustering.stdout",
        stderr = "logs/04_assembly_qc/gene_clustering/{assembler_all}_gene_clustering.stderr"
    params:
        sequence_identity_threshold = config['mmseqs2']['sequence_identity_threshold'],
        alignment_coverage_shorter_sequence = config['mmseqs2']['alignment_coverage_shorter_sequence'],
        minimal_gene_length = config['representative_genes']['minimal_gene_length'],
        tmp = "{assembler_all}_tmp",
        mmseqs2_output = "results/04_assembly_qc/gene_clustering/{assembler_all}/clustering",
        uncompressed_output = "results/04_assembly_qc/gene_clustering/{assembler_all}/non_redundant_gene_catalog.fna",
        uncompressed_output_before_renaming = "results/04_assembly_qc/gene_clustering/{assembler_all}/clustering_rep_seq.fasta",
    threads: config['mmseqs2']['threads']
    shell:
        """
        mmseqs easy-cluster {input} {params.mmseqs2_output} {params.tmp} \
            --min-seq-id {params.sequence_identity_threshold} -c {params.alignment_coverage_shorter_sequence} \
            --alignment-mode 3 \
            --cov-mode 1 \
            --dbtype 2 \
            --threads {threads} \
        > {log.stdout} 2> {log.stderr} \
        && \
        seqkit seq -m {params.minimal_gene_length} {params.uncompressed_output_before_renaming} > {params.uncompressed_output} \
        && \
        rm {params.uncompressed_output_before_renaming} \
        && \
        pigz {params.uncompressed_output}
        """