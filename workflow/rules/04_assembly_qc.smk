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
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
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
    wildcard_constraints:
        assembler_lr = "|".join(ASSEMBLER_LR)
    shell:
        """
        gunzip -c {input} | prodigal -i /dev/stdin -d {output} -p meta \
            > {log.stdout} 2> {log.stderr}
        """

rule concatenating_assembly_genes:
    input:
        expand("results/04_assembly_qc/gene_calling/{{assembler}}/{sample}/genes.fna", sample=SAMPLES)
    output:
        "results/04_assembly_qc/gene_calling/{assembler}/genes.fna.gz"
    conda:
        "../envs/pigz.yaml"
    benchmark:
        "benchmarks/04_assembly_qc/gene_calling/{assembler}.benchmark.txt"
    params:
        uncompressed_output = "results/04_assembly_qc/gene_calling/{assembler}/genes.fna"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        cat {input} > {params.uncompressed_output} && pigz {params.uncompressed_output}
        """

rule concatenating_assembly_genes_long_read:
    input:
        expand("results/04_assembly_qc/gene_calling/{{assembler}}/{sample}/genes.fna", sample=SAMPLES_LR)
    output:
        "results/04_assembly_qc/gene_calling/{assembler}/genes.fna.gz"
    conda:
        "../envs/pigz.yaml"
    benchmark:
        "benchmarks/04_assembly_qc/gene_calling/{assembler}.benchmark.txt"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER_LR)
    shell:
        """
        cat {input} > {output} && pigz {output}
        """

# clustering genes to obtain a non redundant catalog
# we select best gene for cluster based on their length
rule gene_clustering:
    input:
        "results/04_assembly_qc/gene_calling/{assembler}/genes.fna.gz"
    output:
        "results/04_assembly_qc/gene_clustering/{assembler}/non_redundant_gene_catalog.fna.gz"
    conda:
        "../envs/cd-hit.yaml"
    benchmark:
        "benchmarks/04_assembly_qc/gene_clustering/{assembler}_gene_clustering.benchmark.txt"
    log:
        stdout = "logs/04_assembly_qc/gene_clustering/{assembler}_gene_clustering.stdout",
        stderr = "logs/04_assembly_qc/gene_clustering/{assembler}_gene_clustering.stderr"
    params:
        filtering_genes_cluster_script = "workflow/scripts/process_cd_hit_output.py",
        sequence_identity_threshold = config['cdhit']['sequence_identity_threshold'],
        alignment_coverage_shorter_sequence = config['cdhit']['alignment_coverage_shorter_sequence'],
        cdhit_output = "results/04_assembly_qc/gene_clustering/{assembler}/genes_clust",
        cdhit_output_clusters = "results/04_assembly_qc/gene_clustering/{assembler}/genes_clust.clstr",
        minimal_gene_length = config['representative_genes']['minimal_gene_length'],
        clusters_info = "results/04_assembly_qc/gene_clustering/{assembler}/genes_cluster.csv",
        uncompressed_output = "results/04_assembly_qc/gene_clustering/{assembler}/non_redundant_gene_catalog.fna"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR)
    threads: config['cdhit']['threads']
    shell:
        """
        cd-hit-est -i {input} -o {params.cdhit_output} -c {params.sequence_identity_threshold} -aS {params.alignment_coverage_shorter_sequence} \
            -G 0 -d 0 -M 0 -T {threads} \
            > {log.stdout} 2> {log.stderr} \
        && \
        python3 {params.filtering_genes_cluster_script} --minimal_len {params.minimal_gene_length} {params.cdhit_output_clusters} {input} {params.clusters_info} {params.uncompressed_output} \
            >> {log.stdout} 2>> {log.stderr} \
        && \
        pigz {params.uncompressed_output}
        """