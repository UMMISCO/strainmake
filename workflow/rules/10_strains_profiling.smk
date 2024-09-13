from utils import * 

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
ASSEMBLER = config['assembly']['assembler']

HYBRID_ASSEMBLER = config['assembly']['hybrid_assembler'] 

# taking into account the case where we don't have hybrid assembly
if HYBRID_ASSEMBLER == None:
       HYBRID_ASSEMBLER = []

# taking into account the case where we don't have SR assembly
if ASSEMBLER == None:
       ASSEMBLER = []

# rule to concatenate every bins that were dereplicated and filtered into a 
# unique FASTA file
rule creating_ref_genomes_fasta:
    input:
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/bins"
    output:
        "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa"
    log:
        stderr = "logs/10_strain_profiling/refs/{assembler}.concatenate.stderr"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        cat {input}/*.fa > {output} 2> {log.stderr}
        """

rule indexing_ref_genomes:
    input:
        "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa"
    output:
        "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa.fai"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/10_strain_profiling/refs_indexing/{assembler}.stdout",
        stderr = "logs/10_strain_profiling/refs_indexing/{assembler}.stderr"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        samtools faidx {input} > {log.stdout} 2> {log.stderr}
        """

# inStrain will take predicted genes as an input if they have been
# concatenated into a single FASTA file
rule concatenating_predicted_genes:
    input:
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/genes"
    output:
        "results/10_strain_profiling/refs/{assembler}/ref_genes.fa"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        cat {input}/*.fna > {output}
        """

# mapping sample reads on the reference genomes (the bins)
rule reads_mapping_on_reference:
    input:
        # the bins we concatenated into a single FASTA file
        refs = "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa",
        # metagenome reads
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        "results/10_strain_profiling/minimap2/{assembler}/{sample}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        stderr = "logs/10_strain_profiling/minimap2/{assembler}/{sample}.stderr"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER),
        sample="|".join(SAMPLES)
    threads: config['strains_profiling']['minimap2']['threads']
    shell:
        """
        minimap2 -ax sr -t {threads} \
            {input.refs} {input.r1} {input.r2} > {output} 2> {log.stderr}
        """

rule sam_to_bam_strains_profiling:
    input:
        sam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.sam"
    output:
        bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/10_strain_profiling/samtools/{assembler}/{sample}.sam_to_bam.stdout",
        stderr = "logs/10_strain_profiling/samtools/{assembler}/{sample}.sam_to_bam.stderr"
    wildcard_constraints:
        sample = "|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        samtools view -o {output.bam} {input.sam} \
            > {log.stdout} 2> {log.stderr} \
        && rm {input.sam}
        """

rule bam_sorting_strains_profiling:
    input:
        # reads mapped on the assembly
        bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.bam"
    output:
        bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/10_strain_profiling/samtools/{assembler}/{sample}.sorting.stdout",
        stderr = "logs/10_strain_profiling/samtools/{assembler}/{sample}.sorting.stderr"
    wildcard_constraints:
        sample="|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        samtools sort -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr} \
        && \
        rm {input.bam}
        """

rule bam_indexing:
    input:
        "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam"
    output:
        "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """ 
        samtools index {input}
        """

# producing a scaffolds to bin file for inStrain (file with the contig <-> bin link)
rule produce_scaffolds_to_bin_file:
    input:
        bins = "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/bins",
        refs = "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa"
    output:
        "results/10_strain_profiling/refs/{assembler}/ref_genomes.stb"
    conda:
        "../envs/drep.yaml"
    log:
        stdout = "logs/10_strain_profiling/inStrain/{assembler}/stb.stdout",
        stderr = "logs/10_strain_profiling/inStrain/{assembler}/stb.stderr"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    shell:
        """
        parse_stb.py --reverse -f {input.bins}/* -o {output} \
            > {log.stdout} 2> {log.stderr}
        """

# calling variants using the dereplicated bins as reference and the sample reads mapped on it
rule variant_calling:
    input:
        bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam",
        refs = "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa"
    output:
        "results/10_strain_profiling/freebayes/{assembler}/{sample}.vcf"
    conda:
        "../envs/freebayes.yaml"
    log:
        stderr = "logs/10_strain_profiling/freeboys/{assembler}/{sample}.variant_calling.stderr"
    params:
        min_alternate_count = config['strains_profiling']['freebayes']['min_alternate_count'],
        min_alternate_fraction = config['strains_profiling']['freebayes']['min_alternate_fraction']
    shell:
        """
        freebayes -f {input.refs} -F {params.min_alternate_fraction} -C {params.min_alternate_count} \
            --pooled-continuous \
            {input.bam} > {output} 2> {log.stderr}
        """

# strains profiling
rule instrain_profiling:
    input:
        bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam",
        refs = "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa",
        # scaffolds to bin file
        stb = "results/10_strain_profiling/refs/{assembler}/ref_genomes.stb",
        # predicted genes in references
        predicted_genes = "results/10_strain_profiling/refs/{assembler}/ref_genes.fa"
    output:
        directory("results/10_strain_profiling/inStrain/{assembler}/{sample}")
    conda:
        "../envs/instrain.yaml"
    log:
        stdout = "logs/10_strain_profiling/inStrain/{assembler}/{sample}.profile.stdout",
        stderr = "logs/10_strain_profiling/inStrain/{assembler}/{sample}.profile.stderr"
    wildcard_constraints:
        sample = "|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    threads: config['strains_profiling']['instrain']['threads']
    shell:
        """
        inStrain profile --output {output} -p {threads} \
            -s {input.stb} \
            --database_mode \
            -g {input.predicted_genes} \
            {input.bam} {input.refs} > {log.stdout} 2> {log.stderr}
        """

rule instrain_comparing: 
    input:
        instrain_results = expand("results/10_strain_profiling/inStrain/{{assembler}}/{sample}",
                                  sample=SAMPLES),
        stb = "results/10_strain_profiling/refs/{assembler}/ref_genomes.stb"
    output:
        directory("results/10_strain_profiling/inStrain/{assembler}/compare")
    conda:
        "../envs/instrain.yaml"
    log:
        stdout = "logs/10_strain_profiling/inStrain/{assembler}/compare.stdout",
        stderr = "logs/10_strain_profiling/inStrain/{assembler}/compare.stderr"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    threads: config['strains_profiling']['instrain']['threads']
    shell:
        """
        inStrain compare -i {input.instrain_results} --output {output} \
            --database_mode \
            -s {input.stb} \
            > {log.stdout} 2> {log.stderr}
        """

rule floria_profiling:
    input:
        bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam",
        indexed_bam = "results/10_strain_profiling/minimap2/{assembler}/{sample}.sorted.bam.bai",
        refs = "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa",
        # Floria needs the FASTA file to be indexed
        refs_indexed = "results/10_strain_profiling/refs/{assembler}/ref_genomes.fa.fai",
        vcf = "results/10_strain_profiling/freebayes/{assembler}/{sample}.vcf"
    output:
        "results/10_strain_profiling/floria/{assembler}/{sample}/contig_ploidy_info.tsv"
    conda:
        "../envs/floria.yaml"
    log:
        stdout = "logs/10_strain_profiling/floria/{assembler}/{sample}.profile.stdout",
        stderr = "logs/10_strain_profiling/floria/{assembler}/{sample}.profile.stderr"
    params:
        out_dir = "results/10_strain_profiling/floria/{assembler}/{sample}"
    wildcard_constraints:
        sample = "|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER)
    threads: config['strains_profiling']['floria']['threads']
    shell:
        """
        floria -b {input.bam} -v  {input.vcf} -r {input.refs} \
            --output-dir {params.out_dir} -t {threads} \
            --overwrite \
            > {log.stdout} 2> {log.stderr}
        """