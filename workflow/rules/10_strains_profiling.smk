from utils import * 

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
ASSEMBLER = config['assembly']['assembler']
ASSEMBLER_LR = config['assembly']['long_read_assembler']

HYBRID_ASSEMBLER = config['assembly']['hybrid_assembler'] 

# taking into account the case where we don't have hybrid assembly
if HYBRID_ASSEMBLER == None:
       HYBRID_ASSEMBLER = []

# taking into account the case where we don't have SR assembly
if ASSEMBLER == None:
       ASSEMBLER = []

# taking into account the case where we don't have LR
if ASSEMBLER_LR == None:
       ASSEMBLER_LR = []

DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE = str(config['bins_postprocessing']['genes_prediction']['prodigal']['ani'])

# allows a flexibility for the user to use sequences in FASTA or FASTQ format
# (minimap2 can map sequences from FASTA or FASTQ)
seq_format = config["lr_seq_format"]
sequences_file_end = f"_1.{seq_format}.gz"

wildcard_constraints:
    assembler_lr = "|".join(ASSEMBLER_LR) if ASSEMBLER_LR != [] else "none",
    assembler_all = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR) if ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR != [] else "none",
    assembler_sr = "|".join(ASSEMBLER) if ASSEMBLER != [] else "none",
    assembler_hybrid = "|".join(HYBRID_ASSEMBLER) if HYBRID_ASSEMBLER != [] else "none",
    assembler_sr_hybrid = "|".join(ASSEMBLER + HYBRID_ASSEMBLER) if ASSEMBLER + HYBRID_ASSEMBLER != [] else "none",

# rule to concatenate every bins that were dereplicated and filtered into a 
# unique FASTA file
rule creating_ref_genomes_fasta:
    input:
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins"
    output:
        "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa"
    log:
        stderr = "logs/10_strain_profiling/refs/{ani}/{assembler}.concatenate.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/refs/{ani}/{assembler}.concatenate.benchmark.txt"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        cat {input}/*.fa > {output} 2> {log.stderr}
        """

rule indexing_ref_genomes:
    input:
        "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa"
    output:
        "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa.fai"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/10_strain_profiling/refs_indexing/{ani}/{assembler}.stdout",
        stderr = "logs/10_strain_profiling/refs_indexing/{ani}/{assembler}.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/refs_indexing/{ani}/{assembler}.benchmark.txt"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        samtools faidx {input} > {log.stdout} 2> {log.stderr}
        """

# inStrain will take predicted genes as an input if they have been
# concatenated into a single FASTA file
rule concatenating_predicted_genes:
    input:
        "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/genes"
    output:
        "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genes.fa"
    benchmark:
        "benchmarks/10_strain_profiling/refs/{ani}/{assembler}/concatenate_ref_genes.benchmark.txt"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        cat {input}/*.fna > {output}
        """

# mapping sample reads on the reference genomes (the bins)
rule reads_mapping_on_reference:
    input:
        # the bins we concatenated into a single FASTA file
        refs = "results/10_strain_profiling/refs/{ani}/{assembler_sr}/ref_genomes.fa",
        # metagenome reads
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        "results/10_strain_profiling/minimap2/{ani}/{assembler_sr}/{sample}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        stderr = "logs/10_strain_profiling/minimap2/{ani}/{assembler_sr}/{sample}.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/minimap2/{ani}/{assembler_sr}/{sample}.benchmark.txt"
    wildcard_constraints:
        sample="|".join(SAMPLES),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    threads: config['strains_profiling']['minimap2']['threads']
    shell:
        """
        minimap2 -ax sr -t {threads} \
            {input.refs} {input.r1} {input.r2} > {output} 2> {log.stderr}
        """

# there are two possibilities: reads used in hybrid assemblies were downsized, or reads used 
# in hybrid assemblies were not downsized
subsample_hybrid_reads = config["downsizing_for_hybrid"]["lr"] is not None and config["downsizing_for_hybrid"]["sr"] is not None

rule reads_mapping_on_reference_hybrid:
    input:
        # the bins we concatenated into a single FASTA file
        refs = "results/10_strain_profiling/refs/{ani}/{assembler_hybrid}/ref_genomes.fa",
        # Select input files based on whether reads are downsized or not
        r1 = lambda wildcards: f"results/02_preprocess/{'downsized/' if subsample_hybrid_reads else ''}bowtie2/{wildcards.sample}_1.clean{'.downsized' if subsample_hybrid_reads else ''}.fastq.gz",
        r2 = lambda wildcards: f"results/02_preprocess/{'downsized/' if subsample_hybrid_reads else ''}bowtie2/{wildcards.sample}_2.clean{'.downsized' if subsample_hybrid_reads else ''}.fastq.gz",
        long_read = lambda wildcards: f"results/02_preprocess/{'downsized/' if subsample_hybrid_reads else ''}fastp_long_read/{wildcards.sample}{'_downsized' if subsample_hybrid_reads else ''}{sequences_file_end}"
    output:
        "results/10_strain_profiling/minimap2/{ani}/{assembler_hybrid}/{sample}.sorted.bam"
    conda:
        "../envs/minimap2.yaml"
    log:
        sr_stderr = "logs/10_strain_profiling/minimap2/{ani}/{assembler_hybrid}/{sample}.SR.stderr",
        lr_stderr = "logs/10_strain_profiling/minimap2/{ani}/{assembler_hybrid}/{sample}.LR.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/minimap2/{ani}/{assembler_hybrid}/{sample}.benchmark.txt"
    wildcard_constraints:
        sample="|".join(SAMPLES),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    params:
        mapping_sr = "results/10_strain_profiling/minimap2/{ani}/{assembler_hybrid}/{sample}.SR.bam",
        mapping_lr = "results/10_strain_profiling/minimap2/{ani}/{assembler_hybrid}/{sample}.LR.bam",
        method = (
            "map-ont" if config.get('lr_technology', '') == "nanopore"
            else "map-hifi" if config.get('lr_technology', '') == "pacbio-hifi"
            else "map-pb" if config.get('lr_technology', '') == "pacbio"
            else ""
        ),
    threads: config['strains_profiling']['minimap2']['threads']
    shell:
        """
        minimap2 -ax sr -t {threads} \
            {input.refs} {input.r1} {input.r2} | samtools view -Sb - | samtools sort -o {params.mapping_sr} 2> {log.sr_stderr} \
        && \
        minimap2 -ax {params.method} -t {threads} \
            {input.refs} {input.long_read} | samtools view -Sb - | samtools sort -o {params.mapping_lr} 2> {log.lr_stderr} \
        && \
        samtools merge --threads {threads} -o {output} {params.mapping_sr} {params.mapping_lr} \
        && \
        rm -v {params.mapping_sr} {params.mapping_lr}
        """

# mapping sample reads (long reads) on the reference genomes (the bins)
rule reads_LR_mapping_on_reference:
    input:
        # the bins we concatenated into a single FASTA file
        refs = "results/10_strain_profiling/refs/{ani}/{assembler_lr}/ref_genomes.fa",
        # metagenome reads
        long_read = "results/02_preprocess/fastp_long_read/{sample}" + sequences_file_end
    output:
        "results/10_strain_profiling/minimap2/{ani}/{assembler_lr}/{sample}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        stderr = "logs/10_strain_profiling/minimap2/{ani}/{assembler_lr}/{sample}.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/minimap2/{ani}/{assembler_lr}/{sample}.benchmark.txt"
    wildcard_constraints:
        sample="|".join(SAMPLES),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    params:
        method = (
            "map-ont" if config.get('lr_technology', '') == "nanopore"
            else "map-hifi" if config.get('lr_technology', '') == "pacbio-hifi"
            else "map-pb" if config.get('lr_technology', '') == "pacbio"
            else ""
        ),
    threads: config['strains_profiling']['minimap2']['threads']
    shell:
        """
        minimap2 -ax {params.method} -t {threads} \
            {input.refs} {input.long_read} \
            > {output} 2> {log.stderr}
        """

rule sam_to_bam_strains_profiling:
    input:
        sam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sam"
    output:
        bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/10_strain_profiling/samtools/{ani}/{assembler}/{sample}.sam_to_bam.stdout",
        stderr = "logs/10_strain_profiling/samtools/{ani}/{assembler}/{sample}.sam_to_bam.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/samtools/{ani}/{assembler}/{sample}.sam_to_bam.benchmark.txt"
    wildcard_constraints:
        sample = "|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    threads: config['binning']['samtools']['threads']
    shell:
        """
        samtools view --threads {threads} -o {output.bam} {input.sam} \
            > {log.stdout} 2> {log.stderr} \
        && rm {input.sam}
        """

rule bam_sorting_strains_profiling:
    input:
        # reads mapped on the assembly
        bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.bam"
    output:
        bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/10_strain_profiling/samtools/{ani}/{assembler}/{sample}.sorting.stdout",
        stderr = "logs/10_strain_profiling/samtools/{ani}/{assembler}/{sample}.sorting.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/samtools/{ani}/{assembler}/{sample}.sorting.benchmark.txt"
    wildcard_constraints:
        sample="|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        samtools sort -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr} \
        && \
        rm {input.bam}
        """

rule bam_indexing:
    input:
        "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam"
    output:
        "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.bam_indexing.benchmark.txt"
    wildcard_constraints:
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """ 
        samtools index {input}
        """

# producing a scaffolds to bin file for inStrain (file with the contig <-> bin link)
rule produce_scaffolds_to_bin_file:
    input:
        bins = "results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins",
        refs = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa"
    output:
        "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.stb"
    conda:
        "../envs/drep.yaml"
    log:
        stdout = "logs/10_strain_profiling/inStrain/{ani}/{assembler}/stb.stdout",
        stderr = "logs/10_strain_profiling/inStrain/{ani}/{assembler}/stb.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/inStrain/{ani}/{assembler}/stb.benchmark.txt"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        parse_stb.py --reverse -f {input.bins}/* -o {output} \
            > {log.stdout} 2> {log.stderr}
        """

# calling variants using the dereplicated bins as reference and the sample reads mapped on it
rule variant_calling:
    input:
        bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam",
        refs = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa"
    output:
        "results/10_strain_profiling/freebayes/{ani}/{assembler}/{sample}.vcf"
    conda:
        "../envs/freebayes.yaml"
    log:
        stderr = "logs/10_strain_profiling/freeboys/{ani}/{assembler}/{sample}.variant_calling.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/freeboys/{ani}/{assembler}/{sample}.variant_calling.benchmark.txt"
    params:
        min_alternate_count = config['strains_profiling'].get('freebayes', {}).get('min_alternate_count'),
        min_alternate_fraction = config['strains_profiling'].get('freebayes', {}).get('min_alternate_fraction')
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    shell:
        """
        freebayes -f {input.refs} -F {params.min_alternate_fraction} -C {params.min_alternate_count} \
            --pooled-continuous \
            {input.bam} > {output} 2> {log.stderr}
        """

# strains profiling
rule instrain_profiling:
    input:
        bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam",
        refs = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa",
        # scaffolds to bin file
        stb = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.stb",
        # predicted genes in references
        predicted_genes = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genes.fa"
    output:
        directory("results/10_strain_profiling/inStrain/{ani}/{assembler}/{sample}")
    conda:
        "../envs/instrain.yaml"
    log:
        stdout = "logs/10_strain_profiling/inStrain/{ani}/{assembler}/{sample}.profile.stdout",
        stderr = "logs/10_strain_profiling/inStrain/{ani}/{assembler}/{sample}.profile.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/inStrain/{ani}/{assembler}/{sample}.profile.benchmark.txt"
    wildcard_constraints:
        sample = "|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
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
        instrain_results = expand("results/10_strain_profiling/inStrain/{{ani}}/{{assembler}}/{sample}",
                                  sample=SAMPLES),
        stb = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.stb"
    output:
        directory("results/10_strain_profiling/inStrain/{ani}/{assembler}/compare")
    conda:
        "../envs/instrain.yaml"
    log:
        stdout = "logs/10_strain_profiling/inStrain/{ani}/{assembler}/compare.stdout",
        stderr = "logs/10_strain_profiling/inStrain/{ani}/{assembler}/compare.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/inStrain/{ani}/{assembler}/compare.benchmark.txt"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER + ASSEMBLER_LR),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
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
        bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam",
        indexed_bam = "results/10_strain_profiling/minimap2/{ani}/{assembler}/{sample}.sorted.bam.bai",
        refs = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa",
        # Floria needs the FASTA file to be indexed
        refs_indexed = "results/10_strain_profiling/refs/{ani}/{assembler}/ref_genomes.fa.fai",
        vcf = "results/10_strain_profiling/freebayes/{ani}/{assembler}/{sample}.vcf"
    output:
        "results/10_strain_profiling/floria/{ani}/{assembler}/{sample}/contig_ploidy_info.tsv"
    conda:
        "../envs/floria.yaml"
    log:
        stdout = "logs/10_strain_profiling/floria/{ani}/{assembler}/{sample}.profile.stdout",
        stderr = "logs/10_strain_profiling/floria/{ani}/{assembler}/{sample}.profile.stderr"
    benchmark:
        "benchmarks/10_strain_profiling/floria/{ani}/{assembler}/{sample}.profile.benchmark.txt"
    params:
        out_dir = "results/10_strain_profiling/floria/{ani}/{assembler}/{sample}"
    wildcard_constraints:
        sample = "|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER + HYBRID_ASSEMBLER),
        ani = DEREPLICATED_GENOMES_THRESHOLD_TO_PROFILE
    threads: config['strains_profiling'].get('floria', {}).get('threads', 0)
    shell:
        """
        floria -b {input.bam} -v  {input.vcf} -r {input.refs} \
            --output-dir {params.out_dir} -t {threads} \
            --overwrite \
            > {log.stdout} 2> {log.stderr}
        """