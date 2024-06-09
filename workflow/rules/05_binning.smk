from utils import * 

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
ASSEMBLER = config['assembly']['assembler']

# rules for mapping reads on assembly and producing a .BAM file
rule bowtie_assembly_index:
    input:
        # assemblies produced in step 03
        "results/03_assembly/{assembler}/{sample}/assembly.fa.gz"
    output:
        directory("results/05_binning/bowtie2/index/{assembler}/{sample}")
    conda:
        "../envs/bowtie2.yaml"
    log:
        stdout = "logs/05_binning/bowtie2/{assembler}/{sample}.indexing.stdout",
        stderr = "logs/05_binning/bowtie2/{assembler}/{sample}.indexing.stderr"
    params:
        seed = config['binning']['bowtie2']['seed'],
        index_basename = "{sample}",
        assembler = config['assembly']['assembler']
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER)
    threads: config['binning']['bowtie2']['threads']
    shell:
        """
        mkdir -p {output} \
        && \
        bowtie2-build --threads {threads} --seed {params.seed} \
            {input} "{output}/{params.index_basename}" \
            > {log.stdout} 2> {log.stderr}
        """

rule reads_mapping:
    input:
        # metagenome reads
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz",
        # assemblies index produced in rule "bowtie_assembly_index"
        bowtie_index = "results/05_binning/bowtie2/index/{assembler}/{sample}"
    output:
        sam = "results/05_binning/bowtie2/{assembler}/{sample}.sam"
    conda:
        "../envs/bowtie2.yaml"
    log:
        stdout = "logs/05_binning/bowtie2/{assembler}/{sample}.mapping.stdout",
        stderr = "logs/05_binning/bowtie2/{assembler}/{sample}.mapping.stderr"
    params:
        index_basename = "{sample}",
        assembler = config['assembly']['assembler']
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER)
    threads: config['binning']['bowtie2']['threads']
    shell:
        """
        bowtie2 -p {threads} \
            -x "{input.bowtie_index}/{params.index_basename}" \
            -1 {input.r1} -2 {input.r2} \
            -S {output.sam} \
            > {log.stdout} 2> {log.stderr}
        """

rule sam_to_bam:
    input:
        sam = "results/05_binning/bowtie2/{assembler}/{sample}.sam"
    output:
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/{assembler}/{sample}.sam_to_bam.stdout",
        stderr = "logs/05_binning/samtools/{assembler}/{sample}.sam_to_bam.stderr"
    params:
        assembler = config['assembly']['assembler']
    wildcard_constraints:
        sample="|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER)
    shell:
        """
        samtools view -o {output.bam} {input.sam} \
            > {log.stdout} 2> {log.stderr} \
        && rm {input.sam}
        """

rule bam_sorting:
    input:
        # reads mapped on the assembly
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.bam"
    output:
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/{assembler}/{sample}.sorting.stdout",
        stderr = "logs/05_binning/samtools/{assembler}/{sample}.sorting.stderr"
    params:
        assembler = config['assembly']['assembler']
    wildcard_constraints:
        sample="|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER)
    shell:
        """
        samtools sort -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr}
        """


# binning rules

# metabat2
rule get_contigs_depth:
    input:
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.sorted.bam"
    output:
        bam_depth_matrix = "results/05_binning/metabat2/{assembler}/{sample}.depth_matrix.tab"
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{assembler}/{sample}.depth_matrix.stdout",
        stderr = "logs/05_binning/metabat2/{assembler}/{sample}.depth_matrix.stderr"
    params:
        assembler = config['assembly']['assembler']
    wildcard_constraints:
        sample="|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER)
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.bam_depth_matrix} \
            {input.bam} > {log.stdout} 2> {log.stderr}
        """

rule metabat2_binning:
    input:
        assembly = "results/03_assembly/{assembler}/{sample}/assembly.fa.gz",
        bam_depth_matrix = "results/05_binning/metabat2/{assembler}/{sample}.depth_matrix.tab"
    output:
        output = directory("results/05_binning/metabat2/bins/{assembler}/{sample}")
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{assembler}/{sample}.binning.stdout",
        stderr = "logs/05_binning/metabat2/{assembler}/{sample}.binning.stderr",
        stdout_gz = "logs/05_binning/metabat2/{assembler}/{sample}.gzipping.stdout",
        stderr_gz = "logs/05_binning/metabat2/{assembler}/{sample}.gzipping.stderr",
        stdout_mk = "logs/05_binning/metabat2/{assembler}/{sample}.mkdir.stdout",
        stderr_mk = "logs/05_binning/metabat2/{assembler}/{sample}.mkdir.stderr",
        stdout_mv = "logs/05_binning/metabat2/{assembler}/{sample}.moving.stdout",
        stderr_mv = "logs/05_binning/metabat2/{assembler}/{sample}.moving.stderr"
    params:
        min_contig_size = config['binning']['metabat2']['min_contig_size'],
        minimum_mean_coverage = config['binning']['metabat2']['minimum_mean_coverage'],
        min_bin_size = config['binning']['metabat2']['min_bin_size'],
        bin_basename = "{sample}",
        assembler = config['assembly']['assembler']
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER)
    threads: config['binning']['metabat2']['threads']
    shell:
        """
        metabat2 -i {input.assembly} -o "{output.output}/{params.bin_basename}" \
            --abdFile {input.bam_depth_matrix} \
            --minContig {params.min_contig_size} \
            --minCV {params.minimum_mean_coverage} \
            --minClsSize {params.min_bin_size} \
            --numThreads {threads} \
            --verbose \
            > {log.stdout} 2> {log.stderr} \
        && \
        pigz --verbose {output.output}/* > {log.stdout_gz} 2> {log.stderr_gz} \
        && \
        mkdir --verbose {output.output}/bins > {log.stdout_mk} 2> {log.stderr_mk} \
        && \
        mv --verbose {output.output}/*.gz {output.output}/bins > {log.stdout_mv} 2> {log.stderr_mv} \
        """

# semibin2
rule semibin2_binning:
    input:
        assembly = "results/03_assembly/{assembler}/{sample}/assembly.fa.gz",
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.sorted.bam"
    output:
        output = directory("results/05_binning/semibin2/bins/{assembler}/{sample}")
    conda:
        "../envs/semibin.yaml"
    log:
        stdout = "logs/05_binning/semibin2/{assembler}/{sample}.binning.stdout",
        stderr = "logs/05_binning/semibin2/{assembler}/{sample}.binning.stderr",
        stdout_move = "logs/05_binning/semibin2/{assembler}/{sample}.move.stdout",
        stderr_move = "logs/05_binning/semibin2/{assembler}/{sample}.move.stderr"
    params:
        environment = config['binning']['semibin2']['environment'],
        assembler = config['assembly']['assembler']
    threads: config['binning']['semibin2']['threads']
    shell:
        """
        SemiBin2 single_easy_bin \
                --environment {params.environment} \
                -i {input.assembly} \
                -b {input.bam} \
                -o {output.output} \
                --threads {threads} \
                --verbose \
            > {log.stdout} 2> {log.stderr} \
        && \
        mkdir {output.output}/bins \
        && \
        mv --verbose {output.output}/output_bins/* {output.output}/bins \
            > {log.stdout_move} 2> {log.stderr_move}
        """

# vamb
rule vamb_binning:
    input:
        assembly = "results/03_assembly/{assembler}/{sample}/assembly.fa.gz",
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.sorted.bam"
    output:
        output = directory("results/05_binning/vamb/bins/{assembler}/{sample}")
    conda:
        "../envs/vamb.yaml"
    log:
        stdout = "logs/05_binning/vamb/{assembler}/{sample}.binning.stdout",
        stderr = "logs/05_binning/vamb/{assembler}/{sample}.binning.stderr",
    params:
        minfasta = config['binning']['vamb']['minfasta'],
        gpu = config['binning']['vamb']['gpu'],
        epochs = config['binning']['vamb']['epochs'],
        batch_sizes = config['binning']['vamb']['batch_sizes'],
        start_batch_size = config['binning']['vamb']['start_batch_size'],
        assembler = config['assembly']['assembler'],
    threads: config['binning']['vamb']['threads']
    shell:
        """
        vamb --outdir {output.output} \
            --fasta {input.assembly} \
            --bamfiles {input.bam} \
            --minfasta {params.minfasta} \
            -e {params.epochs} \
            -p {threads} \
            -q {params.batch_sizes} \
            -t {params.start_batch_size} \
            {params.gpu} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mkdir -p {output.output}/vamb_files \
        && \
        mv {output.output}/*.npz {output.output}/*.txt {output.output}/*.pt {output.output}/*.tsv {output.output}/vamb_files \
        && \
        pigz --verbose {output.output}/bins/* \
        """