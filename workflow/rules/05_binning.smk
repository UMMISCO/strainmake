# no matter the assembler used next is the assembly name to obtain
assembly_name = "{sample}.final_assembly.fasta"

SAMPLES = config['samples']

# rules for mapping reads on assembly and producing a .SAM file
rule bowtie_assembly_index:
    input:
        # assembly produced in step 03
        assembly = f"results/03_assembly/assembly/{assembly_name}.gz"
    output:
        directory("results/05_binning/bowtie2/index/{sample}")
    conda:
        "../envs/bowtie2.yaml"
    log:
        stdout = "logs/05_binning/bowtie2/{sample}.indexing.stdout",
        stderr = "logs/05_binning/bowtie2/{sample}.indexing.stdout"
    params:
        threads = config['binning']['bowtie2']['threads'],
        seed = config['binning']['bowtie2']['seed'],
        index_basename = "{sample}"
    shell:
        """
        mkdir -p {output} \
        && \
        bowtie2-build --threads {params.threads} --seed {params.seed} \
            {input.assembly} "{output}/{params.index_basename}" \
            > {log.stdout} 2> {log.stderr}
        """

rule reads_mapping:
    input:
        # metagenome reads
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz",
        # assembly index
        bowtie_index = expand("results/05_binning/bowtie2/index/{sample}",
                              sample=SAMPLES)
    output:
        sam = "results/05_binning/bowtie2/{sample}.sam"
    conda:
        "../envs/bowtie2.yaml"
    log:
        stdout = "logs/05_binning/bowtie2/{sample}.mapping.stdout",
        stderr = "logs/05_binning/bowtie2/{sample}.mapping.stdout"
    params:
        threads = config['binning']['bowtie2']['threads'],
        index_basename = "{sample}"
    shell:
        """
        bowtie2 -p {params.threads} \
            -x "{input.bowtie_index}/{params.index_basename}" \
            -1 {input.r1} -2 {input.r2} \
            -S {output.sam} \
            > {log.stdout} 2> {log.stderr}
        """

rule sam_to_bam:
    input:
        sam = "results/05_binning/bowtie2/{sample}.sam"
    output:
        bam = "results/05_binning/bowtie2/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/{sample}.sam_to_bam.stdout",
        stderr = "logs/05_binning/samtools/{sample}.sam_to_bam.stdout"
    shell:
        """
        samtools view -o {output.bam} {input.sam} \
            > {log.stdout} 2> {log.stderr} \
        && rm {input.sam}
        """

rule bam_sorting:
    input:
        # reads mapped on the assembly
        bam = "results/05_binning/bowtie2/{sample}.bam"
    output:
        bam = "results/05_binning/bowtie2/{sample}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/{sample}.sorting.stdout",
        stderr = "logs/05_binning/samtools/{sample}.sorting.stderr"
    shell:
        """
        samtools sort -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr} \
        && \
            rm {input.bam}
        """

# binning rules

# metabat2
rule get_contigs_depth:
    input:
        bam = "results/05_binning/bowtie2/{sample}.sorted.bam"
    output:
        bam_depth_matrix = "results/05_binning/metabat2/{sample}.depth_matrix.tab"
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{sample}.depth_matrix.stdout",
        stderr = "logs/05_binning/metabat2/{sample}.depth_matrix.stderr"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.bam_depth_matrix} \
            {input.bam} > {log.stdout} 2> {log.stderr}
        """

rule metabat2_binning:
    input:
        assembly = f"results/03_assembly/assembly/{assembly_name}.gz",
        bam_depth_matrix = "results/05_binning/metabat2/{sample}.depth_matrix.tab"
    output:
        touch("results/05_binning/{sample}_binning"),
        output = directory("results/05_binning/metabat2/bins/{sample}/")
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{sample}.binning.stdout",
        stderr = "logs/05_binning/metabat2/{sample}.binning.stderr"
    params:
        min_contig_size = config['binning']['metabat2']['min_contig_size'],
        minimum_mean_coverage = config['binning']['metabat2']['minimum_mean_coverage'],
        min_bin_size = config['binning']['metabat2']['min_bin_size'],
        threads = config['binning']['metabat2']['threads'],
        bin_basename = "{sample}"
    shell:
        """
        metabat2 -i {input.assembly} -o "{output.output}/{params.bin_basename}" \
            --abdFile {input.bam_depth_matrix} \
            --minContig {params.min_contig_size} \
            --minCV {params.minimum_mean_coverage} \
            --minClsSize {params.min_bin_size} \
            --numThreads {params.threads} \
            --verbose \
            > {log.stdout} 2> {log.stderr}
        """

# copy bins into a results folder
rule metabat2_move_bins:
    input:
        expand("results/05_binning/metabat2/bins/{sample}", sample=SAMPLES)
    output:
        directory("results/05_binning/bins/{sample}")
    log:
        stdout = "logs/05_binning/bins/{sample}.stdout",
        stderr = "logs/05_binning/bins/{sample}.stderr"
    shell:
       "mkdir -p {output} && cp {input}/*.fa {output}/"