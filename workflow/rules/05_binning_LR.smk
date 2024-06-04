SAMPLES_TABLE = config['samples']
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)

# using minimap2 to map the long reads on the assembly
rule reads_mapping_LR:
    input:
        # metagenome reads
        long_read = "results/02_preprocess/bowtie2/{sample_lr}_1.clean.fastq.gz",
        # assembly
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz"
    output:
        sam = "results/05_binning/minimap2/{assembler_lr}/{sample_lr}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        stdout = "logs/05_binning/minimap2/{assembler_lr}/{sample_lr}.mapping.stdout",
        stderr = "logs/05_binning/minimap2/{assembler_lr}/{sample_lr}.mapping.stderr"
    params:
        assembler_lr = config['assembly']['long_read_assembler'],
        method = "map-ont" if config['assembly']['metaflye']['method'] == "nanopore" else "map-pb"
    wildcard_constraints:
        assembler = "|".join("{params.assembler}")    
    threads: config['binning']['minimap2']['threads']
    shell:
        """
        minimap2 -ax {params.method} -t {threads} \
            {input.assembly} \
            {input.long_read} \
            > {output.sam}
        """

rule sam_to_bam_LR:
    input:
        sam = "results/05_binning/minimap2/{assembler_lr}/{sample_lr}.sam"
    output:
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sam_to_bam.stdout",
        stderr = "logs/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sam_to_bam.stderr"
    params:
        assembler = config['assembly']['long_read_assembler']
    wildcard_constraints:
        sample_lr="|".join(SAMPLES_LR),
        assembler = "|".join("{params.assembler}")
    shell:
        """
        samtools view -o {output.bam} {input.sam} \
            > {log.stdout} 2> {log.stderr} \
        && rm {input.sam}
        """

rule bam_sorting_LR:
    input:
        # reads mapped on the assembly
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.bam"
    output:
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sorting.stdout",
        stderr = "logs/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sorting.stderr"
    params:
        assembler = config['assembly']['long_read_assembler']    
    wildcard_constraints:
        sample_lr="|".join(SAMPLES_LR),
        assembler = "|".join("{params.assembler}")
    shell:
        """
        samtools sort -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr}
        """

# sorting BAM by read name (for VAMB)
rule bam_sorting_by_readname_LR:
    input:
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sorted.bam"
    output:
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sorted_by_readname.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/{assembler_lr}/{sample_lr}.sorting_by_readname.stdout",
        stderr = "logs/05_binning/samtools/{assembler_lr}/{sample_lr}.sorting_by_readname.stderr"
    shell:
        """
        samtools sort -n -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr}
        """
        
# semibin2
rule semibin2_binning_LR:
    input:
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz",
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sorted.bam"
    output:
        output = directory("results/05_binning/LR/semibin2/bins/{assembler_lr}/{sample_lr}")
    conda:
        "../envs/semibin.yaml"
    log:
        stdout = "logs/05_binning/semibin2/{assembler_lr}/{sample_lr}.binning.stdout",
        stderr = "logs/05_binning/semibin2/{assembler_lr}/{sample_lr}.binning.stderr",
        stdout_move = "logs/05_binning/semibin2/{assembler_lr}/{sample_lr}.move.stdout",
        stderr_move = "logs/05_binning/semibin2/{assembler_lr}/{sample_lr}.move.stderr"
    params:
        environment = config['binning']['semibin2']['environment'],
        assembler_lr = config['assembly']['long_read_assembler']
    threads: config['binning']['semibin2']['threads']
    shell:
        """
        SemiBin2 single_easy_bin \
                --environment {params.environment} \
                -i {input.assembly} \
                -b {input.bam} \
                -o {output.output} \
                --threads {threads} \
                --sequencing-type=long_read \
                --verbose \
            > {log.stdout} 2> {log.stderr} \
        && \
        mkdir {output.output}/bins \
        && \
        mv --verbose {output.output}/output_bins/* {output.output}/bins \
            > {log.stdout_move} 2> {log.stderr_move}
        """