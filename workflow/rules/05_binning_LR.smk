SAMPLES_TABLE = config['samples']
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)
ASSEMBLER_LR = config['assembly']['long_read_assembler']

# taking into account the case where we don't have LR
if ASSEMBLER_LR == None:
       ASSEMBLER_LR = []

# using minimap2 to map the long reads on the assembly
rule reads_mapping_LR:
    input:
        # metagenome reads
        long_read = "results/02_preprocess/bowtie2/{sample_lr}_1.clean.fastq.gz",
        # assembly
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz"
    output:
        sam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        stdout = "logs/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.mapping.stdout",
        stderr = "logs/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.mapping.stderr"
    benchmark:
        "benchmarks/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.mapping.benchmark.txt"
    params:
        method = "map-ont" if config['assembly']['metaflye']['method'] == "nanopore" else "map-pb"
    wildcard_constraints:
        assembler = "|".join(ASSEMBLER_LR)    
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
        sam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sam"
    output:
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sam_to_bam.stdout",
        stderr = "logs/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sam_to_bam.stderr"
    benchmark:
        "benchmarks/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sam_to_bam.benchmark.txt"
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
    benchmark:
        "benchmarks/05_binning/samtools/LR/{assembler_lr}/{sample_lr}.sorting.benchmark.txt"
    params:
        assembler = config['assembly']['long_read_assembler']    
    wildcard_constraints:
        sample_lr="|".join(SAMPLES_LR),
        assembler = "|".join("{params.assembler}")
    shell:
        """
        samtools sort -o {output.bam} {input.bam} \
            > {log.stdout} 2> {log.stderr} \
        && \
        rm {input.bam}
        """
        
# binning rules

# metabat2
rule get_contigs_depth_LR:
    input:
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sorted.bam"
    output:
        bam_depth_matrix = "results/05_binning/LR/metabat2/{assembler_lr}/{sample_lr}.depth_matrix.tab"
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.depth_matrix.stdout",
        stderr = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.depth_matrix.stderr"
    benchmark:
        "benchmarks/05_binning/metabat2/{assembler_lr}/{sample_lr}.depth_matrix.benchmark.txt"
    wildcard_constraints:
        sample_lr = "|".join(SAMPLES_LR),
        assembler_lr   = "|".join(ASSEMBLER_LR)
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.bam_depth_matrix} \
            {input.bam} > {log.stdout} 2> {log.stderr}
        """

rule metabat2_binning_LR:
    input:
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz",
        bam_depth_matrix = "results/05_binning/LR/metabat2/{assembler_lr}/{sample_lr}.depth_matrix.tab"
    output:
        output = directory("results/05_binning/LR/metabat2/bins/{assembler_lr}/{sample_lr}")
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.binning.stdout",
        stderr = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.binning.stderr",
        stdout_gz = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.gzipping.stdout",
        stderr_gz = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.gzipping.stderr",
        stdout_mk = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.mkdir.stdout",
        stderr_mk = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.mkdir.stderr",
        stdout_mv = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.moving.stdout",
        stderr_mv = "logs/05_binning/metabat2/{assembler_lr}/{sample_lr}.moving.stderr"
    benchmark:
        "benchmarks/05_binning/metabat2/{assembler_lr}/{sample_lr}.binning.benchmark.txt"
    params:
        min_contig_size = config['binning']['metabat2']['min_contig_size'],
        minimum_mean_coverage = config['binning']['metabat2']['minimum_mean_coverage'],
        min_bin_size = config['binning']['metabat2']['min_bin_size'],
        bin_basename = "{sample_lr}"
    wildcard_constraints:
        assembler_lr = "|".join(ASSEMBLER_LR)
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
    benchmark:
        "benchmarks/05_binning/semibin2/{assembler_lr}/{sample_lr}.binning.benchmark.txt"
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

# vamb
rule vamb_binning_LR:
    input:
        assembly = "results/03_assembly/LR/{assembler_lr}/{sample_lr}/assembly.fa.gz",
        bam = "results/05_binning/minimap2/LR/{assembler_lr}/{sample_lr}.sorted.bam"
    output:
        output = directory("results/05_binning/LR/vamb/bins/{assembler_lr}/{sample_lr}")
    conda:
        "../envs/vamb.yaml"
    log:
        stdout = "logs/05_binning/vamb/{assembler_lr}/{sample_lr}.binning.stdout",
        stderr = "logs/05_binning/vamb/{assembler_lr}/{sample_lr}.binning.stderr",
    benchmark:
        "benchmarks/05_binning/vamb/{assembler_lr}/{sample_lr}.binning.benchmark.txt"
    params:
        minfasta = config['binning']['vamb']['minfasta'],
        gpu = config['binning']['vamb']['gpu'],
        epochs = config['binning']['vamb']['epochs'],
        batch_sizes = config['binning']['vamb']['batch_sizes'],
        start_batch_size = config['binning']['vamb']['start_batch_size'],
        assembler_lr = config['assembly']['assembler'],
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
        && \
        for file in {output.output}/bins/*.fna.gz; do mv "$file" "${{file%.fna.gz}}.fa.gz"; done
        """