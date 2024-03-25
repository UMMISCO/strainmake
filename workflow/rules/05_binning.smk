SAMPLES = config['samples']

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
        threads = config['binning']['bowtie2']['threads'],
        seed = config['binning']['bowtie2']['seed'],
        index_basename = "{sample}",
        assembler = config['assembly']['assembler']
    shell:
        """
        mkdir -p {output} \
        && \
        bowtie2-build --threads {params.threads} --seed {params.seed} \
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
        threads = config['binning']['bowtie2']['threads'],
        index_basename = "{sample}",
        assembler = config['assembly']['assembler']
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
        sam = "results/05_binning/bowtie2/{assembler}/{sample}.sam"
    output:
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        stdout = "logs/05_binning/samtools/{assembler}/{sample}.sam_to_bam.stdout",
        stderr = "logs/05_binning/samtools/{assembler}/{sample}.sam_to_bam.stderr"
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
        bam = "results/05_binning/bowtie2/{assembler}/{sample}.sorted.bam"
    output:
        bam_depth_matrix = "results/05_binning/metabat2/{assembler}/{sample}.depth_matrix.tab"
    conda:
        "../envs/metabat.yaml"
    log:
        stdout = "logs/05_binning/metabat2/{assembler}/{sample}.depth_matrix.stdout",
        stderr = "logs/05_binning/metabat2/{assembler}/{sample}.depth_matrix.stderr"
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
        stderr_gz = "logs/05_binning/metabat2/{assembler}/{sample}.gzipping.stderr"
    params:
        min_contig_size = config['binning']['metabat2']['min_contig_size'],
        minimum_mean_coverage = config['binning']['metabat2']['minimum_mean_coverage'],
        min_bin_size = config['binning']['metabat2']['min_bin_size'],
        threads = config['binning']['metabat2']['threads'],
        bin_basename = "{sample}",
        assembler = config['assembly']['assembler']
    shell:
        """
        metabat2 -i {input.assembly} -o "{output.output}/{params.bin_basename}" \
            --abdFile {input.bam_depth_matrix} \
            --minContig {params.min_contig_size} \
            --minCV {params.minimum_mean_coverage} \
            --minClsSize {params.min_bin_size} \
            --numThreads {params.threads} \
            --verbose \
            > {log.stdout} 2> {log.stderr} \
        && \
        pigz --verbose {output}/* > {log.stdout_gz} 2> {log.stderr_gz}
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
        threads = config['binning']['semibin2']['threads'],
        assembler = config['assembly']['assembler']
    shell:
        """
        SemiBin2 single_easy_bin \
                --environment {params.environment} \
                -i {input.assembly} \
                -b {input.bam} \
                -o {output.output} \
                --threads {params.threads} \
                --verbose \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv --verbose {output.output}/output_bins/* {output.output} \
            > {log.stdout_move} 2> {log.stderr_move}
        """

# basalt binning
rule basalt_binning:
    input:
        assembly = expand("results/03_assembly/{assembler}/{sample}/assembly.fa.gz",
                          assembler=config['assembly']['assembler'],
                          sample=config['samples']),
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        output = touch("results/05_binning/basalt/bins/{sample}")
    conda:
        "../envs/basalt.yaml"
    log:
        stdout = "logs/05_binning/basalt/{sample}.binning.stdout",
        stderr = "logs/05_binning/basalt/{sample}.binning.stderr",
        stdout_unzip = "logs/05_binning/basalt/{sample}.unzip.stdout",
        stderr_unzip = "logs/05_binning/basalt/{sample}.unzip.stderr",
    params:
        threads = config['binning']['basalt']['threads'],
        memory = config['binning']['basalt']['memory'],
        # joining the assemblies' path with using a comma
        assemblies = lambda wildcards, input: ','.join(input.assembly),
        basalt_zip_path = "BASALT_script.zip",
        basalt_exe_path = "BASALT_script/BASALT",
        r1_copy = "r1.fa.gz",
        r2_copy = "r2.fa.gz"
    shell:
        # https://github.com/EMBL-PKU/BASALT/issues/11
        """
        wget -O {params.basalt_zip_path} \
            https://github.com/EMBL-PKU/BASALT/raw/master/BASALT_script.zip \
        && \
        unzip -o {params.basalt_zip_path} \
            > {log.stdout_unzip} 2> {log.stderr} \
        && \
        chmod 777 {params.basalt_exe_path} \
        && \
        mkdir -p /root/.cache/BASALT \
        && \
        cp {input.r1} {params.r1_copy} \
        && \
        cp {input.r2} {params.r2_copy} \
        && \
        ./{params.basalt_exe_path} --assemblies {params.assemblies} \
            --shortreads '{params.r1_copy}'/'{params.r2_copy}' \
            --threads {params.threads} \
            --ram {params.memory} \
            > {log.stdout} 2> {log.stderr}
        """