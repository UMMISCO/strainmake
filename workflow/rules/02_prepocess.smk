rule fastp:
    input: lambda wildcards: get_fastq_pair(SAMPLES_DF, wildcards.sample)
    output:
        r1 = "results/02_preprocess/fastp/{sample}_1.fastq.gz",
        r2 = "results/02_preprocess/fastp/{sample}_2.fastq.gz",
        html_report = "results/02_preprocess/fastp/{sample}_report.html",
        json_report = "results/02_preprocess/fastp/{sample}_report.json"
    conda: 
        "../envs/fastp.yaml"
    log:
        stdout = "logs/02_preprocess/fastp/{sample}.stdout",
        stderr = "logs/02_preprocess/fastp/{sample}.stderr"
    benchmark:
        "benchmarks/02_preprocess/fastp/{sample}.benchmark.txt"
    params:
        compression_level = config['fastp']['compression'],
        min_phred = config['fastp']['qualified_quality_phred'],
        min_read_length = config['fastp']['minimal_read_length']
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            --length_required {params.min_read_length} \
            --qualified_quality_phred {params.min_phred} \
            --compression {params.compression_level} \
            --json {output.json_report} \
            --html {output.html_report} \
            > {log.stdout} 2> {log.stderr}
        """

# keeps reads that don't map on human genome to decontaminate the metagenome
rule host_decontamination:
    input:
        r1 = "results/02_preprocess/fastp/{sample}_1.fastq.gz",
        r2 = "results/02_preprocess/fastp/{sample}_2.fastq.gz",
        index = "results/02_preprocess/bowtie2/index"
    output:
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    conda:
        "../envs/bowtie2.yaml"
    log:
        stderr = "logs/02_preprocess/bowtie2/{sample}.stderr"
    benchmark:
        "benchmarks/02_preprocess/bowtie2/{sample}.benchmark.txt"
    params:
        organism_name = config['bowtie2']['index_name'],
        bowtie_output_name = "results/02_preprocess/bowtie2/{sample}_%.clean.fastq.gz"
    threads: config['bowtie2']['threads']
    shell:
        """
        bowtie2 -p {threads} -x "{input.index}/{params.organism_name}" \
            -1 {input.r1} -2 {input.r2} \
            --un-conc-gz {params.bowtie_output_name} \
            > /dev/null 2> {log.stderr}
        """

# get the bowtie2 index from the internet
rule get_bowtie_index:
    output:
        directory("results/02_preprocess/bowtie2/index")
    log:
        wget_stdout = "logs/02_preprocess/bowtie2/get_bowtie_index.wget.stdout",
        wget_stderr = "logs/02_preprocess/bowtie2/get_bowtie_index.wget.stderr",
        unzip_stdout = "logs/02_preprocess/bowtie2/get_bowtie_index.unzip.stdout",
        unzip_stderr = "logs/02_preprocess/bowtie2/get_bowtie_index.unzip.stderr",
        mv_stdout = "logs/02_preprocess/bowtie2/get_bowtie_index.mv.stdout",
        mv_stderr = "logs/02_preprocess/bowtie2/get_bowtie_index.mv.stderr",
        rm_stdout = "logs/02_preprocess/bowtie2/get_bowtie_index.rm.stdout",
        rm_stderr = "logs/02_preprocess/bowtie2/get_bowtie_index.rm.stderr"
    conda:
        "../envs/get_bowtie_index.yaml"
    params:
        organism_name = config['bowtie2']['index_name']
    benchmark:
        "benchmarks/02_preprocess/bowtie2/get_bowtie_index.benchmark.txt"
    shell:
    # downloading the already made index file and unzipping it. The indexes 
    # will be in a folder of name "index"
        """
        mkdir -p results/02_preprocess/bowtie2 \
        && wget -O results/02_preprocess/bowtie2/{params.organism_name}.zip \
            https://genome-idx.s3.amazonaws.com/bt/{params.organism_name}.zip \
            > {log.wget_stdout} 2> {log.wget_stderr} \
        && unzip -d results/02_preprocess/bowtie2/index \
            results/02_preprocess/bowtie2/{params.organism_name}.zip \
            > {log.unzip_stdout} 2> {log.unzip_stderr} \
        && mv results/02_preprocess/bowtie2/index/{params.organism_name}/* \
            results/02_preprocess/bowtie2/index/ \
            > {log.mv_stdout} 2> {log.mv_stderr} \
        && rm -r results/02_preprocess/bowtie2/index/{params.organism_name} \
            > {log.rm_stdout} 2> {log.rm_stderr}
        """

rule fastqc_after_preprocessing:
    input:
        "results/02_preprocess/bowtie2/{sample}_{read}.clean.fastq.gz"
    output:
        html_report="results/02_preprocess/fastqc/{sample}_{read}.clean_fastqc.html",
        zip_report="results/02_preprocess/fastqc/{sample}_{read}.clean_fastqc.zip"
    conda: 
        "../envs/fastqc.yaml"
    log:
        stdout = "logs/02_preprocess/fastqc/{sample}_{read}.clean.stdout",
        stderr = "logs/02_preprocess/fastqc/{sample}_{read}.clean.stderr"
    benchmark:
        "benchmarks/02_preprocess/fastqc/{sample}.{read}.benchmark.txt"
    shell:
        """
        fastqc {input} -o results/02_preprocess/fastqc/ \
            > {log.stdout} 2> {log.stderr}
        """

# allows a flexibility for the user to use sequences in FASTA or FASTQ format
seq_format = config["lr_seq_format"]
sequences_file_end = f"_1.{seq_format}.gz"

rule downsize_reads_for_hybrid_parts_sr:
    input:
        sr_1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        sr_2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz",
        lr = "results/02_preprocess/fastp_long_read/{sample}" + sequences_file_end
    output:
        sr_1 = "results/02_preprocess/downsized/bowtie2/{sample}_1.clean.downsized.fastq.gz",
        sr_2 = "results/02_preprocess/downsized/bowtie2/{sample}_2.clean.downsized.fastq.gz",
        lr = "results/02_preprocess/downsized/fastp_long_read/{sample}_downsized" + sequences_file_end
    conda:
        "../envs/seqkit.yaml"
    log:
        stdout = "logs/02_preprocess/downsizing/{sample}.stdout",
        stderr = "logs/02_preprocess/downsizing/{sample}.stderr"
    benchmark:
        "benchmarks/02_preprocess/downsizing/{sample}.benchmark.txt"
    params:
        prop_lr = config["downsizing_for_hybrid"]["lr"],
        prop_sr = config["downsizing_for_hybrid"]["sr"]
    shell:
        """
        (
            seqkit sample -p {params.prop_sr} {input.sr_1} --rand-seed 11111 > {output.sr_1} \
            && \
            seqkit sample -p {params.prop_sr} {input.sr_2} --rand-seed 11111 > {output.sr_2} \
            && \
            seqkit sample -p {params.prop_lr} {input.lr} --rand-seed 11111 > {output.lr}
        ) > {log.stdout} 2> {log.stderr}
        """
