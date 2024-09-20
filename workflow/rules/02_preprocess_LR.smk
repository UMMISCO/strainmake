rule fastp_long_read:
    input: lambda wildcards: get_fastq_long_read(SAMPLES_DF, wildcards.sample_lr)
    output:
        r1 = "results/02_preprocess/fastp_long_read/{sample_lr}_1.fastq.gz",
        html_report = "results/02_preprocess/fastp_long_read/{sample_lr}_report.html",
        json_report = "results/02_preprocess/fastp_long_read/{sample_lr}_report.json"
    conda: 
        "../envs/fastp.yaml"
    log:
        stdout = "logs/02_preprocess/fastp_long_read/{sample_lr}.stdout",
        stderr = "logs/02_preprocess/fastp_long_read/{sample_lr}.stderr"
    benchmark:
        "benchmarks/02_preprocess/fastp_long_read/{sample_lr}.benchmark.txt"
    params:
        compression_level = config['fastp_long_read']['compression'],
        min_phred = config['fastp_long_read']['qualified_quality_phred'],
        min_read_length = config['fastp_long_read']['minimal_read_length']
    shell:
        """
        fastp -i {input} -o {output.r1} \
            --length_required {params.min_read_length} \
            --qualified_quality_phred {params.min_phred} \
            --compression {params.compression_level} \
            --json {output.json_report} \
            --html {output.html_report} \
            > {log.stdout} 2> {log.stderr}
        """