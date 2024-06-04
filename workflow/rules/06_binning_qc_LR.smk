SAMPLES_TABLE = config['samples']
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)
ASSEMBLER_LR = config['assembly']['long_read_assembler']

rule checkm2_assessment_LR:
    input:
        # folder with bins created in step 05. One folder per binning program
        bins = "results/05_binning/LR/{long_read_binner}/bins/{assembler_lr}/{sample_lr}", 
        diamond_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        "results/06_binning_qc/checkm2/LR/{long_read_binner}/{assembler_lr}/{sample_lr}/quality_report.tsv",
        out_dir = directory("results/06_binning_qc/checkm2/LR/{long_read_binner}/{assembler_lr}/{sample_lr}")
    conda:
        "../envs/checkm2.yaml"
    log:
        stdout = "logs/06_binning_qc/checkm2/{long_read_binner}/{assembler_lr}/{sample_lr}.assessment.stdout",
        stderr = "logs/06_binning_qc/checkm2/{long_read_binner}/{assembler_lr}/{sample_lr}.assessment.stderr"
    params:
        long_read_binner = config['binning']['long_read_binner'],
        assembler = config['assembly']['long_read_assembler']
    threads: config['checkm2']['threads']
    wildcard_constraints:
        sample="|".join(SAMPLES),
        assembler = "|".join(ASSEMBLER_LR)
    shell:
        """
        echo {input.bins} \
        && \
        checkm2 predict --input {input.bins}/bins --threads {threads} \
            -x .gz \
            --database_path {input.diamond_database} \
            --output-directory {output.out_dir} > {log.stdout} 2> {log.stderr}
        """    