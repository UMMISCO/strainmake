SAMPLES = config['samples']

# download the database for CheckM2
rule checkm2_database:
    output:
        "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    conda:
        "../envs/checkm2.yaml"
    log:
        stdout = "logs/06_binning_qc/checkm2/checkm2.db.stdout",
        stderr = "logs/06_binning_qc/checkm2/checkm2.db.stdout"
    params:
        output_path = "results/06_binning_qc/checkm2/database"
    shell:
        """
        checkm2 database --download --path {params.output_path} \
            > {log.stdout} 2> {log.stderr}
        """

rule checkm2_assessment:
    input:
        # folder with bins created in step 05
        bins = expand("results/05_binning/bins/{sample}", sample=SAMPLES),
        diamond_database = "results/06_binning_qc/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        out_dir = directory("results/06_binning_qc/checkm2/{sample}")
    conda:
        "../envs/checkm2.yaml"
    log:
        stdout = "logs/06_binning_qc/checkm2/{sample}.assessment.stdout",
        stderr = "logs/06_binning_qc/checkm2/{sample}.assessment.stdout"
    params:
        threads = config['checkm2']['threads']
    shell:
        """
        echo {input.bins} \
        && \
        checkm2 predict --input {input.bins} --threads {params.threads} \
            -x .fa \
            --database_path {input.diamond_database} \
            --output-directory {output.out_dir} > {log.stdout} 2> {log.stderr}
        """