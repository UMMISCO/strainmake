# no matter the assembler used next is the assembly name to obtain
assembly_name = "{sample}.final_assembly.fasta"

rule megahit_assembly:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        assembly = "results/03_assembly/megahit/{sample}/assembly.fa"
    conda:
        "../envs/megahit.yaml"
    log:
        stdout = "logs/03_assembly/megahit/{sample}.stdout",
        stderr = "logs/03_assembly/megahit/{sample}.stdout"
    params:
        out_dir = "results/03_assembly/megahit/{sample}",
        threads = config['assembly']['megahit']['threads'],
        tmp_dir = "tmp/",
        tmp_output = "{sample}_tmp_megahit_output/"
    shell:
        """
        mkdir -p {params.tmp_dir} \
        && \
        megahit -1 {input.r1} -2 {input.r2} \
            --num-cpu-threads {params.threads} \
            --tmp-dir {params.tmp_dir} \
            --out-dir {params.tmp_output} > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.tmp_output}/* {params.out_dir} \
        && \
        rm -r {params.tmp_output} \
        && \
        mv {params.out_dir}/final.contigs.fa {params.out_dir}/assembly.fa \
        """

# for VAMB for example. It will replace the spaces in the headers an gzip the asssembly
rule megahit_fasta_headers_renaming:
    input: "results/03_assembly/megahit/{sample}/assembly.fa"
    output: "results/03_assembly/megahit/{sample}/assembly.fa.gz"
    conda:
        "../envs/pigz.yaml"
    log:
        stdout = "logs/03_assembly/megahit/{sample}.rename.stdout",
        stderr = "logs/03_assembly/megahit/{sample}.rename.stdout",
        stdout_pigz = "logs/03_assembly/megahit/{sample}.rename.pigz.stdout",
        stderr_pigz = "logs/03_assembly/megahit/{sample}.rename.pigz.stdout"
    params:
        rename_script = "workflow/scripts/megahit_fasta_header_rename.py",
        intermediate_output = "results/03_assembly/megahit/{sample}/assembly.new.fa"
    shell:
        """
        python3 {params.rename_script} {input} {params.intermediate_output} \
        > {log.stdout} 2> {log.stderr} \
        && \
        rm {input} \
        && \
        mv {params.intermediate_output} {input} \
        && \
        pigz --verbose {input} > {log.stdout_pigz} 2> {log.stderr_pigz}
        """

rule metaspades_assembly:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        assembly = "results/03_assembly/metaspades/{sample}/assembly.fa.gz"
    conda:
        "../envs/spades.yaml"
    log:
        stdout = "logs/03_assembly/metaspades/{sample}.stdout",
        stderr = "logs/03_assembly/metaspades/{sample}.stdout"
    params:
        out_dir = "results/03_assembly/metaspades/{sample}",
        threads = config['assembly']['metaspades']['threads'],
        memory_limit = config['assembly']['metaspades']['memory_limit']
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} \
            --threads {params.threads} \
            -o {params.out_dir} \
            -m {params.memory_limit} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.out_dir}/scaffolds.fasta {params.out_dir}/assembly.fa \
        && \
        pigz {params.out_dir}/assembly.fa
        """