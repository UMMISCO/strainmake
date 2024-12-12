# no matter the assembler used next is the assembly name to obtain
assembly_name = "{sample}.final_assembly.fasta"

# allows a flexibility for the user to use sequences in FASTA or FASTQ format
# (Flye can assembly sequences using FASTA or FASTQ)
seq_format = config["lr_seq_format"]
sequences_file_end = f"_1.{seq_format}.gz"

rule megahit_assembly:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        assembly = "results/03_assembly/megahit/{sample}/assembly.fa",
    conda:
        "../envs/megahit.yaml"
    log:
        stdout = "logs/03_assembly/megahit/{sample}.stdout",
        stderr = "logs/03_assembly/megahit/{sample}.stderr"
    benchmark:
        "benchmarks/03_assembly/megahit/{sample}.benchmark.txt"
    params:
        out_dir = "results/03_assembly/megahit/{sample}",
        tmp_dir = "tmp/",
        tmp_output = "{sample}_tmp_megahit_output",
        min_contig_len = config['assembly']['megahit']['min_contig_len']
    threads: config['assembly']['megahit']['threads']
    shell:
        """
        mkdir -p {params.tmp_dir} \
        && \
        megahit -1 {input.r1} -2 {input.r2} \
            --min-contig-len {params.min_contig_len} \
            --num-cpu-threads {threads} \
            --tmp-dir {params.tmp_dir} \
            --out-dir {params.tmp_output} > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.tmp_output}/* {params.out_dir} \
        && \
        rm -r {params.tmp_output} \
        && \
        mv {params.out_dir}/final.contigs.fa {params.out_dir}/assembly.fa
        """

# for VAMB for example. It will replace the spaces in the FASTA headers and gzip the asssembly
rule megahit_fasta_headers_renaming:
    input: "results/03_assembly/megahit/{sample}/assembly.fa"
    output: 
        assembly = "results/03_assembly/megahit/{sample}/assembly.fa.gz",
        other_files = "results/03_assembly/megahit/{sample}/other_files.tar.gz"
    conda:
        "../envs/pigz.yaml"
    log:
        stdout = "logs/03_assembly/megahit/{sample}.rename.stdout",
        stderr = "logs/03_assembly/megahit/{sample}.rename.stderr",
        stdout_pigz = "logs/03_assembly/megahit/{sample}.rename.pigz.stdout",
        stderr_pigz = "logs/03_assembly/megahit/{sample}.rename.pigz.stderr"
    benchmark:
        "benchmarks/03_assembly/megahit/{sample}.rename.benchmark.txt"
    params:
        rename_script = "workflow/scripts/megahit_fasta_header_rename.py",
        intermediate_output = "results/03_assembly/megahit/{sample}/assembly.new.fa",
        compressing_files_script = "workflow/scripts/compress_spades_megahit_results.sh",
        out_dir = "results/03_assembly/megahit/{sample}"
    shell:
        """
        python3 {params.rename_script} {input} {params.intermediate_output} \
        > {log.stdout} 2> {log.stderr} \
        && \
        rm {input} \
        && \
        mv {params.intermediate_output} {input} \
        && \
        pigz --verbose {input} > {log.stdout_pigz} 2> {log.stderr_pigz} \
        && \
        bash {params.compressing_files_script} {params.out_dir}
        """

rule metaspades_assembly:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        assembly = "results/03_assembly/metaspades/{sample}/assembly.fa.gz",
        other_files = "results/03_assembly/metaspades/{sample}/other_files.tar.gz"
    conda:
        "../envs/spades.yaml"
    log:
        stdout = "logs/03_assembly/metaspades/{sample}.stdout",
        stderr = "logs/03_assembly/metaspades/{sample}.stderr"
    benchmark:
        "benchmarks/03_assembly/metaspades/{sample}.benchmark.txt"
    params:
        out_dir = "results/03_assembly/metaspades/{sample}",
        memory_limit = config['assembly']['metaspades']['memory_limit'],
        min_contig_len = config['assembly']['metaspades']['min_contig_len'],
        compressing_files_script = "workflow/scripts/compress_spades_megahit_results.sh",
        intermediate_assembly = "{sample}_metaspades_tmp_assembly.fa"
    threads: config['assembly']['metaspades']['threads']
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            -o {params.out_dir} \
            -m {params.memory_limit} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.out_dir}/scaffolds.fasta {params.out_dir}/assembly.fa \
        && \
        pigz {params.out_dir}/assembly.fa \
        && \
        seqkit seq -m {params.min_contig_len} {output.assembly} > {params.intermediate_assembly} \
        && \
        pigz {params.intermediate_assembly} \
        && \
        mv {params.intermediate_assembly}.gz {output.assembly} \
        && \
        bash {params.compressing_files_script} {params.out_dir}
        """

# hybrid assembly using hybridSPADes
rule hybridspades_assembly:
    input:
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz",
        long_read = "results/02_preprocess/fastp_long_read/{sample}" + sequences_file_end
    output:
        assembly = "results/03_assembly/hybridspades/{sample}/assembly.fa.gz",
        other_files = "results/03_assembly/hybridspades/{sample}/other_files.tar.gz"
    conda:
        "../envs/spades.yaml"
    log:
        stdout = "logs/03_assembly/hybridspades/{sample}.stdout",
        stderr = "logs/03_assembly/hybridspades/{sample}.stderr"
    benchmark:
        "benchmarks/03_assembly/hybridspades/{sample}.benchmark.txt"
    params:
        out_dir = "results/03_assembly/hybridspades/{sample}",
        memory_limit = config['assembly']['hybridspades']['memory_limit'],
        method_flag = "--nanopore" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio",
        min_contig_len = config['assembly']['hybridspades']['min_contig_len'],
        compressing_files_script = "workflow/scripts/compress_spades_megahit_results.sh",
        intermediate_assembly = "{sample}_hybridspades_tmp_assembly.fa"
    threads: config['assembly']['hybridspades']['threads']
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} \
            {params.method_flag} {input.long_read} \
            --threads {threads} \
            -o {params.out_dir} \
            -m {params.memory_limit} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.out_dir}/scaffolds.fasta {params.out_dir}/assembly.fa \
        && \
        pigz {params.out_dir}/assembly.fa \
        && \
        seqkit seq -m {params.min_contig_len} {output.assembly} > {params.intermediate_assembly} \
        && \
        pigz {params.intermediate_assembly} \
        && \
        mv {params.intermediate_assembly}.gz {output.assembly} \
        && \
        bash {params.compressing_files_script} {params.out_dir}
        """

# metaFlye for long read assembly
rule metaflye_assembly:
    input:
        # files produced by fastp (should have already been decontaminated before used in the pipeline)
        long_read = "results/02_preprocess/fastp_long_read/{sample}" + sequences_file_end,
    output:
        assembly = "results/03_assembly/metaflye/{sample}/assembly.fa.gz"
    conda:
        "../envs/flye.yaml"
    log:
        stdout = "logs/03_assembly/metaflye/{sample}.stdout",
        stderr = "logs/03_assembly/metaflye/{sample}.stderr"
    benchmark:
        "benchmarks/03_assembly/metaflye/{sample}.benchmark.txt"
    params:
        out_dir = "results/03_assembly/metaflye/{sample}",
        method_flag = "--nano-hq" if config['assembly']['metaflye']['method'] == "nanopore" else "--pacbio-hifi",
        min_contig_len = config['assembly']['metaflye']['min_contig_len'],
        intermediate_assembly = "{sample}_metaflye_tmp_assembly.fa"
    threads: config['assembly']['metaflye']['threads']
    shell:
        """
        flye {params.method_flag} {input.long_read} --out-dir {params.out_dir} \
            --meta --threads {threads} \
            > {log.stdout} 2> {log.stderr} \
        && \
        mv {params.out_dir}/assembly.fasta {params.out_dir}/assembly.fa \
        && \
        pigz {params.out_dir}/assembly.fa \
        && \
        seqkit seq -m {params.min_contig_len} {output.assembly} > {params.intermediate_assembly} \
        && \
        pigz {params.intermediate_assembly} \
        && \
        mv {params.intermediate_assembly}.gz {output.assembly}
        """