import os
from utils import convert_to_si_units, convert_from_si_units_to_int

rule metaphlan_profiling:
    input:
        reads = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz", # on short reads only
        database = METAPHLAN_DB_MARKER
    output: "results/09_taxonomic_profiling/metaphlan/{sample}.profile.txt"
    conda:
        "../envs/metaphlan.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/metaphlan/{sample}.profile.stdout",
        stderr = "logs/09_taxonomic_profiling/metaphlan/{sample}.profile.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/metaphlan/{sample}.profile.benchmark.txt"
    params:
        database_dir = METAPHLAN_DB
    threads: config.get('taxonomic_profiling', {}).get('metaphlan', {}).get('threads', 0)
    shell:
    # --db_dir is mandatory here: left to itself MetaPhlAn downloads its markers
    # into its own conda environment, which needs network access and writes into
    # an environment that may be shared with other projects
        """
        metaphlan --input_type fastq --nproc {threads} \
            --db_dir {params.database_dir} \
            {input.reads} {output} \
        > {log.stdout} 2> {log.stderr}
        """

# we create a folder with the short reads METEOR will work on (we use symbolic links)
rule meteor_prepare_fastq:
    input:
        # files produced by fastp and decontaminated using bowtie2
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz",
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz"
    output:
        # same fastq files, but in another folder
        r1 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_1.fastq.gz",
        r2 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_2.fastq.gz"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_prepare_fastq_folder.stdout",
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.prepare_fastq_folder.benchmark.txt"
    shell:
        """
        ln -sv $(realpath {input.r1}) {output.r1} > {log} 2>&1
        ln -sv $(realpath {input.r2}) {output.r2} >> {log} 2>&1
        """

rule meteor_fastq_indexing:
    input:
        r1 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_1.fastq.gz",
        r2 = "results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_2.fastq.gz"
    output:
        # same fastq files, but in another folder
        directory("results/09_taxonomic_profiling/meteor/{sample}/fastq_index")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_fastq_indexing.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_fastq_indexing.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.fastq_indexing.benchmark.txt"
    params:
        fastq_data = lambda wildcards: os.path.dirname("results/09_taxonomic_profiling/meteor/{sample}/fastq/{sample}_1.fastq.gz".format(sample=wildcards.sample))
    shell:
        """
        meteor fastq -i $(realpath {params.fastq_data}) -p -o {output} > {log.stdout} 2> {log.stderr}
        """

rule meteor_mapping:
    input:
       fastq_index = "results/09_taxonomic_profiling/meteor/{sample}/fastq_index",
       # the catalogue is only needed here: the profiling rules that also read it
       # already depend on this rule's output
       database = METEOR_DB_MARKER
    output:
        directory("results/09_taxonomic_profiling/meteor/{sample}/mapping")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_mapping.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_mapping.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.fastq_mapping.benchmark.txt"
    threads:
        config.get('taxonomic_profiling', {}).get('meteor', {}).get('threads', 0)
    params:
        indexed_fastq_file_with_sample = lambda wildcards: os.path.join(f"results/09_taxonomic_profiling/meteor/{wildcards.sample}/fastq_index", f"{wildcards.sample}"),
        reference = METEOR_REFERENCE
    shell:
        """
        meteor mapping -i {params.indexed_fastq_file_with_sample} -o {output} -r {params.reference} \
            --trim 0 --ka --kf -t {threads} \
            > {log.stdout} 2> {log.stderr} 
        """

rule meteor_profiling:
    input:
        "results/09_taxonomic_profiling/meteor/{sample}/mapping"
    output:
        directory("results/09_taxonomic_profiling/meteor/{sample}/profiling")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_profiling.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_profiling.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.profile.benchmark.txt"
    params:
        mapping_with_sample = lambda wildcards: os.path.join(f"results/09_taxonomic_profiling/meteor/{wildcards.sample}/mapping", f"{wildcards.sample}"),
        reference = METEOR_REFERENCE
    shell:
        """ 
        meteor profile -i {params.mapping_with_sample} -o {output} -r {params.reference} \
            -n coverage \
            --seed 100 \
            > {log.stdout} 2> {log.stderr}
        """

# taxonomically profile samples using Meteor, but at a downsized level
# we use the same mapping produced by Meteor in previous steps
rule meteor_profiling_downsized:
    input:
        "results/09_taxonomic_profiling/meteor/{sample}/mapping"
    output:
        directory("results/09_taxonomic_profiling/meteor/{sample}/profiling_downsized_{downsize}")
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/meteor/{sample}_profiling_downsized_{downsize}M.stdout",
        stderr = "logs/09_taxonomic_profiling/meteor/{sample}_profiling_downsized_{downsize}M.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/meteor/{sample}.profile_downsized_{downsize}M.benchmark.txt"
    params:
        mapping_with_sample = lambda wildcards: os.path.join(f"results/09_taxonomic_profiling/meteor/{wildcards.sample}/mapping", f"{wildcards.sample}"),  
        downsize_int = lambda wildcards: convert_from_si_units_to_int(wildcards.downsize),
        reference = METEOR_REFERENCE
    wildcard_constraints:
        downsize = "|".join([convert_to_si_units(int(size)) for size in config.get('taxonomic_profiling', {}).get('meteor', {}).get('downsize', {})])
    shell:
        """
        meteor profile -i {params.mapping_with_sample} -o {output} -r {params.reference} \
            -n coverage \
            --seed 100 \
            -l {params.downsize_int} \
            > {log.stdout} 2> {log.stderr}
        """

rule strainscan_profiling:
    input: 
        r1 = "results/02_preprocess/bowtie2/{sample}_1.clean.fastq.gz", # on short reads only
        r2 = "results/02_preprocess/bowtie2/{sample}_2.clean.fastq.gz" # on short reads only
    output: "results/09_taxonomic_profiling/strainscan/{sample}/final_report.txt"
    conda:
        "../envs/strainscan.yaml"
    log:
        stdout = "logs/09_taxonomic_profiling/strainscan/{sample}.profile.stdout",
        stderr = "logs/09_taxonomic_profiling/strainscan/{sample}.profile.stderr"
    benchmark:
        "benchmarks/09_taxonomic_profiling/strainscan/{sample}.profile.benchmark.txt"
    params:
        output_dir = lambda wildcards: f"results/09_taxonomic_profiling/strainscan/{wildcards.sample}",
        reference = STRAINSCAN_DB
    wildcard_constraints:
        sample = "|".join(SAMPLES)
    shell:
        """
        strainscan \
            -i {input.r1} \
            -j {input.r2} \
            -d {params.reference} \
            -o {params.output_dir} \
        > {log.stdout} 2> {log.stderr}
        """
