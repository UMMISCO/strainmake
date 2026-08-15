################################################################################
#                             Reference databases                              #
################################################################################
#
# Every rule that reaches out to the internet lives here.
#
# That grouping intends to make StrainMake usable on clusters whose compute nodes
# have no route to the internet: the whole set can be materialised in one go
# from a login node with
#
#     strainmake fetch-databases --config <config.yaml>
#
# and the pipeline proper then runs offline against the resulting
# directory. Database locations come from the `databases:` config section (see
# `database_path` in rules/utils.py), so a single downloaded copy can be shared
# between projects and users instead of being re-fetched per run.
#
# The rules below are only pulled into the DAG when their output is missing, so
# an "ordinary" online run still fetches what it needs on demand, exactly as
# before. When `deployment: offline:` is set they abort immediately with a
# message naming the database and the command that would provide it.

import os

from utils import *

# resolved once, at parse time, so every consuming rule agrees on the locations
BOWTIE2_INDEX_DB = database_path(config, "bowtie2_index")
CHECKM2_DB = database_path(config, "checkm2")
BAKTA_DB = database_path(config, "bakta")
METAPHLAN_DB = database_path(config, "metaphlan")
GTDBTK_DB = database_path(config, "gtdbtk")
METEOR_DB = database_path(config, "meteor")
STRAINSCAN_DB = database_path(config, "strainscan")

# the catalogue subdirectory `meteor download` produces, and what `meteor -r` wants
METEOR_REFERENCE = meteor_reference(config)

# the `db` subdirectory `bakta --db` expects inside the download
BAKTA_DB_DIR = bakta_db_dir(config)

# what each download rule produces, and what every consuming rule depends on.
# These are marker files rather than directories: a `directory()` output would
# let Snakemake delete a multi-hundred-gigabyte download whenever the rule is
# re-triggered. `database_target` is shared with the `strainmake` CLI, so that
# `strainmake fetch-databases` asks Snakemake for exactly these paths.
CHECKM2_DB_MARKER = database_target(config, "checkm2")
BAKTA_DB_MARKER = database_target(config, "bakta")
METAPHLAN_DB_MARKER = database_target(config, "metaphlan")
GTDBTK_DB_MARKER = database_target(config, "gtdbtk")
METEOR_DB_MARKER = database_target(config, "meteor")


# the bowtie2 index of the host genome, used to decontaminate the metagenome
rule get_bowtie_index:
    output:
        directory(BOWTIE2_INDEX_DB)
    log:
        wget_stdout = "logs/00_databases/bowtie2/get_bowtie_index.wget.stdout",
        wget_stderr = "logs/00_databases/bowtie2/get_bowtie_index.wget.stderr",
        unzip_stdout = "logs/00_databases/bowtie2/get_bowtie_index.unzip.stdout",
        unzip_stderr = "logs/00_databases/bowtie2/get_bowtie_index.unzip.stderr",
        mv_stdout = "logs/00_databases/bowtie2/get_bowtie_index.mv.stdout",
        mv_stderr = "logs/00_databases/bowtie2/get_bowtie_index.mv.stderr",
        rm_stdout = "logs/00_databases/bowtie2/get_bowtie_index.rm.stdout",
        rm_stderr = "logs/00_databases/bowtie2/get_bowtie_index.rm.stderr"
    conda:
        "../envs/get_bowtie_index.yaml"
    params:
        organism_name = config['bowtie2']['index_name'],
        offline_guard = offline_guard(config, "bowtie2 index", "bowtie2_index"),
        archive = os.path.join(BOWTIE2_INDEX_DB, "..", config['bowtie2']['index_name'] + ".zip")
    benchmark:
        "benchmarks/00_databases/bowtie2/get_bowtie_index.benchmark.txt"
    shell:
    # downloading the already made index file and unzipping it. The indexes
    # will be in a folder of name "index"
        """
        {params.offline_guard} \
        mkdir -p {output} \
        && wget -O {params.archive} \
            https://genome-idx.s3.amazonaws.com/bt/{params.organism_name}.zip \
            > {log.wget_stdout} 2> {log.wget_stderr} \
        && unzip -d {output} {params.archive} \
            > {log.unzip_stdout} 2> {log.unzip_stderr} \
        && mv {output}/{params.organism_name}/* {output}/ \
            > {log.mv_stdout} 2> {log.mv_stderr} \
        && rm -r {output}/{params.organism_name} {params.archive} \
            > {log.rm_stdout} 2> {log.rm_stderr}
        """


# the DIAMOND database CheckM2 scores bins against
rule checkm2_database:
    output:
        CHECKM2_DB_MARKER
    conda:
        "../envs/checkm2.yaml"
    log:
        stdout = "logs/00_databases/checkm2/checkm2.db.stdout",
        stderr = "logs/00_databases/checkm2/checkm2.db.stderr"
    benchmark:
        "benchmarks/00_databases/checkm2/checkm2.db.benchmark.txt"
    params:
        output_path = CHECKM2_DB,
        offline_guard = offline_guard(config, "CheckM2", "checkm2")
    shell:
        """
        {params.offline_guard} \
        checkm2 database --download --path {params.output_path} \
            > {log.stdout} 2> {log.stderr}
        """


# the Bakta database used to annotate bacterial MAGs
rule bakta_get_database:
    output:
        BAKTA_DB_MARKER
    conda:
        "../envs/bakta.yaml"
    log:
        stdout = "logs/00_databases/bakta/database.stdout",
        stderr = "logs/00_databases/bakta/database.stderr"
    benchmark:
        "benchmarks/00_databases/bakta/database.benchmark.txt"
    params:
        output_path = BAKTA_DB,
        offline_guard = offline_guard(config, "Bakta", "bakta")
    shell:
        """
        {params.offline_guard} \
        bakta_db download --type full --output {params.output_path} \
            > {log.stdout} 2> {log.stderr}
        """


# the MetaPhlAn marker gene database. Without an explicit --db_dir MetaPhlAn
# downloads this into its own conda environment, which fails offline and, with
# a shared conda prefix, writes into environments other projects also use
rule metaphlan_database:
    output:
        METAPHLAN_DB_MARKER
    conda:
        "../envs/metaphlan.yaml"
    log:
        stdout = "logs/00_databases/metaphlan/database.stdout",
        stderr = "logs/00_databases/metaphlan/database.stderr"
    benchmark:
        "benchmarks/00_databases/metaphlan/database.benchmark.txt"
    params:
        output_path = METAPHLAN_DB,
        offline_guard = offline_guard(config, "MetaPhlAn", "metaphlan")
    threads: config.get('taxonomic_profiling', {}).get('metaphlan', {}).get('threads', 1)
    shell:
        """
        {params.offline_guard} \
        mkdir -p {params.output_path} \
        && metaphlan --install --nproc {threads} --db_dir {params.output_path} \
            > {log.stdout} 2> {log.stderr}
        """


# the GTDB reference data used to assign taxonomy to MAGs.
# this is by far the largest download of the pipeline (~101 GiB for r220)
rule gtdbtk_database:
    output:
        GTDBTK_DB_MARKER
    log:
        stdout = "logs/00_databases/gtdb_tk/database.stdout",
        stderr = "logs/00_databases/gtdb_tk/database.stderr"
    benchmark:
        "benchmarks/00_databases/gtdb_tk/database.benchmark.txt"
    conda:
        "../envs/get_bowtie_index.yaml"  # provides wget and tar
    params:
        output_path = GTDBTK_DB,
        release_url = database_setting(
            config,
            "gtdbtk_release_url",
            "https://data.gtdb.ecogenomic.org/releases/release220/220.0/"
            "auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz",
        ),
        archive = os.path.join(GTDBTK_DB, "gtdbtk_data.tar.gz"),
        offline_guard = offline_guard(config, "GTDB-Tk", "gtdbtk")
    shell:
    # --continue lets an interrupted transfer of this size be resumed rather
    # than restarted, and the tarball is removed once unpacked
        """
        {params.offline_guard} \
        mkdir -p {params.output_path} \
        && wget --continue -O {params.archive} {params.release_url} \
            > {log.stdout} 2> {log.stderr} \
        && tar xzf {params.archive} -C {params.output_path} --strip-components=1 \
            >> {log.stdout} 2>> {log.stderr} \
        && rm -f {params.archive}
        """


# the Meteor gene catalogue used for taxonomic profiling
rule meteor_database:
    output:
        directory(METEOR_DB_MARKER)
    conda:
        "../envs/meteor.yaml"
    log:
        stdout = "logs/00_databases/meteor/database.stdout",
        stderr = "logs/00_databases/meteor/database.stderr"
    benchmark:
        "benchmarks/00_databases/meteor/database.benchmark.txt"
    params:
        output_path = METEOR_DB,
        catalogue = database_setting(config, "meteor_catalogue", "hs_10_4_gut"),
        offline_guard = offline_guard(config, "Meteor", "meteor")
    shell:
        """
        {params.offline_guard} \
        mkdir -p {params.output_path} \
        && meteor download -i {params.catalogue} -o {params.output_path} -c \
            > {log.stdout} 2> {log.stderr}
        """


# NOTE: there is deliberately no `strainscan_database` rule. StrainScan's
# prebuilt databases are only published as per-species Google Drive links
# (https://github.com/liaoherui/StrainScan#pre-built-databases-download), which
# cannot be fetched unattended. `strainmake fetch-databases` reports it as a
# manual step and the preflight check verifies the configured path exists.
