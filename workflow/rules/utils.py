import os

import pandas as pd

################################################################################
#                             Reference databases                              #
################################################################################

# Where each reference database lives when the configuration does not say.
DATABASE_DEFAULT_PATHS = {
    "bowtie2_index": "data/databases/bowtie2",
    "checkm2": "data/databases/checkm2",
    "bakta": "data/databases/bakta",
    "metaphlan": "data/databases/metaphlan",
    "gtdbtk": "data/databases/gtdb_tk",
    "meteor": "data/databases/meteor",
    "strainscan": "data/databases/strainscan",
}

# human-readable name of each database, used in error messages and reports
DATABASE_LABELS = {
    "bowtie2_index": "bowtie2 host index",
    "checkm2": "CheckM2",
    "bakta": "Bakta",
    "metaphlan": "MetaPhlAn",
    "gtdbtk": "GTDB-Tk",
    "meteor": "Meteor",
    "strainscan": "StrainScan",
}

# databases StrainMake knows how to download itself, and the tool that does it.
# StrainScan is deliberately absent: its prebuilt databases are only distributed
# as per-species Google Drive links, which cannot be fetched unattended
DOWNLOADABLE_DATABASES = {
    "bowtie2_index": "wget",
    "checkm2": "checkm2 database --download",
    "bakta": "bakta_db download",
    "metaphlan": "metaphlan --install",
    "gtdbtk": "GTDB release tarball",
    "meteor": "meteor download",
}

# rough download size of each database, so that a run that is about to fetch one
# can say so before it starts rather than after an unexplained hour
DATABASE_APPROX_SIZES = {
    "bowtie2_index": "~4 GiB",
    "checkm2": "~3 GiB",
    "bakta": "~70 GiB",
    "metaphlan": "~6 GiB",
    "gtdbtk": "~105 GiB",
    "meteor": "~10 GiB",
}

# where to obtain the databases StrainMake cannot fetch on its own
MANUAL_DATABASE_SOURCES = {
    "strainscan": (
        "https://github.com/liaoherui/StrainScan"
        "#pre-built-databases-download"
    ),
}

# config section whose presence means the pipeline will use a given database.
# `None` marks a database the pipeline always needs.
DATABASE_REQUIRED_BY = {
    "bowtie2_index": None,  # host decontamination feeds nearly every other step
    "checkm2": ("binning",),
    "bakta": ("bins_postprocessing", "bakta"),
    "gtdbtk": ("bins_postprocessing", "gtdbtk"),
    "metaphlan": ("taxonomic_profiling", "metaphlan"),
    "meteor": ("taxonomic_profiling", "meteor"),
    "strainscan": ("taxonomic_profiling", "strainscan"),
}


def database_path(config: dict, key: str) -> str:
    """
    Resolve where a reference database lives on disk.

    Each database is located independently, from its own `databases: <key>:`
    entry, because in practice they are not always kept together

    Parameters:
    config (dict): the pipeline configuration
    key (str): a key of DATABASE_DEFAULT_PATHS, e.g. 'checkm2'

    Returns:
    str: path to the database directory
    """

    if key not in DATABASE_DEFAULT_PATHS:
        raise ValueError(
            f"Unknown database '{key}'. "
            f"Known databases: {sorted(DATABASE_DEFAULT_PATHS)}"
        )

    databases = config.get("databases") or {}

    configured = databases.get(key)

    # `false` marks a database as unused, not as living at a path named "False"
    if configured and configured is not True:
        return str(configured)

    return DATABASE_DEFAULT_PATHS[key]


def database_setting(config: dict, key: str, default):
    """
    Read a non-path entry of the `databases:` config section.

    Tolerates configs written before the section existed, so an old config file
    keeps working instead of failing with a KeyError at parse time.

    Parameters:
    config (dict): the pipeline configuration
    key (str): entry to read, e.g. 'meteor_catalogue'
    default: value to fall back on when the entry is absent or empty

    Returns:
    the configured value, or `default`
    """

    return (config.get("databases") or {}).get(key) or default


def bakta_db_dir(config: dict) -> str:
    """
    Resolve the directory `bakta --db` expects.

    `bakta_db download --type full` unpacks the database into a `db`
    subdirectory of the download location.

    Returns:
    str: path to the Bakta database directory
    """

    return os.path.join(database_path(config, "bakta"), "db")


def database_target(config: dict, key: str) -> str:
    """
    The path Snakemake builds in order to materialise a database.

    For most databases this is a file inside the download whose presence proves
    the transfer completed. Pointing a rule at the enclosing directory instead
    would let Snakemake delete a multi-hundred-gigabyte download whenever the
    rule is re-triggered.

    Parameters:
    config (dict): the pipeline configuration
    key (str): a key of DATABASE_DEFAULT_PATHS, e.g. 'checkm2'

    Returns:
    str: the path to request from Snakemake to obtain this database
    """

    root = database_path(config, key)

    markers = {
        # the bowtie2 index has no single canonical file name: its contents
        # depend on which index the user selected, so the directory is the target
        "bowtie2_index": root,
        "checkm2": os.path.join(root, "CheckM2_database", "uniref100.KO.1.dmnd"),
        "bakta": os.path.join(bakta_db_dir(config), "version.json"),
        "metaphlan": os.path.join(root, "mpa_latest"),
        "gtdbtk": os.path.join(root, "metadata", "metadata.txt"),
        "meteor": meteor_reference(config),
        # StrainScan is provided by the user, never built by the pipeline
        "strainscan": root,
    }

    if key not in markers:
        raise ValueError(
            f"Unknown database '{key}'. Known databases: {sorted(markers)}"
        )

    return markers[key]


def database_is_required(config: dict, key: str) -> bool:
    """
    Whether the pipeline, as configured, will actually need a database.

    The rule is deliberately simple, so that it is predictable from reading the
    config alone: a database is required when the config section naming its tool
    is *present*, whatever its contents. Deleting the section is how a user
    declares they do not use that tool, which stops `strainmake fetch-databases`
    from downloading it and stops the preflight check from asking for it.

    Parameters:
    config (dict): the pipeline configuration
    key (str): a key of DATABASE_DEFAULT_PATHS, e.g. 'checkm2'

    Returns:
    bool: True when a configured step uses this database
    """

    # an explicit `false` marks a database this run does not use. Needed for
    # databases whose tool section cannot simply be deleted: `bowtie2:` holds
    # settings the preprocessing rules read at parse time, so it stays in the
    # config even for a run that starts from already-preprocessed reads and
    # never touches the host index.
    databases = config.get("databases") or {}
    if databases.get(key, None) is False:
        return False

    section = DATABASE_REQUIRED_BY.get(key, None)

    # databases with no gating section are always needed
    if section is None:
        return True

    node = config
    for level in section:
        if not isinstance(node, dict) or level not in node:
            return False
        node = node[level]

    return True


def meteor_reference(config: dict) -> str:
    """
    Resolve the Meteor catalogue directory passed to `meteor -r`.

    `meteor download` unpacks the catalogue into a subdirectory of the download
    location named after the catalogue itself, and that subdirectory is what
    Meteor expects as its reference.

    Returns:
    str: path to the catalogue directory
    """

    catalogue = database_setting(config, "meteor_catalogue", "hs_10_4_gut")

    return os.path.join(database_path(config, "meteor"), str(catalogue))


def is_offline(config: dict) -> bool:
    """
    Whether the pipeline was told it cannot reach the internet.

    Returns:
    bool: True when `deployment: offline:` is set
    """

    return bool((config.get("deployment") or {}).get("offline", False))


def offline_guard(config: dict, database_label: str, key: str) -> str:
    """
    Build a shell prelude that aborts a download rule when running offline.

    Without it, a compute node with no route to the internet spends minutes in
    urllib3/wget retries before failing with a stack trace that says nothing
    about what the user should do. With it, the rule stops immediately and says
    which command to run on a machine that does have network access.

    Parameters:
    config (dict): the pipeline configuration
    database_label (str): human-readable database name, e.g. 'CheckM2'
    key (str): the `databases:` config key holding this database's path

    Returns:
    str: shell code to prefix the download command with (empty when online)
    """

    if not is_offline(config):
        return ""

    path = database_path(config, key)

    message = (
        f"StrainMake is running in offline mode and the {database_label} database "
        f"is missing from '{path}'."
    )
    hint = (
        "Run 'strainmake fetch-databases' on a machine with internet access "
        f"or point 'databases: {key}:' at an existing copy."
    )

    # `>&2` so the message lands in the rule's stderr log next to the failure
    return f'echo "ERROR: {message}" >&2; echo "{hint}" >&2; exit 1; '


def gurobi_license_path(config: dict) -> str:
    """
    Resolve the Gurobi license file used by CarveMe.

    Returns:
    str: path to the license file (may not exist; CarveMe falls back to SCIP)
    """

    carveme = ((config.get("bins_postprocessing") or {}).get("carveme") or {})

    return str(carveme.get("gurobi_license") or "config/gurobi.lic")


# get samples name
def read_table(table_path: str):
       """
       This functions returns the list of sample to use in the experiment
       (the ones with small-read experiments)
       """

       df = pd.read_csv(table_path, sep="\t")
       filtered_df = df[df['type'] == 'R1']

       sample_id_list = filtered_df['sample_id'].unique().tolist()

       return sample_id_list

def read_table_long_reads(table_path: str):
       """
       This functions returns the list of sample to use in the experiment that
       have a long-read experiment
       """

       df = pd.read_csv(table_path, sep="\t")
       filtered_df = df[df['type'] == 'long']

       sample_id_list = filtered_df['sample_id'].unique().tolist()

       return sample_id_list

def get_fastq_pair(df, sample_id):
       """
       This function returns the pair of FASTQ corresponding to a sample (small reads)
       """
       sample_df = df[df['sample_id'] == sample_id]

       forward_read = sample_df[sample_df['type'] == 'R1']['sample'].values[0]
       reverse_read = sample_df[sample_df['type'] == 'R2']['sample'].values[0]

       return (forward_read, reverse_read)

def get_fastq_long_read(df, sample_id):
       """
       This function returns the long read FASTQ corresponding to a sample
       """
       sample_df = df[df['sample_id'] == sample_id]

       long_read = sample_df[sample_df['type'] == 'long']['sample'].values[0]

       return long_read

def get_all_fastq(df):
       """
       This function returns all the FASTQ
       """

       sample = df['sample'].tolist()

       return set(sample)

def validate_assemblers(df: pd.DataFrame, assemblers: list):
    """
    This function validates the availability of input files required for the specified
    assemblers
    
    Parameters:
    df (pd.DataFrame): Metadata TSV used by the pipeline, read into a DataFrame.
                       It should have three columns: 'sample_id', 'type', and 'sample'
    assemblers (list): A list of assemblers to be used. Valid options are 'metaflye', 
                       'metaspades', 'megahit', and 'hybridspades'.
    
    Raises:
    ValueError: If the required input files are missing for any of the specified assemblers
    """

    # check if metaflye is selected and there are no long reads
    if 'metaflye' in assemblers:
        if not any(df['type'] == 'long'):
            raise ValueError("Metaflye assembler requires at least one 'long' read file.")
    
    # check if hybridspades is selected and there is no combination of R1, R2, and long reads for any sample
    if 'hybridspades' in assemblers:
        hybrid_valid = df.groupby('sample_id').apply(
            lambda x: all(type_ in x['type'].values for type_ in ['R1', 'R2', 'long'])
        ).any()
        
        if not hybrid_valid:
            raise ValueError("Hybridspaes assembler requires 'R1', 'R2', and 'long' read files for at least one sample.")
    
    # check if megahit or metaspades are selected and there are no R1 and R2 reads for any sample
    if any(asm in assemblers for asm in ['megahit', 'metaspades']):
        sr_valid = df.groupby('sample_id').apply(
            lambda x: all(type_ in x['type'].values for type_ in ['R1', 'R2'])
        ).any()
        
        if not sr_valid:
            raise ValueError("Megahit and Metaspades assemblers require 'R1' and 'R2' read files for at least one sample.")

def convert_to_si_units(value):
    """
    Convert a numerical value to a string representation using SI units

    Parameters:
    value (float): The numerical value to convert

    Returns:
    str: The value converted to a string with an appropriate SI unit suffix
    """

    units = ["", "K", "M", "G", "T", "P", "E", "Z", "Y"]
    power = 1000
    n = 0

    while value >= power and n < len(units) - 1:
        value /= power
        n += 1

    return f"{int(value)}{units[n]}"

def convert_from_si_units_to_int(value):
    """
    Convert a string representation of a numerical value using SI units to an integer

    Parameters:
    value (str): The string representation of the value with an SI unit suffix

    Returns:
    int: The value converted to an integer
    """

    units = {"": 1, "K": 10**3, "M": 10**6, "G": 10**9, "T": 10**12, "P": 10**15, "E": 10**18, "Z": 10**21, "Y": 10**24}
    unit = value[-1]

    if unit.isdigit() or unit == '.':
        return float(value)

    num = float(value[:-1])

    return int(num * units[unit])
