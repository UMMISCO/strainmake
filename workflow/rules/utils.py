import pandas as pd

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

