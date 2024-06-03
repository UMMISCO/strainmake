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
