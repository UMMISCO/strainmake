"""
Script for merging CheckM2 quality reports and adding a new column with binning
program name and a new column with assembly program, extracted from the 
table path
"""

import argparse
import pandas as pd

def identify_assembly_program(report_path: str):
    """
    Returns the assembly program it identified from the report path
    """

    print(f'Extracting assembly program name from: {report_path}')

    if 'metaspades' in report_path:
        return 'metaspades'
    elif 'megahit' in report_path:
        return 'megahit'
    elif 'metaflye' in report_path:
        return 'metaflye'
    elif 'hybridspades' in report_path:
        return 'hybridspades'
    else:
        return None
    
def identify_binning_program(report_path: str):
    """
    Returns the binning program it identified from the report path
    """

    print(f'Extracting binning program name from: {report_path}')

    if 'semibin2' in report_path:
        return 'sembin2'
    elif 'metabat2' in report_path:
        return 'metabat2'
    elif 'vamb' in report_path:
        return 'vamb'
    else:
        return None
    
def merge_reports(input_files: list):
    """
    Returns a dataframe being the CheckM2 merged reports
    """

    # loading the tables
    df_list = [pd.read_csv(file, sep="\t") for file in input_files]

    # identifying the binning program
    df_with_binning_program = []
    for df, path_to_report in zip(df_list, input_files):
        df['assembly'] = identify_assembly_program(path_to_report)
        df['binning'] = identify_binning_program(path_to_report)

        print(f'Will be merged:\n\n{df}')

        df_with_binning_program.append(df)

    # merging the reports
    all_reports_merged = pd.concat(df_with_binning_program, axis=0,
                                   ignore_index=True) 

    return all_reports_merged
    
def main():
    parser = argparse.ArgumentParser(description="Merge multiple CheckM2 reports into one.")
    parser.add_argument("--input", nargs="+", 
                        help="Input TSV report to merge.", required=True)
    parser.add_argument("--output", help="Output merged TSV report.", 
                        required=True)
    
    args = parser.parse_args()
    
    merged_reports = merge_reports(args.input)

    # saving the merged reports
    merged_reports.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()