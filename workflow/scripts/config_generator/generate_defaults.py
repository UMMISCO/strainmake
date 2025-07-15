""" 
This script reads a YAML configuration template, extracts specific sub-dictionaries corresponding to 
different pipeline sections (preprocessing, assembly, gene catalog, binning, and taxonomic profiling), 
and writes them as Python dictionaries to a 'defaults.py' file for import. The full defaults 
dictionary is also exported.

Update this script when the YAML template changes.
"""

import yaml
from pprint import pformat

with open('../../../config/template_config.yaml', 'r') as f:
    defaults = yaml.safe_load(f)

# extracting sub-dictionaries for export
PREPROCESSING = {k: v for k, v in defaults.items() if k in ['fastp', 'fastp_long_read', 'bowtie2', 'downsizing_for_hybrid']}
ASSEMBLY = {k: v for k, v in defaults.items() if k in ['assembly', 'quast']}
GENE_CATALOG = {k: v for k, v in defaults.items() if k in ['mmseqs2', 'representative_genes']}
BINNING = {k: v for k, v in defaults.items() if k in ['binning', 'checkm2', 'bins_refinement', 'bins_postprocessing']}
TAXO_PROFILING = {k: v for k, v in defaults.items() if k in ['taxonomic_profiling']}

# writing the extracted configuration sections and the full defaults dictionary
# to a Python file ('defaults.py') for easy import and reuse in other scripts.
with open('defaults.py', 'w') as out:
    for section in ['PREPROCESSING', 'ASSEMBLY', 'GENE_CATALOG', 'BINNING', 'TAXO_PROFILING']:
        out.write(f"{section} = ")
        out.write(pformat(locals()[section], sort_dicts=False))
        out.write('\n')
        out.write('\n')

    # exporting the full defaults as well if needed
    out.write('defaults = ')
    out.write(pformat(defaults, sort_dicts=False))
    out.write('\n')