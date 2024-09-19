configfile: "config/template_config.yaml"

include: "rules/02_prepocess.smk"
include: "rules/01_qc.smk"
include: "rules/02_preprocess_LR.smk"
include: "rules/03_assembly.smk"
include: "rules/03_assembly_LR.smk"
include: "rules/04_assembly_qc.smk"
include: "rules/04_assembly_qc_LR.smk"
include: "rules/05_binning.smk"
include: "rules/05_binning_LR.smk"
include: "rules/06_binning_qc.smk"
include: "rules/06_binning_qc_LR.smk"
include: "rules/07_bins_refinement.smk"
include: "rules/07_bins_refinement_LR.smk"
include: "rules/08_bins_postprocessing.smk"
include: "rules/09_taxonomic_profiling.smk"
include: "rules/10_strains_profiling.smk"

from rules.utils import *

################################################################################
#                                Config options                                #
################################################################################

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
SAMPLES_LR = read_table_long_reads(SAMPLES_TABLE)
SAMPLES_DF = pd.read_csv(SAMPLES_TABLE, sep="\t")
FASTQ_FILES = SAMPLES_DF['sample'].tolist()
FASTQ_FILES = [f[:-9] for f in FASTQ_FILES]

READS = [1, 2]
ASSEMBLER = config['assembly']['assembler']
LONG_READ_ASSEMBLER = config['assembly']['long_read_assembler']
HYBRID_ASSEMBLER = config['assembly']['hybrid_assembler']

# taking into account the case where we don't have SR
if ASSEMBLER == None:
       ASSEMBLER = []

# taking into account the case where we don't have LR
if LONG_READ_ASSEMBLER == None:
       LONG_READ_ASSEMBLER = []

# taking into account the case where we don't have hybrid assembler
if HYBRID_ASSEMBLER == None:
       HYBRID_ASSEMBLER = []

# validating that we can run the pipeline using the given FASTQ
samples_table_df = pd.read_csv(SAMPLES_TABLE, sep="\t")
validate_assemblers(samples_table_df, 
                    ASSEMBLER + LONG_READ_ASSEMBLER + HYBRID_ASSEMBLER)

################################################################################
#                                   Workflow                                   #
################################################################################

rule all:
        input:
            # assembly part
            expand("results/03_assembly/{assembler}/{sample}/assembly.fa.gz", 
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES),
            expand("results/03_assembly/LR/{assembler_long_read}/{sample_lr}/assembly.fa.gz", 
                   assembler_long_read=LONG_READ_ASSEMBLER, sample_lr=SAMPLES_LR),
            # assembly qc
            expand("results/04_assembly_qc/quast/{assembler}/{sample}/report.html", 
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES),
            expand("results/04_assembly_qc/quast_long_read/{assembler_long_read}/{sample_lr}/report.html", 
                   assembler_long_read=LONG_READ_ASSEMBLER, sample_lr=SAMPLES_LR),
            # binning (step 05) + bins qc (step 06)
            expand("results/06_binning_qc/checkm2/samples/{sample}/all_quality_reports.pdf",
                   sample=SAMPLES),
            # bins refinement
            expand("results/07_bins_refinement/binette/{assembler}/{sample}", 
                   assembler=ASSEMBLER + LONG_READ_ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES + SAMPLES_LR),
            # bins post-processing
            expand("results/08_bins_postprocessing/dRep/{assembler}",
                   assembler=ASSEMBLER + LONG_READ_ASSEMBLER + HYBRID_ASSEMBLER),
            expand("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{assembler}/bins",
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER + LONG_READ_ASSEMBLER),
                               expand("results/08_bins_postprocessing/gtdb_tk/{assembler}", 
                   assembler=ASSEMBLER + LONG_READ_ASSEMBLER + HYBRID_ASSEMBLER),
            expand("results/08_bins_postprocessing/checkm1/{assembler}/{sample}/profile.processed.tsv",
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES),
            # taxonomic profiling
            expand("results/09_taxonomic_profiling/metaphlan/{sample}.profile.txt",
                   sample=SAMPLES),
            # strains profiling
            #      -> inStrain
            expand("results/10_strain_profiling/inStrain/{assembler}/compare",
                   assembler = ASSEMBLER + HYBRID_ASSEMBLER),
            #      -> Floria
            expand("results/10_strain_profiling/floria/{assembler}/{sample}/contig_ploidy_info.tsv",
                   assembler = ASSEMBLER + HYBRID_ASSEMBLER, sample = SAMPLES)
