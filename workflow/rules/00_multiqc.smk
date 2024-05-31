from rules.utils import *

SAMPLES_TABLE = config['samples']
SAMPLES = read_table(SAMPLES_TABLE)
SAMPLES_DF = pd.read_csv(SAMPLES_TABLE, sep="\t")
FASTQ_FILES = SAMPLES_DF['sample'].tolist()
FASTQ_FILES = [f[:-9] for f in FASTQ_FILES]

READS = [1, 2]
ASSEMBLER = config['assembly']['assembler']

# produces a MultiQC report from supported tools' result
rule multiqc:
  input: 
    expand("results/01_qc/fastqc/{fastq_files}/", fastq_files=FASTQ_FILES),
    expand("results/02_preprocess/fastqc/{sample}_{read}.clean_fastqc.html", 
    sample=SAMPLES, read=READS),
    expand("results/03_assembly/{assembler}/{sample}/assembly.fa.gz", 
      assembler=ASSEMBLER, sample=SAMPLES),
    expand("results/04_assembly_qc/quast/{assembler}/{sample}/report.html", 
      assembler=ASSEMBLER, sample=SAMPLES),
    expand("results/06_binning_qc/checkm2/samples/{sample}/all_quality_reports.pdf",
      sample=SAMPLES),
    expand("results/07_bins_refinement/binette/{assembler}/{sample}", 
      assembler=ASSEMBLER, sample=SAMPLES),
    expand("results/08_bins_postprocessing/gtdb_tk/{assembler}/{sample}", 
      assembler=ASSEMBLER, sample=SAMPLES),
    expand("results/08_bins_postprocessing/dRep/{assembler}",
      assembler=ASSEMBLER),
  output: directory("results/00_multiqc")
  conda: "../envs/multiqc.yaml"
  log:
    stdout = "logs/00_multiqc/mutiqc.stdout",
    stderr = "logs/00_multiqc/mutiqc.stderr"
  shell:
    """
    multiqc --outdir {output} results/
    """