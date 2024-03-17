FROM condaforge/mambaforge:latest
LABEL version="0.1.0"

# installing seqkit
RUN mamba install -y bioconda::seqkit=2.8.*

# installing snakemake
RUN mamba install -y bioconda::snakemake=7.32.*
