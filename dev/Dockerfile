FROM condaforge/mambaforge:latest
LABEL version="0.2.0"

# installing seqkit
RUN mamba install -y bioconda::seqkit=2.8.*

# installing snakemake
RUN mamba install -y bioconda::snakemake=7.32.*

# installing screen
RUN mamba install -y conda-forge::screen=4.8.*