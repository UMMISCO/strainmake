FROM condaforge/mambaforge:latest
LABEL version="0.3.0"

# installing seqkit
RUN mamba install -y bioconda::seqkit=2.8.*

# installing snakemake
RUN mamba install -y bioconda::snakemake=7.32.*

# installing screen
RUN mamba install -y conda-forge::screen=4.8.*

# adding an alias for always timing Snakemake runtime
RUN echo "alias snakemake='time snakemake'" >> ~/.bashrc