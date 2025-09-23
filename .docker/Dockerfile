FROM continuumio/miniconda3:25.3.1-1

# metadata
LABEL org.opencontainers.image.url="https://github.com/UMMISCO/strainmake"

ARG PIPELINE_COMMIT_SHA="f533a52ed9b"
ARG PIPELINE_REPO="https://github.com/UMMISCO/strainmake.git"
ARG SNAKEMAKE_VERSION="8.24.1"
ARG CONDA_ENV_NAME="snakemake${SNAKEMAKE_VERSION}"

# installing Snakemake
RUN conda create -n ${CONDA_ENV_NAME} -c bioconda -c conda-forge -y snakemake=${SNAKEMAKE_VERSION}
# making Snakemake conda env the default one
RUN echo "source activate ${CONDA_ENV_NAME}" > ~/.bashrc

# getting pipeline code
RUN git clone ${PIPELINE_REPO} /opt/strainmake
WORKDIR /opt/strainmake
# checkout the desired commit
RUN git checkout ${PIPELINE_COMMIT_SHA}
# where the results, logs, benchmarks, will be stored
# the user should mount a volume to /opt/strainmake when running the container in order to keep the data
RUN mkdir -p /opt/strainmake/results /opt/strainmake/logs /opt/strainmake/benchmarks

ENTRYPOINT ["/opt/conda/envs/snakemake8.24.1/bin/snakemake"]