Ongoing implementation of several metagenomic tools into a single pipeline.

# Roadmap

## Quality control, preprocessing

| Tool   | Original year                                         | Conda available?                            | Link                                                      | Implemented? |
| :----- | :---------------------------------------------------- | :------------------------------------------ | :-------------------------------------------------------- | :----------- |
| fastp  | [2018](https://doi.org/10.1093/bioinformatics/bty560) | [Yes](https://anaconda.org/bioconda/fastp)  | https://github.com/OpenGene/fastp                         | Yes          |
| fastQC | 2010                                                  | [Yes](https://anaconda.org/bioconda/fastqc) | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ | Yes          |

### Human decontamination

| Tool    | Original year                              | Conda available?                             | Link                                   | Implemented? |
| :------ | :----------------------------------------- | :------------------------------------------- | :------------------------------------- | :----------- |
| bowtie2 | [2012](https://doi.org/10.1038/nmeth.1923) | [Yes](https://anaconda.org/bioconda/bowtie2) | https://github.com/BenLangmead/bowtie2 | Yes          |

Human assembly for mapping: 

* GRCh38. [Link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)

## Asssembly

| Tool         | Original year                                                  | Conda available?                                | Link                                                          | Implemented? |
| :----------- | :------------------------------------------------------------- | :---------------------------------------------- | :------------------------------------------------------------ | :----------- |
| MEGAHIT      | [2015](https://doi.org/10.1093/bioinformatics/btv033)          | [Yes](https://anaconda.org/bioconda/megahit)    | https://github.com/voutcn/megahit                             | Yes          |
| (Meta)SPAdes | [2017](https://doi.org/10.1101/gr.213959.116)                  | [Yes](https://anaconda.org/bioconda/spades)     | https://github.com/ablab/spades                               | Yes          |
| IDBA-UD      | [2012](https://doi.org/10.1093/bioinformatics/bts174)          | [Yes](https://anaconda.org/bioconda/idba)       | https://github.com/loneknightpy/idba                          | No           |
| MetaVelvet   | [2012](https://doi.org/10.1093/nar/gks678)                     | [Yes](https://anaconda.org/bioconda/metavelvet) | https://github.com/hacchy/MetaVelvet                          | No           |
| MetaCompass  | [2024](https://arxiv.org/ftp/arxiv/papers/2403/2403.01578.pdf) | No                                              | https://gitlab.umiacs.umd.edu/mpop/metacompass (can't access) | No           |

## Assembly quality assessment

| Tool  | Original year                                         | Conda available?                           | Link                           | Implemented? |
| :---- | :---------------------------------------------------- | :----------------------------------------- | :----------------------------- | :----------- |
| QUAST | [2013](https://doi.org/10.1093/bioinformatics/btt086) | [Yes](https://anaconda.org/bioconda/quast) | https://github.com/ablab/quast | No           |

## Binning

## Gene prediction, functional annotation and taxonomy classification