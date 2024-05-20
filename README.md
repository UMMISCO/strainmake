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
| QUAST | [2013](https://doi.org/10.1093/bioinformatics/btt086) | [Yes](https://anaconda.org/bioconda/quast) | https://github.com/ablab/quast | Yes          |

## Binning

| Tool      | Original year                                          | Conda available?                              | Link                                                  | Implemented? |
|:--------- |:------------------------------------------------------ |:--------------------------------------------- |:----------------------------------------------------- |:------------ |
| MetaBAT 2 | [2019](https://doi.org/10.7717%2Fpeerj.7359)           | [Yes](https://anaconda.org/bioconda/metabat2) | https://bitbucket.org/berkeleylab/metabat/src/master/ | Yes          |
| SemiBin2  | [2023](https://doi.org/10.1093/bioinformatics/btad209) | [Yes](https://anaconda.org/bioconda/semibin)  | https://github.com/BigDataBiology/SemiBin             | Yes          |
| VAMB      | [2019](https://doi.org/10.1038/s41587-020-00777-4)     | [Yes](https://anaconda.org/bioconda/vamb)     | https://github.com/RasmussenLab/vamb                  | Yes          |
| COMEBin   | [2024](https://doi.org/10.1038/s41467-023-44290-z)     | [Yes](https://anaconda.org/bioconda/comebin)  | https://github.com/ziyewang/COMEBin                   | No           |

## Bins quality assessment

| Tool    | Original year                                      | Conda available?                             | Link                                 | Implemented? |
| :------ | :------------------------------------------------- | :------------------------------------------- | :----------------------------------- | :----------- |
| CheckM2 | [2023](https://doi.org/10.1038/s41592-023-01940-w) | [Yes](https://anaconda.org/bioconda/checkm2) | https://github.com/chklovski/CheckM2 | Yes          |

## Bins refinement

| Tool    | Original year                                     | Conda available?                             | Link                                        | Implemented? |
|:------- |:------------------------------------------------- |:-------------------------------------------- |:------------------------------------------- |:------------ |
| Binette | [2024](https://doi.org/10.1101/2024.04.20.585171) | [Yes](https://anaconda.org/bioconda/binette) | https://github.com/genotoul-bioinfo/Binette | Yes           |

## Gene prediction, functional annotation and taxonomy classification