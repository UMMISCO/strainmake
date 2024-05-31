#!/usr/bin/env bash 

# retrieving the reads from the ENA
RUN=ERR1855257
RUN2=ERR1855268
RUN3=ERR1855286
RUN4=ERR1855267

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/007/${RUN}/${RUN}_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/007/${RUN}/${RUN}_2.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/008/${RUN2}/${RUN2}_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/008/${RUN2}/${RUN2}_2.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/006/${RUN3}/${RUN3}_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/006/${RUN3}/${RUN3}_2.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/007/${RUN4}/${RUN4}_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR185/007/${RUN4}/${RUN4}_2.fastq.gz

mv ${RUN}_1.fastq.gz data/${RUN}_1.fastq.gz
mv ${RUN}_2.fastq.gz data/${RUN}_2.fastq.gz
mv ${RUN2}_1.fastq.gz data/${RUN2}_1.fastq.gz
mv ${RUN2}_2.fastq.gz data/${RUN2}_2.fastq.gz
mv ${RUN3}_1.fastq.gz data/${RUN3}_1.fastq.gz
mv ${RUN3}_2.fastq.gz data/${RUN3}_2.fastq.gz
mv ${RUN4}_1.fastq.gz data/${RUN4}_1.fastq.gz
mv ${RUN4}_2.fastq.gz data/${RUN4}_2.fastq.gz