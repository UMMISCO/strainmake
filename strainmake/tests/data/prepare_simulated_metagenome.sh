#!/usr/bin/env bash

mkdir -p genomes
mkdir -p data
wget -P genomes -i list_genomes.txt

python3 make_abundance_table.py
python3 make_genome_list_table.py

iss generate --genomes genomes/* --compress \
    --seed 999 --model novaseq --cpus 1 -n 0.4M \
    --output data/fake_illumina

wget https://github.com/bcgsc/NanoSim/raw/master/pre-trained_models/metagenome_ERR3152364_Even.tar.gz
tar -xvzf metagenome_ERR3152364_Even.tar.gz
rm metagenome_ERR3152364_Even.tar.gz

simulator.py metagenome --abun abundance_table_for_nanosim.tsv \
    -gl genomes_list_for_nanosim -dl dna_type_list_for_nanosim.tsv \
    --fastq --seed 999 -c metagenome_ERR3152364_Even/training \
    -b guppy -t 1 --perfect \
    --output data/fake_nanopore

gzip -v data/*.fastq

mv data/fake_illumina_R1.fastq.gz data/fake_illumina_R1.SAMPLE1.fastq.gz
mv data/fake_illumina_R2.fastq.gz data/fake_illumina_R2.SAMPLE1.fastq.gz
mv data/fake_nanopore_sample0_aligned_reads.fastq.gz data/fake_nanopore_sample0_aligned_reads.SAMPLE1.fastq.gz
seqkit fq2fa data/fake_nanopore_sample0_aligned_reads.SAMPLE1.fastq.gz -o data/fake_nanopore_sample0_aligned_reads.SAMPLE1.fasta.gz
