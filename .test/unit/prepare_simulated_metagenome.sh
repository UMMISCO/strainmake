# downloading genomes of 5 marine bacteria
mkdir -p genomes
mkdir -p data
wget -P genomes -i list_genomes.txt

# generate a table with bacteria abundance in the metagenome for NanoSim
python3 make_abundance_table.py
# generate genomes list for simulators
python3 make_genome_list_table.py

# generating short reads for a metagenome containing those 5 species
iss generate --genomes genomes/* --compress \
    --seed 999 --model novaseq --cpus 1 -n 0.4M \
    --output data/fake_illumina

# downloading the NanoSim model for metagenomic data
wget https://github.com/bcgsc/NanoSim/raw/master/pre-trained_models/metagenome_ERR3152364_Even.tar.gz
tar -xvzf metagenome_ERR3152364_Even.tar.gz
rm metagenome_ERR3152364_Even.tar.gz

# generating long reads for a metagenome containing those 5 species
simulator.py metagenome --abun abundance_table_for_nanosim.tsv \
    -gl genomes_list_for_nanosim -dl dna_type_list_for_nanosim.tsv \
    --fastq --seed 999 -c metagenome_ERR3152364_Even/training \
    -b guppy -t 1 --perfect \
    --output data/fake_nanopore

gzip -v data/*.fastq

mv data/fake_illumina_R1.fastq.gz data/fake_illumina_R1.SAMPLE1.fastq.gz
mv data/fake_illumina_R2.fastq.gz data/fake_illumina_R2.SAMPLE1.fastq.gz
mv data/fake_nanopore_sample0_aligned_reads.fastq.gz data/fake_nanopore_sample0_aligned_reads.SAMPLE1.fastq.gz

# ensure generated files are correct
echo "MD5 sum checking"
md5sum --check fastq.md5

rm data/*.vcf
rm data/*.txt
rm data/*.fasta
rm data/.iss*
rm data/*error_profile
rm -r metagenome_ERR3152364_Even/
rm *.tsv genomes_list_for_nanosim

# prepare folders for testing case
mkdir -p "1._SR_only"
mkdir -p "2._LR_only"
mkdir -p "3._ALL_seq"
mkdir -p "4._SR_fail_if_LR_only"
mkdir -p "5._LR_fail_if_SR_only"
mkdir -p "6._SR_not_paired"

python3 generate_metadata_table.py data/ SR && mv metadata.tsv "1._SR_only"
python3 generate_metadata_table.py data/ LR && mv metadata.tsv "2._LR_only"
python3 generate_metadata_table.py data/ all && mv metadata.tsv "3._ALL_seq"

cp "1._SR_only/metadata.tsv" "5._LR_fail_if_SR_only"
cp "2._LR_only/metadata.tsv" "4._SR_fail_if_LR_only"

grep -v R2 "3._ALL_seq/metadata.tsv" > "6._SR_not_paired/metadata.tsv"