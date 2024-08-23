This folder contains scripts to test parts of the pipeline.

In order to have Illumina short reads data we generated them using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) and in order to have
Nanopore long reads data we generated them using [NanoSim](https://github.com/bcgsc/NanoSim).
Data can be found in [`data`](data/).
To reproduce the generated data: 

```sh
conda install -c bioconda -c conda-forge insilicoseq==2.0.1 nanosim==3.1.0 scikit-learn==0.22.1 numpy==1.21.5

sh prepare_simulated_metagenome.sh
```

If we already have the sequencing data but need to generate the metadata tables used by the the pipeline to identify where are the FASTQ, run;

```sh
python3 generate_metadata_table.py data/ SR && mv metadata.tsv "1._SR_only"
python3 generate_metadata_table.py data/ LR && mv metadata.tsv "2._LR_only"
python3 generate_metadata_table.py data/ all && mv metadata.tsv "3._ALL_seq"

cp "1._SR_only/metadata.tsv" "5._LR_fail_if_SR_only"
cp "2._LR_only/metadata.tsv" "4._SR_fail_if_LR_only"```

Then, to run the tests:

```sh
python3 test.py
```