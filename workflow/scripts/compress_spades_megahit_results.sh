#!/bin/bash

# checking if the path to the SPAdes/MEGAHIT results directory is provided
if [ -z "$1" ]; then
  echo "Usage: bash compress_spades_megahit_results.sh path/to/assembly/results"
  exit 1
fi

# assigning the first argument (the path) to a variable
SPADES_DIR="$1"

# checking if the directory exists
if [ ! -d "$SPADES_DIR" ]; then
  echo "Error: The directory $SPADES_DIR does not exist."
  exit 1
fi

# moving to the SPAdes results directory
cd "$SPADES_DIR"

echo "Directory content before:"
ls

# creating a tar.gz archive with all files except "assembly.fa.gz"
time tar --use-compress-program="pigz --recursive" -cvf other_files.tar.gz --exclude="assembly.fa.gz" *

# checking if the archive was successfully created
if [ -f "other_files.tar.gz" ]; then
    # if it was successful, deleting all files except "assembly.fa.gz" and "other_files.tar.gz"
    find . ! -name "assembly.fa.gz" ! -name "other_files.tar.gz" -type f,d -print -delete
    echo "Compression complete. Only assembly.fa.gz and other_files.tar.gz remain."
else
    echo "Error: Failed to create other_files.tar.gz."
    exit 1
fi

echo "Directory content after:"
ls
