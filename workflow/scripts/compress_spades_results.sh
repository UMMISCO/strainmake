#!/bin/bash

# checking if the path to the SPAdes results directory is provided
if [ -z "$1" ]; then
  echo "Usage: bash compress_spades_results.sh path/to/spades/results"
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
cd "$SPADES_DIR" || exit 1

# creating a tar.gz archive with all files except "assembly.fa.gz"
tar --use-compress-program="pigz --recursive" -cvf other_files.tar.gz --exclude="assembly.fa.gz" *

# checking if the archive was successfully created
if [ -f "other_files.tar.gz" ]; then
    # if it was successful, deleting all files except "assembly.fa.gz" and "other_files.tar.gz"
    find . ! -name "assembly.fa.gz" ! -name "other_files.tar.gz" -type f -delete
    echo "Compression complete. Only assembly.fa.gz and other_files.tar.gz remain."
else
    echo "Error: Failed to create other_files.tar.gz."
    exit 1
fi
