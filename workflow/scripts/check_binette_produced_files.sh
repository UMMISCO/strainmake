#!/bin/bash

# check if at least two arguments are passed (check_folder and at least one bins_folder)
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 check_folder bins_folder1 [bins_folder2...]"
    exit 1
fi

# get the check folder path
check_folder=$1
shift

# check if the check_folder contains any .fa files
fa_files=$(find "$check_folder" -type f -name "*.fa")

if [ -z "$fa_files" ]; then
    echo "No .fa files found in $check_folder. Copying and extracting .gz files from bins folders."

    # iterate over each bins folder provided as arguments
    for bins_folder in "$@"; do
        if [ -d "$bins_folder" ]; then
            # find and copy .gz files to check_folder
            find "$bins_folder" -type f -name "*.gz" -exec cp {} "$check_folder" \;
        else
            echo "Warning: $bins_folder is not a valid directory"
        fi
    done

    # extract all .gz files in the check_folder
    for gz_file in "$check_folder"/*.gz; do
        if [ -f "$gz_file" ]; then
            gunzip "$gz_file"
        fi
    done

    # rename all files that are not .fa (.fna...) to .fa
    for file in "$check_folder"/*; do
        if [ -f "$file" ] && [[ "$file" != *.fa ]]; then
            mv "$file" "${file%.*}.fa"
        fi
    done

    echo "All .gz files have been copied and extracted in $check_folder."
else
    echo "Found .fa files in $check_folder. No need to copy and extract .gz files produced by binner(s)."
fi
