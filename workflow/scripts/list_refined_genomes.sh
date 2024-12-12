#!/bin/bash

# check if at least two arguments are passed (output file and at least one input path)
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 output_file folder1 [folder2 ... folderN]"
    exit 1
fi

# get the output file path
output_file=$1
# remove the first argument to get the input file paths
shift

# initialize an empty array to store the .fa and .fa.gz file paths
fa_files=()

# iterate over each input path provided as arguments
for path in "$@"; do
    # check if the path is a directory
    if [ -d "$path" ]; then
        # find .fa and .fa.gz files in the directory and add them to the array
        while IFS= read -r -d $'\0' file; do
            fa_files+=("$file")
        done < <(find "$path" -type f \( -name "*.fa" -o -name "*.fa.gz" \) -print0)
    elif [ -f "$path" ]; then
        # if the path is a file and ends with .fa or .fa.gz, add it to the array
        if [[ "$path" == *.fa || "$path" == *.fa.gz ]]; then
            fa_files+=("$path")
        fi
    else
        echo "Warning: $path is not a valid file or directory"
    fi
done

# write the paths of .fa files to the output file
printf "%s\n" "${fa_files[@]}" > "$output_file"

echo "Paths of .fa files have been written to $output_file"