#!/bin/bash

# Path to the folder where the files and the list are located
directory="/path/to/folder"

# Reads the names from the list and stores them in an array
mapfile -t names < "$directory/loci_list.txt"

# Lists all files in the directory (without full paths)
files=($(ls -A1 $directory))

# Prepares the output file, clearing it if it already exists
> "$directory/missing_loci.txt"

# Loop to check each name in the list of names
for name in "${names[@]}"; do
    # Assumes the name was not found
    found=0
    for file in "${files[@]}"; do
        if [[ "$name" == "$file" ]]; then
            found=1
            break
        fi
    done

    # If the name was not found in any file, prints the name and saves it in the file
    if [[ $found -eq 0 ]]; then
        echo "$name" >> "$directory/missing_loci.txt"
    fi
done
