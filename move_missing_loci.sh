#!/bin/bash

# Define the path for the missing loci file
list="missing_loci.txt"

# Loop to read each line from the list file
while IFS= read -r file; do
    # Search for the file in the directories above the current folder
    found=$(find ../ -name "$file" -print -quit)
    if [ -n "$found" ]; then
        # If the file is found, copy it to the current folder
        cp "$found" ./
    else
        # If the file is not found, display a message
        echo "File $file not found in the directories above."
    fi
done < "$list"
