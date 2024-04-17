#!/bin/bash

# Specify the directory containing the text files
directory="/home/mc462/basilisk-tx/velocity_data"

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Directory $directory does not exist."
    exit 1
fi

# Loop through each text file in the directory
for file in "$directory"/*.err; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Remove the first 22 lines from the file and save the result to a temporary file
        #tail -n +23 "$file" > "$file.tmp"
	# Remove the last line from the file
	head -n -1 "$file" > "$file.tmp"
        # Replace the original file with the temporary file
        mv "$file.tmp" "$file"
    fi
done

echo "Operation completed."

