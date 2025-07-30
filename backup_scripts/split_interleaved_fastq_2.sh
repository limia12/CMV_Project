#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_dir>"
    exit 1
fi

# Assign arguments to variables
input_dir="$1"

# Create the output directory 'split_interleaved' within the input directory
output_dir="$input_dir/split_interleaved"
mkdir -p "$output_dir"

# Function to split a single fastq file
split_fastq() {
    fastq_file="$1"
    output_dir="$2"

    # Get the base name of the file (without extension)
    base_name=$(basename "$fastq_file" .fastq)

    # Create the R1 and R2 fastq files
    r1_file="$output_dir/${base_name}_R1.fastq"
    r2_file="$output_dir/${base_name}_R2.fastq"

    # Initialize empty R1 and R2 files
    > "$r1_file"
    > "$r2_file"

    # Read the interleaved fastq file and split every 4 lines
    while read -r line1 && read -r line2 && read -r line3 && read -r line4; do
        # Write the first 4 lines (R1) to the R1 file
        echo -e "$line1\n$line2\n$line3\n$line4" >> "$r1_file"
        
        # Read the next 4 lines (R2) and write them to the R2 file
        read -r line1 && read -r line2 && read -r line3 && read -r line4
        echo -e "$line1\n$line2\n$line3\n$line4" >> "$r2_file"
    done < "$fastq_file"

    echo "Finished splitting $fastq_file into $r1_file and $r2_file"
}

export -f split_fastq

# Find all the fastq files and process them with xargs in parallel
find "$input_dir" -name "*.fastq" | xargs -n 1 -P 16 -I {} bash -c 'split_fastq "$@"' _ {} "$output_dir"
