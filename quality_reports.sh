#!/bin/bash

# Function to create necessary directories if they don't exist
create_directories() {
    # Prompt the user for the desired folder name for FastQC reports
    echo "Enter a name for the folder where the FastQC reports will be saved:"
    read -e -r FOLDER_NAME  # Read user input for folder name

    # If no name is provided, default to "FastQC_Reports_<timestamp>"
    if [ -z "$FOLDER_NAME" ]; then
        FOLDER_NAME="FastQC_Reports_$(date +%Y%m%d%H%M%S)"  # Generate a timestamp-based folder name
    fi

    # Define the directories
    OUTPUT_DIR="/home/stp/Limia/CMV_Project/CMV_Data/fastqc_reports/$FOLDER_NAME"
    TEMP_DIR="/home/stp/Limia/CMV_Project/CMV_Data/temp_fastq"

    # Create directories if they don't exist
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$TEMP_DIR"

    echo "Directories verified or created: $OUTPUT_DIR, $TEMP_DIR"
}

# Function to process each file with FastQC
process_file() {
    # Get the file path from the first argument
    local file=$1
    # Extract the base name of the file without the .gz extension
    local base_name=$(basename "$file" .gz)
    
    # Unzip the file and store the content in TEMP_DIR with the original file name (minus the .gz extension)
    gunzip -c "$file" > "$TEMP_DIR/$base_name"

    # Run FastQC on the unzipped file, outputting the results to OUTPUT_DIR
    fastqc -o "$OUTPUT_DIR" "$TEMP_DIR/$base_name"
    
    # Remove the temporary unzipped file to free up space
    rm "$TEMP_DIR/$base_name"
    
    echo "Processed and cleaned: $base_name"
}

# Main script logic
main() {
    # Ask the user for the input directory
    echo "Where is your directory located?"
    read -e -r INPUT_DIR  # Read the user input for the directory

    # Check if the input directory was provided
    if [ -z "$INPUT_DIR" ]; then
        echo "No directory specified. Exiting..."
        exit 1  # Exit the script with an error status
    fi

    # Create necessary directories
    create_directories  # Call the function to create OUTPUT_DIR and TEMP_DIR

    # Export the functions and variables for parallel processing
    export -f process_file  # Export the process_file function so it can be used with parallel
    export TEMP_DIR OUTPUT_DIR  # Export the directories for use in parallel processes

    # Run file processing in parallel
    echo "Starting parallel processing of files in $INPUT_DIR"
    find "$INPUT_DIR" -name "*.fastq.gz" | parallel process_file {}  # Find all .fastq.gz files and process them in parallel

    # Clean up temporary directory
    rmdir "$TEMP_DIR"  # Remove the TEMP_DIR after processing

    # Inform the user that the quality check is complete and where the reports are saved
    echo "Quality check complete. Reports are saved in $OUTPUT_DIR."
}

# Run the main function and pass any arguments provided to this script
main
