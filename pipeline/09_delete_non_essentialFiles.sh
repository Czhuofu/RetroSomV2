#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 -t TARGET_DIR"
    echo
    echo "  -t, --target-dir   Specify the target directory"
}

# Parse command-line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--target-dir) TARGET_DIR="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check if TARGET_DIR is specified
if [ -z "$TARGET_DIR" ]; then
    echo "Error: Target directory not specified."
    show_help
    exit 1
fi


# Check if TARGET_DIR exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' does not exist."
    exit 1
fi

# Change to the target directory
cd "$TARGET_DIR" || exit
SUB_DIRS=($(find . -maxdepth 1 -mindepth 1 -type d))


# Loop through the array of subdirectories
for dir in "${SUB_DIRS[@]}"; do
    echo "Processing directory: $dir"
    
    # step-1: delete exonerate files"
    ls $dir/script/exonerate*
    rm $dir/*/script/exonerate*
    
    # step-2: delete seg files"
    ls $dir/retro_v1_*/seg*
    rm $dir/*/retro_v1_*/seg*
    
    # step-3: delete bam files"
    ls $dir/align/*.final.bam
    rm $dir/*/align/*.final.bam
    

done;

