#!/bin/bash

# ---
# Title: filter_merged_vcf.sh
# Date: 2025
# Author: Adam, Carolina de Lima
# Purpose: filter multi-sample VCF file
# ---

# Default values for parameters
MAX_MISSING=1
SD_THRESHOLD=3
AP_THRESHOLD=0.6
OUTPUT_DIR="."

# Function to display usage information
usage() {
    echo "Usage: $0 [-i input_vcf] [-m max_missing] [-s sd_threshold] [-a ap_threshold] [-o output_dir]"
    echo "  -i  Input VCF file (required)"
    echo "  -m  Max missing data fraction (default: 1 - no missing data)"
    echo "  -s  Minimum spanning depth (default: 3)"
    echo "  -a  Minimum purity score (default: 0.6)"
    echo "  -o  Output directory (default: current directory)"
    exit 1
}

# Parse command-line arguments
while getopts "i:m:s:a:o:" opt; do
    case "$opt" in
        i) INPUT_VCF="$OPTARG" ;;
        m) MAX_MISSING="$OPTARG" ;;
        s) SD_THRESHOLD="$OPTARG" ;;
        a) AP_THRESHOLD="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

# Ensure required arguments are provided
if [ -z "$INPUT_VCF" ]; then
    echo "Error: Input VCF file is required."
    usage
fi

# Get the base name of the input file (without path and extension)
BASE_NAME=$(basename "$INPUT_FILE" .vcf)

# Filter by missing data
vcftools --vcf "$INPUT_FILE" --max-missing "$MAX_MISSING" --recode --recode-INFO-all --out "${OUTPUT_FILE}_m${MAX_MISSING}"

# Filter alleles by spanning depth
bcftools view -i "MIN(FORMAT/SD) >= $SD_THRESHOLD" "${OUTPUT_FILE}_m${MAX_MISSING}.recode.vcf" -o "${OUTPUT_FILE}_m${MAX_MISSING}_${SD_THRESHOLD}.vcf"

# Filter alleles by purity score
bcftools view -i "MIN(FORMAT/AP) >= $AP_THRESHOLD" "${OUTPUT_FILE}_m${MAX_MISSING}_${SD_THRESHOLD}.vcf" -o "${OUTPUT_FILE}_m${MAX_MISSING}_${SD_THRESHOLD}_${AP_THRESHOLD}.vcf"

echo "Processing completed. Final output file: ${OUTPUT_FILE}_m${MAX_MISSING}_${SD_THRESHOLD}_${AP_THRESHOLD}.vcf"
