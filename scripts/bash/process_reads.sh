#!/bin/bash

# Script to process ZCH and UCMC raw microbiome reads using kraken_tools (https://github.com/haslamdb/kraken_tools)
# see the link above for installation of kraken2 and associated tools for this analysis

# Usage: ./process_microbiome.sh <input_dir> <output_dir> <metadata_file> <kreport_dir> <bracken_dir>

# Set paths to required databases and tools
KNEADDATA_DB="/path/to/kneaddata_db"  # Replace with your kneaddata reference database path
KRAKEN_DB="/path/to/kraken_db"        # Replace with your Kraken2 database path
BRACKEN_DB="${KRAKEN_DB}/database150mers.kmer_distrib"  # Modify if your Bracken DB is elsewhere

# Parse input arguments
INPUT_DIR=$1  # Directory containing raw FASTQ files
OUTPUT_DIR=$2  # Directory for output
METADATA_FILE=$3  # Your metadata file (AllNICUSampleKeyRevised*.csv)
KREPORT_DIR=${4:-"${OUTPUT_DIR}/kraken_reports"}  # Optional: Custom directory for Kraken reports
BRACKEN_DIR=${5:-"${OUTPUT_DIR}/bracken_output"}  # Optional: Custom directory for Bracken output

# Create output directories
mkdir -p "${OUTPUT_DIR}"

# Run the full pipeline
echo "Starting microbiome analysis pipeline..."
kraken-tools full-pipeline \
    --input-fastq "${INPUT_DIR}"/*.fastq.gz \
    --paired \
    --kneaddata-dbs "${KNEADDATA_DB}" \
    --kraken-db "${KRAKEN_DB}" \
    --bracken-db "${BRACKEN_DB}" \
    --sample-key "${METADATA_FILE}" \
    --output-dir "${OUTPUT_DIR}" \
    --group-col "Location" \
    --threads 8 \
    --log-file "${OUTPUT_DIR}/pipeline.log" \
    --log-level INFO

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "Pipeline encountered errors. Check the log file: ${OUTPUT_DIR}/pipeline.log"
fi