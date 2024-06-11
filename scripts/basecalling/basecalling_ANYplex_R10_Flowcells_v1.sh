#!/bin/bash

### Script for basecalling with Bonito
### Author: Yuri Malina
### Contact: ymalina@berkeley.edu
### Organization: Meyer Lab at UC Berkeley

### REQUIREMENTS:
# Install Dorado instructions:
# https://community.nanoporetech.com/posts/how-can-install-dorado-via
# https://github.com/nanoporetech/dorado
# Remember to download basecalling models if running for the first time with following command:
# dorado download --model all

### Define variables:
INDIR="/Data1/seq_data/BN_96DPY27Deg_Fiber_Hia5_MCviPI_05_24_24/no_sample/combined_pod5"
OUTFOLDER="/Data1/git/meyer-nanopore/scripts/basecalling/"
OUTFILE="${OUTFOLDER}mod_mappings.bam"
SORTED_OUTFILE="${OUTFOLDER}mod_mappings.sorted.bam"
REF="/Data1/reference/ws235.mmi" # this is c. elegans reference file, change to appropriate reference
COMPLETED_OUTFILE="${OUTFOLDER}mod_mappings.sorted.completed.bam" # filename for completed modification encoding
COMPLETED_SORTED_OUTFILE="${OUTFOLDER}mod_mappings.sorted.completed.sorted.bam" # filename for completed modification encoding

# Set parameters
DEMULTIPLEX=false
RESUME=false
KIT_NAME="SQK-NBD114-24"

# Create output folder if it doesn't exist
if [ ! -d "$OUTFOLDER" ]; then
    mkdir -p "$OUTFOLDER"
    echo "Output folder created: $OUTFOLDER"
else
    echo "Output folder already exists: $OUTFOLDER"
fi

echo "Input Directory: $INDIR"
echo "Output File: $OUTFILE"
echo "Sorted Output File: $SORTED_OUTFILE"
echo "Reference: $REF"

### PICK APPROPRIATE MODEL
# 5khz-data
# dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mC@v2 (5 kHz)
# dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v2 (5 kHz)
# res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2 (4 KHz)
# res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2 (4 KHz)

cd /Data1/software/dorado/bin

BASECALLER_COMMAND="./dorado basecaller \
/Data1/software/rerio/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
$INDIR \
--reference $REF \
--modified-bases-threshold 0 \
--device cuda:2,3 \
--modified-bases 5mC 6mA"

if [ "$RESUME" = true ]; then
    BASECALLER_COMMAND="$BASECALLER_COMMAND \
    --resume-from $OUTFILE \
    > $COMPLETED_OUTFILE"
    OUTFILE="$COMPLETED_OUTFILE"
    SORTED_OUTFILE="$COMPLETED_SORTED_OUTFILE"
    echo "Completed Output File: $OUTFILE"
    echo "Completed Sorted Output File: $SORTED_OUTFILE"
else
    BASECALLER_COMMAND="$BASECALLER_COMMAND \
    > $OUTFILE"
fi

if [ "$DEMULTIPLEX" = true ]; then
    BASECALLER_COMMAND="$BASECALLER_COMMAND \
    --kit-name $KIT_NAME"
fi

eval $BASECALLER_COMMAND

if [ "$DEMULTIPLEX" = true ]; then
    /home/grid/Documents/Basecalling/dorado-0.6.0-linux-x64/bin/dorado demux \
    --no-trim \
    -o $OUTFOLDER \
    --no-classify $OUTFILE
fi

# Iterate over each .bam file in the folder and sort it. This will handle singleplex and multiplex
for bam_file in "$OUTFOLDER"/*.bam; do
  # Check if there are .bam files in the folder
  if [ -f "$bam_file" ]; then
    # Get the filename without extension
    filename=$(basename -- "$bam_file")
    filename_no_ext="${filename%.*}"
    
    # Sort the .bam file
    samtools sort -o "$OUTFOLDER/$filename_no_ext.sorted.bam" "$bam_file"
    
    # Optional: Index the sorted .bam file
    samtools index "$OUTFOLDER/$filename_no_ext.sorted.bam"
    
    echo "Sorted and indexed $filename"
  else
    echo "No .bam files found in $OUTFOLDER"
  fi
done

# Correct .bam file if unecessary '.' characters are found in modification .bam file encodings.
# This is likely unecessary. Leaving commented out.
# samtools view -h $SORTED_OUTFILE | \
#             sed 's/A+a./A+a/g' | \
#             sed 's/C+m./C+m/g' | \
#             samtools view -o $COMPLETED_OUTFILE'
            
# Index .bam file
# samtools index $COMPLETED_OUTFILE
