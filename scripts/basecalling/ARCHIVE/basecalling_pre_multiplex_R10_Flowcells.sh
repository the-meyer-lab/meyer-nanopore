### Script for basecalling with Bonito
### Author: Yuri Malina
### Contact: ymalina@berkeley.edu
### Organization: Meyer Lab at UC Berkeley

#!/bin/bash

### REQUIREMENTS:
# Install Dorado instructions:
# https://community.nanoporetech.com/posts/how-can-install-dorado-via
# https://github.com/nanoporetech/dorado
# Remember to download basecalling models if running for the first time with following command:
# dorado download --model all

# Define variables:
#INDIR="/Data1/seq_data/AI_N2_dimelo_antiDPY27_mpx_8_19_23/pod5_pass"
INDIR="/Data1/seq_data/AH_N2_SDC2aid_AuxRem_fiberseq_8_19_23/combined_pod5"
REF="/Data1/reference/ws235.mmi" # this is c. elegans reference file, change to appropriate reference
# Corrected file name, uncomment if uncommenting the corrected bam asociated lines below. This should not be necessary.
#CORRECTED_OUTFILE="${OUTFOLDER}mod_mappings2.sorted.corrected.bam" # filename for corrected modification encoding

echo "Input Directory: $INDIR"
echo "Reference: $REF"

# For each folder starting with barcode in INDIR, create a new folder within called "basecalls"
# and run dorado basecaller
for folder in $INDIR/barcode*; do
  # set variable to barcode number
  BARCODE=$(basename $folder)

  echo "Starting on folder: $folder"
  # if $folder/baecalls does not exist
  if [ ! -d "$folder/basecalls" ]; then
    # create $folder/basecalls
    mkdir -p "$folder/basecalls"
    echo "Output folder created: $folder/basecalls"
  else
    echo "Output folder already exists: $folder/basecalls"
  fi
  # Define output file name
  OUTFOLDER="$folder/basecalls/"
  # Define output file path and name as $OUTFOLDER + $BARCODE + "mod_mappings.bam"
  OUTFILE="${OUTFOLDER}${BARCODE}.mod_mappings.bam"
  SORTED_OUTFILE="${OUTFOLDER}${BARCODE}.mod_mappings.sorted.bam"

  # <<4kHzData
  #dorado basecaller \
  #/Data1/software/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 \
  #$INDIR \
  #--reference $REF \
  #--modified-bases-threshold 0 \
  #--modified-bases-models \
  #  /Data1/software/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2,/Data1/software/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2\
  #> $OUTFILE
  #4kHzData

  # <<5khz-data
  # For each folder starting with barcode in INDIR run the following code:
  dorado basecaller \
  /Data1/software/rerio/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
  $folder \
  --device cuda:0,1,2,3,4,5,6,7 \
  --reference $REF \
  --modified-bases-threshold 0 \
  --modified-bases 5mC 6mA \
  > $OUTFILE
  # 5khz-data

  # Sort .bam file
  samtools sort -o $SORTED_OUTFILE $OUTFILE

  # Index .bam file
  samtools index $SORTED_OUTFILE
done

# Note: The above basecalles 5mC and m6A base modofications. Either can be removed if desired.
# Recommend using modified base threshold of 0, and setting base modification thresholds in downstream analysis pipeline
# In order to avoid having to re-basecall.
# dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mC@v2 (5 kHz)
# dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v2 (5 kHz)
# res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2 (4 KHz)
# res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2 (4 KHz)



# Correct .bam file if unecessary '.' characters are found in modification .bam file encodings.
# This is likely unecessary. Leaving commented out.
<<optional
samtools view -h $SORTED_OUTFILE | \
            sed 's/A+a./A+a/g' | \
            sed 's/C+m./C+m/g' | \
            samtools view -o $CORRECTED_OUTFILE'
            
# Index .bam file
samtools index $CORRECTED_OUTFILE
optional
