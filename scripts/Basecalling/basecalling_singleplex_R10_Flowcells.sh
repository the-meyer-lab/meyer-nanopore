### Script for basecalling with Bonito
#!/bin/bash

# Install Dorado instructions:
# https://community.nanoporetech.com/posts/how-can-install-dorado-via
# https://github.com/nanoporetech/dorado

# Define variables:
INDIR="/Data1/seq_data/TubeAD1_N2_fiberseq_6_13_23/pod5/"
OUTFOLDER="/Data1/seq_data/TubeAD1_N2_fiberseq_6_13_23/mod_basecalls/"
OUTFILE="${OUTFOLDER}mod_mappings2.bam"
SORTED_OUTFILE="${OUTFOLDER}mod_mappings2.sorted.bam"
CORRECTED_OUTFILE="${OUTFOLDER}mod_mappings2.sorted.corrected.bam" # filename for corrected modification encoding
REF="/Data1/reference/ws235.mmi" # this is c. elegans reference file, change to appropriate reference

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

  > calls.bam

dorado basecaller \
/Data1/software/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 \
$INDIR \
--reference $REF \
--modified-bases-threshold 0 \
--modified-bases-models \
  /Data1/software/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2,/Data1/software/rerio/dorado_models/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2 \
> $OUTFILE
# Note: The above basecalles 5mC and m6A base modofications. Either can be removed if desired.

# Sort .bam file
samtools sort -o $SORTED_OUTFILE $OUTFILE

# Correct .bam file if unecessary '.' characters are found in modification .bam file encodings.
# This is likely unecessary. Leaving commented out.
<<optional
samtools view -h $SORTED_OUTFILE | \
            sed 's/A+a./A+a/g' | \
            sed 's/C+m./C+m/g' | \
            samtools view -o $CORRECTED_OUTFILE'
optional

# Index .bam file
samtools index $CORRECTED_OUTFILE