#!/bin/bash

# Path to the software tools
MODKIT_ADJUST="/Data1/software/modkit/modkit adjust-mods"
MODKIT_UPDATE="/Data1/software/modkit/modkit update-tags"

# Number of threads
THREADS=256
MAX_BG_PROCS=4

# Array of input files
declare -a filepaths=(
  "/Data1/seq_data/TubeY9B_N2_fiberseq_timec_mpx_3_21_2023/demux/mod_mappings_barcode_07.bam"
  "/Data1/seq_data/TubeY9B_N2_fiberseq_timec_mpx_3_21_2023/demux/mod_mappings_barcode_08.bam"
  "/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/basecalls/m6A/mod_mappings.sorted.bam"
  "/Data1/seq_data/TubeAB_FiberSeq_TimeC_N2_021Aux_4_10_23/basecalls/m6A_full/demux/mod_mappings_barcode_05.bam"
  "/Data1/seq_data/TubeAB_FiberSeq_TimeC_N2_021Aux_4_10_23/basecalls/m6A_full/demux/mod_mappings_barcode_04.bam"
  "/Data1/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/basecalls/mod_mappings.sorted.m6Aonly.bam"
  "/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/basecalls/m6A/mod_mappings.sorted.bam"
  "/Data1/seq_data/TubeAB_FiberSeq_TimeC_N2_021Aux_4_10_23/basecalls/m6A_full/demux/mod_mappings_barcode_10.bam"
  "/Data1/seq_data/TubeAB_FiberSeq_TimeC_N2_021Aux_4_10_23/basecalls/m6A_full/demux/mod_mappings_barcode_09.bam"
  "/Data1/seq_data/TubeAB_FiberSeq_TimeC_N2_021Aux_4_10_23/basecalls/m6A_full/demux/mod_mappings_barcode_06.bam"
  "/Data1/seq_data/TubeAB_FiberSeq_TimeC_N2_021Aux_4_10_23/basecalls/m6A_full/demux/mod_mappings_barcode_01.bam"
  "/Data1/seq_data/TubeY9B_N2_fiberseq_timec_mpx_3_21_2023/demux/mod_mappings_barcode_01.bam"
)

# Check that all input files exist
for filepath in "${filepaths[@]}"; do
    if [[ ! -f "$filepath" ]]; then
        echo "Error: File '$filepath' not found. Exiting."
        exit 1
    fi
done

# Initialize counter for background processes
bg_procs=0

# Loop through each file in the array
for filepath in "${filepaths[@]}"; do
    {
        # Extract the file name and directory
        filename=$(basename -- "$filepath")
        dir=$(dirname -- "$filepath")
        extension="${filename##*.}"
        filename_noext="${filename%.*}"

        # Create output file names
        tmp_output="$dir/${filename_noext}_tmp.$extension"
        final_output="$dir/${filename_noext}_Ya_ambig.$extension"

        # Execute COMMAND 1
        $MODKIT_ADJUST --convert Y a --convert Z m --threads $THREADS "$filepath" "$tmp_output"

        # Execute COMMAND 2
        $MODKIT_UPDATE --threads $THREADS --mode ambiguous "$tmp_output" "$final_output"

        # Optional: Remove temporary file
        rm "$tmp_output"
    } &

    # Increment counter and check for max background processes
    ((++bg_procs))
    if [[ $bg_procs -ge $MAX_BG_PROCS ]]; then
        wait -n  # Wait for any background process to complete
        ((--bg_procs))
    fi
done

# Wait for all remaining background processes to complete
wait


# Initialize an empty array to hold output files for QC
declare -a output_files=()

# Quality Control and Reporting
echo "Quality Control and Reporting:"
echo "---------------------------------------"

# Loop through each file in the array to collect output file names
for filepath in "${filepaths[@]}"; do
    # Extract the file name and directory
    filename=$(basename -- "$filepath")
    dir=$(dirname -- "$filepath")
    extension="${filename##*.}"
    filename_noext="${filename%.*}"

    # Create final output file names
    final_output="$dir/${filename_noext}_Ya_ambig.$extension"

    # Add to the list of output files for QC
    output_files+=("$final_output")
done

# Loop through each output file for QC
for filepath in "${output_files[@]}"; do
    # Check if file was created
    if [[ -f "$filepath" ]]; then
        # File size
        filesize=$(ls -lh "$filepath" | awk '{print $5}')

        # Number of rows (assuming BAM format)
        num_rows=$(samtools view -c "$filepath")

        echo "File: $filepath"
        echo "  - Size: $filesize"
        echo "  - Rows: $num_rows"
    else
        echo "File: $filepath was not created."
    fi
    echo "---------------------------------------"
done