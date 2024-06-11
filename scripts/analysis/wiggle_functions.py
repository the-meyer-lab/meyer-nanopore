### Used to handle .wig file for nucleosome repeat length calculations.
# Must be paired with NRL_Analysis.ipynb
import sys
import io #used for interacting with in memory text file
from collections import defaultdict
import os.path
import pysam
import pandas as pd

# if output folder does not exist, create it.
#if not os.path.exists("output/"):
#    os.makedirs("output/")

'''m6a_fwd_H_name = "/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/wig/H_modified_bases.6mA.fwd_strand.wig"
m6a_rev_H_name = "/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/wig/H_modified_bases.6mA.rev_strand.wig"
m6a_fwd_D_name = "/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/wig/D_modified_bases.6mA.fwd_strand.wig"
m6a_rev_D_name = "/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/wig/D_modified_bases.6mA.rev_strand.wig"


# Take in 4 file names
wig_file_name1 = m6a_fwd_H_name
wig_file_name2 = m6a_rev_H_name
wig_file_name3 = m6a_fwd_D_name
wig_file_name4 = m6a_rev_D_name

wig_file1_obj = open(wig_file_name1,"r")
wig_file1 = io.StringIO(wig_file1_obj.read())
wig_file2_obj = open(wig_file_name2,"r")
wig_file2 = io.StringIO(wig_file2_obj.read())
wig_file3_obj = open(wig_file_name3,"r")
wig_file3 = io.StringIO(wig_file3_obj.read())
wig_file4_obj = open(wig_file_name4,"r")
wig_file4 = io.StringIO(wig_file4_obj.read())

# Take in output file name
out_file_name = "D_modified_bases.6mA.merged_strand_normH.wig"'''

# Define function that normalizes wig file to max m6A probability
def norm_wig(wig_file1):
  # Wig file 2 correction factor
  correction_factor = 1/256

  norm_file = io.StringIO()
  for line in wig_file1:
    if line[0].isdigit():
      position, value = line.strip().split()
      if (float(value) * correction_factor) > 1:
        normalized_value = str(1)
      else:
        normalized_value = str(float(value) * correction_factor)
      norm_file.write(position + ' ' + normalized_value + '\n')
    else:
      norm_file.write(line)
  norm_file.seek(0)
  return norm_file
  
# Define function that takes in two wig file names, and returns a iostring object of the normalized second file.
def normalize_wigs(wig_file1, wig_file2):
  print("NORMALIZING", wig_file1, " TO ", wig_file2)
  #with open(wig_file_name1, 'r') as wig_file1, open(wig_file_name2, 'r') as wig_file2:
  # Uncomment below to process limited number of lines
  #lines1 = [next(wig_file1) for _ in range(200)]
  #lines2 = [next(wig_file2) for _ in range(200)]
  #wig_file1 = io.StringIO(''.join(lines1))
  #wig_file2 = io.StringIO(''.join(lines2))
  #wig_file1.seek(0)
  #wig_file2.seek(0)

  #Determine average value in first wig file
  values1 = []
  for line in wig_file1:
    if line[0].isdigit():
      position, value = line.strip().split()
      values1.append(float(value))
  avg_value1 = sum(values1) / len(values1)
  wig_file1.seek(0)
  print(wig_file_name1, " average meth value = ", avg_value1)

  # Determine average value in second wig file
  values2 = []
  for line in wig_file2:
    columns = line.split()
    if line[0].isdigit():
      position, value = line.strip().split()
      values2.append(float(value))
  avg_value2 = sum(values2) / len(values2)
  wig_file2.seek(0)
  print(wig_file_name2, " average meth value = ", avg_value2)

  # Wig file 2 correction factor
  correction_factor = avg_value1/avg_value2

  norm_file = io.StringIO()
  for line in wig_file2:
    if line[0].isdigit():
      position, value = line.strip().split()
      if (float(value) * correction_factor) > 1:
        normalized_value = str(1)
      else:
        normalized_value = str(float(value) * correction_factor)
      norm_file.write(position + ' ' + normalized_value + '\n')
    else:
      norm_file.write(line)
  norm_file.seek(0)
  return norm_file

  # Determine average of wig 1
  # Determine average of wig 2
  # Multiply all values of wig2
  # Output normalized file.

# Define function that takes in two wig file iostring objects and returns an iostring object of a merged wig
def merge_wig(wig_file1, wig_file2):
  # Uncomment below to process limited number of lines
  #lines1 = [next(wig_file1) for _ in range(200)]
  #lines2 = [next(wig_file2) for _ in range(200)]
  #wig_file1 = io.StringIO(''.join(lines1))
  #wig_file2 = io.StringIO(''.join(lines2))
  wig_file1.seek(0)
  wig_file2.seek(0)

  print("MERGING", wig_file1, "INTO ", wig_file2)
  # Initialize a dictionary to store the positions and values for each chromosome
  chrom_data = defaultdict(list)

  # Initialize the current chromosome
  chrom = None
  # Iterate over each line in the file
  print("GATHERING DATA FROM FIRST FILE")
  for line in wig_file1:
    # If the line starts with "track" or "variableStep", extract the chromosome information
    if line.startswith('track'):
      track_line = line
    elif line.startswith('variableStep'):
      chrom = None
      # Split the line into key-value pairs
      kv_pairs = line.split()
      # Iterate over the key-value pairs
      for kv in kv_pairs:
        if kv != 'variableStep':
          # Split the key-value pair into a key and value
          k, v = kv.split('=')
          # If the key is "chrom", set the current chromosome
          if k == 'chrom':
            chrom = v
            print("GATHERING DATA FROM: ", chrom)
    elif chrom is not None and line.split()[0].isdigit():
      chrom_data[chrom].append((int(line.split()[0]), float(line.split()[1])))


  print("GATHERING DATA FROM SECOND FILE")
  for line in wig_file2:
    # If the line starts with "track" or "variableStep", extract the chromosome information
    if line.startswith('track'):
      pass
    elif line.startswith('variableStep'):
      chrom = None
      # Split the line into key-value pairs
      kv_pairs = line.split()
      # Iterate over the key-value pairs
      for kv in kv_pairs:
        if kv != 'variableStep':
          # Split the key-value pair into a key and value
          k, v = kv.split('=')
          # If the key is "chrom", set the current chromosome
          if k == 'chrom':
            chrom = v
            print("GATHERING DATA FROM: ", chrom)
        # If the line has at least two elements and the first column is a number, store the position and value for the chromosome
    elif chrom is not None and line.split()[0].isdigit():
      chrom_data[chrom].append((int(line.split()[0]), float(line.split()[1])))

  # Sort chrom positions for each chrom
  for chrom, positions in chrom_data.items():
    chrom_data[chrom] = sorted(positions, key=lambda x: x[0])

  merged_file = io.StringIO()
  # Write the "track" line
  merged_file.write(track_line)
  # Iterate over the chromosomes in the dictionary
  for chrom, positions in chrom_data.items():
    # Write the "variableStep" line
    merged_file.write(f'variableStep chrom={chrom} span=1\n')
    # Iterate over the positions and values for the chromosome
    for pos, val in positions:
      # Write the position and value to the output file
      merged_file.write(f'{pos}\t{val}\n')
  merged_file.seek(0)
  return merged_file

# Define function that returns a wig file with values equal to 1 - the current value
def invert_wig(wig_file,max_v):
  # Read in .wig file from command line argument
  # Get the total number of lines in the file
  num_lines = sum(1 for line in wig_file)
  # Seek back to the beginning of the file
  wig_file.seek(0)
  wig_file_lines = wig_file.readlines()
  wig_file.seek(0)
  # Initialize a counter to track progress
  line_counter = 0

  out_file = io.StringIO()
  for i in range(0,len(wig_file_lines)):
    line = wig_file_lines[i]
    # Split the line into columns
    columns = line.split()

    if columns[0].isdigit():
      # If the line is a data line, print it
      current_value = int(columns[0])
      value = round((1 - float(columns[1]))*max_v,2)
      out_file.write(f'{columns[0]}\t{value}\n')

    # If the line is a metadata line, print it
    else:
      out_file.write(f'{line}')

    # Increment the counter and display progress
    line_counter += 1
    progress = line_counter / num_lines * 100
    if line_counter % int(num_lines/10) == 0:
      print(f'Progress: {progress:.2f}%', end='\r')
      # Print a newline after the progress message
      print()
  out_file.seek(0)
  return out_file

# Define function that fills in any gaps in a wig file with a max_value
def fill_wig(wig_file1,max_v):
  # Fill in gaps
  print(f'Filling wig: {wig_file1}')

  # Get the total number of lines in the file
  num_lines = sum(1 for line in wig_file)
  # Seek back to the beginning of the file
  wig_file.seek(0)
  wig_file_lines = wig_file.readlines()
  wig_file.seek(0)
  # Initialize a counter to track progress
  line_counter = 0

  out_file = io.StringIO()
  for i in range(0,len(wig_file_lines)):
    line = wig_file_lines[i]
    # Split the line into columns
    columns = line.split()

    if columns[0].isdigit():
      # If the line is a data line, print it
      current_value = int(columns[0])
      value = columns[1]
      out_file.write(f'{columns[0]}\t{value}\n')

      # Read the next line
      if i+1 < len(wig_file_lines):
        # Split the next line into columns
        next_line = wig_file_lines[i+1]
        next_columns = next_line.split()
        # If the next line's first column is a number, check if it is consecutive with the current line's first column

        if next_columns[0].isdigit():
          # Get the next line's first column value
          next_value = int(next_columns[0])
          # If the values are not consecutive, insert "n" rows where n is the difference between the two values
          # and where the second column is equal to 1
          if (next_value - current_value) > 1:
            for j in range(current_value + 1, next_value):
              out_file.write(f'{j}\t{max_v}\n')

    # If the line is a metadata line, print it
    else:
      out_file.write(f'{line}')

    # Increment the counter and display progress
    line_counter += 1
    progress = line_counter / num_lines * 100
    if line_counter % int(num_lines/10) == 0:
      print(f'Progress: {progress:.2f}%', end='\r')
      # Print a newline after the progress message
      print()
  return out_file

# Takes in a iostring object and saves to a file
def write_wig(wig_file1,out_file_path):
  print("SAVING OUTPUT TO:", out_file_path)
  with open(out_file_path, 'w') as out_file:
    out_file.write(wig_file1.read())


def bam_mod_dict_to_wig(bam_file,output_file,bed_file):
  # Read in .bam file; If bam file is null or incorrect, return error message.
  # If bedfile, for each region in bedfile get each 
  # Write track line
  track_line = 'track type=wiggle_0 name="6mA" description="6mA" visibility=full color=0,0,0 altColor=0,0,0 priority=20\n'
  #


### Caculate average m6A per region based on input bed and bam files
def modbam_to_wig(bam_file, bed_file, threshold,output_stem):
  # Load the BED file
  regions = pysam.TabixFile(bed_file)
  max_regions = len(list(regions.fetch()))
  print("Loading BED file:", bed_file)
  print("BED file", bed_file, "has", max_regions, "regions")
  print("Loading BAM file:", bam_file)
  # Load the BAM file
  bam_ext = pysam.AlignmentFile(bam_file, "rb")
  print("BAM file", bam_file, "has", bam_ext.count(), "reads")

  # Initialize a list to store the results
  results_df = pd.DataFrame(columns=['chr', 'start', 'end', 'smt_pos', 'read_id'])

  # Iterate over the regions in the BED file
  # methylation_status = pd.DataFrame(columns=["chr","start","end","strand","region_type","chr_type",
  #                               "read","1","2","3","4","5","6","7","8","9","10","11","12"])
  region_count = 0
  for region in regions.fetch():
    region_count += 1
    # Split the region string into the chromosome, start, and end positions
    print("Starting on region:", region_count, " out of ", max_regions)
    chromosome, start, end, strand, region_type, chr_type = region.split()
    start = int(start)
    end = int(end)

    for read in bam_ext.fetch(chromosome, start, end):
      if not read.is_unmapped and read.is_forward == True:
        # Load up the modified base probilities
        m6A_dict = read.modified_bases[('A', 0, 'Y')]  # ('A', 0, 'a')

        # Convert list of tuples to dataframe in .wig format
        read_wig = pd.DataFrame(m6A_dict, columns=['variableStep', 'chrom=' + chromosome])
        read_wig['span=1'] = ""
        read_wig['variableStep'] = read_wig['variableStep'] + read.reference_start

        # Save temp .wig file
        read_wig.to_csv(output_stem + 'temp.wig', sep='\t', index=False)

        # open .wig as file in memory
        wig_file_obj = open(os.path.join(output_stem, "temp.wig"), "r")
        wig_file = io.StringIO(wig_file_obj.read())
        wig_file.seek(0)

        # Normalize wig file (divides by 256 so value is between 0-1)
        wig_file_norm = norm_wig(wig_file)
        # print("wig_file_norm",wig_file_norm.readlines())
        # wig_file_norm.seek(0)

        """# Invert wig file, since high methylation indicates no protein occupancy. Also puts on 0-10 scale.
        wig_file_inv = invert_wig(wig_file_norm, 10)
        # print("wig_file_inv",wig_file_inv.readlines())
        # wig_file_inv.seek(0)"""

        #set wigfile name to bamfilename but replacing .bam with .wig
        wig_filepath = os.path.join(output_stem, bam_file.replace(".bam", "_m6A_frac.wig"))

        # Output updated temp.wig file
        write_wig(wig_file_norm, wig_filepath)
  return wig_filepath

'''Hmerged_wig = merge_wig(wig_file1,wig_file2)
Dmerged_wig = merge_wig(wig_file3,wig_file4)

# Normalize last 3 wigs to first
wig_file1_norm = Hmerged_wig
wig_file2_norm = normalize_wigs(Hmerged_wig,Dmerged_wig)

# Output normalized files
#write_wig(wig_file1_norm,"m6A_fwd_norm.wig")
write_wig(wig_file1_norm,os.path.join("/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/wig/output/", "H_modified_bases.6mA.merged_strand_NOnorm.wig"))
write_wig(wig_file2_norm,os.path.join("/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/wig/output/", "D_modified_bases.6mA.merged_strand_Hnorm.wig"))
wig_file1_norm.seek(0)
wig_file2_norm.seek(0)

# Merge last 3 wigs into first
#merged_wig = merge_wig(wig_file1_norm,wig_file2_norm)
#merged_wig = merge_wig(merged_wig,wig_file3_norm)
#merged_wig = merge_wig(merged_wig,wig_file4_norm)

# Invert merged wig
inverted_wig1 = invert_wig(wig_file1_norm,10)
inverted_wig2 = invert_wig(wig_file2_norm,10)

# Fill wig
#python3 /Data1/software/DANPOS3/danpos.py dpos /Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/wig/output/D_modified_bases.6mA.merged_strand_Hnorm_10SCALE.wig -q 9 -o result_m6A_9
#python3 /Data1/software/DANPOS3/danpos.py dpos /Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/wig/output/H_modified_bases.6mA.merged_strand_NOnorm_10SCALE.wig -q 9 -o result_m6A_9

# Write wig
write_wig(inverted_wig1,os.path.join("/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/wig/output/", "H_modified_bases.6mA.merged_strand_NOnorm_10SCALE.wig"))
write_wig(inverted_wig2,os.path.join("/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/wig/output/", "D_modified_bases.6mA.merged_strand_Hnorm_10SCALE.wig"))
'''
# DPOS Commands

"""with open(wig_file_name, 'r') as wig_file:
  # Get the total number of lines in the file
  num_lines = sum(1 for line in wig_file)
  # Seek back to the beginning of the file
  wig_file.seek(0)
  wig_file_lines = wig_file.readlines()
  wig_file.seek(0)
  # Initialize a counter to track progress
  line_counter = 0

  with open(out_file_name, 'a') as out_file:
    # Iterate over each line in the file
    print("Converting value to 1-value:")
    for i in range(0,len(wig_file_lines)):
      line = wig_file_lines[i]
      # Split the line into columns
      columns = line.split()

      if columns[0].isdigit():
        # If the line is a data line, invert the value and print it
        current_value = int(columns[0])
        # Read the next line
        if i+1 < len(wig_file_lines):
          next_line = wig_file_lines[i+1]
          # Split the next line into columns
          next_columns = next_line.split()
          # If the next line's first column is a number, check if it is consecutive with the current line's first column
          value = round((1 - float(columns[1]))*10,2)
          out_file.write(f'{columns[0]}\t{value}\n')
          if next_columns[0].isdigit():
            # Get the next line's first column value
            next_value = int(next_columns[0])
            # If the values are not consecutive, insert "n" rows where n is the difference between the two values
            # and where the second column is equal to 1
            if (next_value - current_value) > 1:
              for j in range(current_value + 1, next_value):
                out_file.write(f'{j}\t10\n')

      # If the line is a data line, print it
      else:
        out_file.write(f'{line}')

      # Increment the counter and display progress
      line_counter += 1
      progress = line_counter / num_lines * 100
      if line_counter % int(num_lines/10) == 0:
        print(f'Progress: {progress:.2f}%', end='\r')
        # Print a newline after the progress message
        print()"""


# bash command to create new conda env called "DPOS" with python 3.7, rpy2 and all software dependencies
# conda create -n DPOS python=3.7 rpy2

# https://github.com/epi2me-labs/modbam2bed
# conda activate modbam2bed
# modbam2bed --extended -m 6mA -r CHROMOSOME_I:10000-100000 --threshold 0.5 /Data1/reference/c_elegans.WS235.genomic.fa ["/Data1/seq_data/AH_N2_SDC2aid_AuxRem_fiberseq_8_19_23/combined_pod5/barcode07/basecalls/barcode07.mod_mappings.sorted.bam","/Data1/seq_data/AH_N2_SDC2aid_AuxRem_fiberseq_8_19_23/combined_pod5/barcode08/basecalls/barcode08.mod_mappings.sorted.bam"]
#
# ["/Data1/seq_data/AH_N2_SDC2aid_AuxRem_fiberseq_8_19_23/pod5_pass/barcode07/basecalls/barcode07.mod_mappings.sorted.bam","/Data1/seq_data/AH_N2_SDC2aid_AuxRem_fiberseq_8_19_23/pod5_pass/barcode08/basecalls/barcode08.mod_mappings.sorted.bam"]

