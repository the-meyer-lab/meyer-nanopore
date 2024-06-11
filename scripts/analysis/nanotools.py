"""
nanotools.py

A collection of Python functions for analysis of nanopore datasets.

Created: July 28, 2023
Author: Yuri Malina (ymalina@berkeley.edu)
"""

# Standard library imports
import os
import random
import subprocess
from io import StringIO
import multiprocessing
from multiprocessing import Pool, Lock
import inspect
from typing import Optional

# Third-party imports
import numpy as np
import pandas as pd
from tqdm import tqdm

import plotly.graph_objects as go
import plotly.express as px # Used for plotting
import pyBigWig
import pysam
import seaborn as sns
from scipy import stats
from scipy.signal import find_peaks
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Set seaborn style for all plots
sns.set(style="whitegrid")

# Function to process each type and bam file
# Function to preprocess each type
def preprocess_type_for_enrichment(each_type, combined_bed_df):
    each_type_df = combined_bed_df[combined_bed_df["type"] == each_type]
    each_type_df.drop(each_type_df.columns[-2:], axis=1, inplace=True)

    each_type_df.insert(len(each_type_df.columns) - 1, len(each_type_df.columns) - 1, ".", allow_duplicates=True)
    each_type_df.insert(len(each_type_df.columns) - 1, len(each_type_df.columns) - 1, ".", allow_duplicates=True)

    each_type_df.columns = range(len(each_type_df.columns))
    return each_type_df

# Function to process each type and bam file
def process_type_bam_for_enrichment(each_type_df, each_type, each_bam, each_cond, each_thresh, modkit_path):
    temp_filename = f"temp_files/temp_{each_type}.bed"
    each_type_df.to_csv(temp_filename, sep="\t", header=False, index=False)

    summary_df = get_summary_from_bam(sampling_frac=0.1, a_threshold=each_thresh,
                                                modkit_path=modkit_path, bam_path=each_bam,
                                                each_condition=each_type, each_exp_id=each_cond,
                                                thread_ct=1000, bed_file=temp_filename)
    return summary_df


def plot_enrichment_results(results_df, numerator, denominator):
    """
    Plots the results of methylation analysis for 'a' and 'm' codes and calculates fold enrichment
    between two specified conditions.

    Parameters:
    - results_df: DataFrame containing the results of the analysis.
    - numerator: The condition to use as the numerator in the fold enrichment calculation.
    - denominator: The condition to use as the denominator in the fold enrichment calculation.
    """
    # Define the specific codes to plot
    codes_to_plot = ['a', 'm']

    for code in codes_to_plot:
        # Filter DataFrame for the current code
        filtered_df = results_df[results_df["code"] == code]

        # Plot 1: Average pass fraction by condition and experiment ID
        avg_pass_frac_df = filtered_df.groupby(["condition", "exp_id"])["pass_frac"].mean().reset_index()
        fig1 = px.bar(avg_pass_frac_df, x="exp_id", y="pass_frac", color="condition", barmode="group", text="pass_frac",
                      title=f"Average Pass Fraction for Code: {code.upper()}")
        fig1.update_traces(texttemplate='%{text:.1%}', textposition='outside')
        fig1.update_yaxes(title_text=f"m6A/{code.upper()}", tickformat=".0%", nticks=10)
        fig1.update_layout(width=900, height=600, template="plotly_white")
        fig1.show()

        # Enrichment Analysis for the specified numerator and denominator conditions
        numerator_df = filtered_df[filtered_df["condition"] == numerator].groupby(["exp_id"])[
            "pass_frac"].mean().reset_index()
        denominator_df = filtered_df[filtered_df["condition"] == denominator].groupby(["exp_id"])[
            "pass_frac"].mean().reset_index()

        # Merge the numerator and denominator DataFrames on exp_id
        enrichment_df = pd.merge(numerator_df, denominator_df, on="exp_id", suffixes=('_num', '_den'))
        enrichment_df["fold_enrichment"] = enrichment_df["pass_frac_num"] / enrichment_df["pass_frac_den"]

        # Display enrichment DataFrame for each code
        print(f"Fold Enrichment Data for Code: {code} (Numerator: {numerator}, Denominator: {denominator})")
        display(enrichment_df)

        # Plot 2: Fold enrichment as a bar plot for each code
        fig2 = px.bar(enrichment_df, y="fold_enrichment", x="exp_id", color="exp_id", text="fold_enrichment",
                      title=f"Fold Enrichment for Code: {code.upper()} (Numerator: {numerator}, Denominator: {denominator})")
        fig2.update_traces(texttemplate='%{text:.2f}X', textposition='outside')
        fig2.update_layout(width=900, height=600, template="plotly_white")
        fig2.show()


def create_lookup_bed(new_bed_files):
    """
    Converts standard bed files into a bed file formatted for modkit input.

    Parameters:
    - new_bed_files: List of paths to new bed files.

    Returns:
    - combined_bed_df: DataFrame containing combined bed data.
    """
    combined_bed_df = pd.DataFrame()
    for each_bed in new_bed_files:
        bed_path = each_bed[:-3]  # Remove .gz extension
        bed_df = pd.read_csv(bed_path, sep="\t", header=None)
        combined_bed_df = combined_bed_df.append(bed_df)

    combined_bed_df.columns = ['chrom', 'bed_start', 'bed_end', 'bed_strand', 'type', 'chr_type']
    combined_bed_df.sort_values(by=['chrom', 'bed_start'], ascending=True, inplace=True)
    combined_bed_df.dropna(inplace=True)
    combined_bed_df.drop_duplicates(inplace=True)
    combined_bed_df.reset_index(drop=True, inplace=True)

    return combined_bed_df

def filter_bed_file(bed_file, sample_source, selection, chromosome_selected, chr_type_selected, type_selected,
                    strand_selected, max_regions, bed_window, intergenic_window):
    """
    Filters a bed file based on various criteria and saves 1 bed file per sample source type.

    Parameters:
    - bed_file: Path to the complete bed file.
    - sample_source: Type of sample source to filter by ("type", "chr_type", "chromosome").
    - selection: List of types to select.
    - chromosome_selected: List of chromosomes to select.
    - chr_type_selected: List of chromosome types to select.
    - type_selected: List of types to select.
    - strand_selected: List of strands to select.
    - max_regions: Maximum number of regions to select per type.
    - bed_window: Number of basepairs to add to each side of the bed file.
    - intergenic_window: Window size for intergenic regions.

    Returns:
    - List of paths to saved bed files.
    """
    print("Filtering bed file...")
    print("Configs:",sample_source, selection, chromosome_selected, chr_type_selected, type_selected, strand_selected, max_regions, bed_window, intergenic_window)
    
    ### Select bed file
    full_bed = pd.read_csv(bed_file,sep='\t')
    bed=[]

    # Select rows where type == "whole_chr" and store end value to chromosome_ends
    chromosome_ends = full_bed[full_bed["type"] == "whole_chr"]
    print("Chromosome ends:",chromosome_ends)

    for each_type in selection:
    # REGION CONFIGURATION
        if sample_source == "type":
            temp_bed = full_bed[full_bed["chromosome"].isin(chromosome_selected) &
                                full_bed["chr-type"].isin(chr_type_selected) &
                                #full_bed["type"] is the same as each_type
                                full_bed["type"].__eq__(each_type) &
                                full_bed["strand"].isin(strand_selected)]
        if sample_source == "chr_type":
            temp_bed = full_bed[full_bed["chromosome"].isin(chromosome_selected) &
                                full_bed["chr-type"].str.contains(each_type) &
                                full_bed["type"].isin(type_selected) &
                                full_bed["strand"].isin(strand_selected)]
        if sample_source == "chromosome":
            temp_bed = full_bed[full_bed["chromosome"].__eq__(each_type) &
                                full_bed["chr-type"].isin(chr_type_selected) &
                                full_bed["type"].isin(type_selected) &
                                full_bed["strand"].isin(strand_selected)]
        
        # Drop random regions to match max_regions
        drop_count = len(temp_bed)-max_regions
        # If max regions > selected regions, do not drop any.
        if(drop_count<0):
            drop_count=0
        # If max_regions = 0, do not drop any.
        if (max_regions == 0):
            drop_count = 0
        # Drop random regions to match max_regions
        temp_bed = temp_bed.copy()
        drop_indices = np.random.choice(temp_bed.index, drop_count, replace=False)
        temp_bed.drop(drop_indices,inplace=True)

        # Sort by chromosome and start
        temp_bed.sort_values(by=["chromosome","start"],ascending=True,inplace=True)
        temp_bed.reset_index(drop=True, inplace=True)

        # Adjust start for rows where type != "intergenic_control" by bed_window, otherwise ajust by intergenic_window
        temp_bed["start"] = np.where(temp_bed["type"] != "intergenic_control", temp_bed["start"] - bed_window, temp_bed["start"] - intergenic_window)
        # Adjust end for rows where type != "intergenic_control" by bed_window, otherwise ajust by intergenic_window
        temp_bed["end"] = np.where(temp_bed["type"] != "intergenic_control", temp_bed["end"] + bed_window, temp_bed["end"] + intergenic_window)

        # Set start for all rows where start < 0 to 0
        temp_bed["start"] = np.where(temp_bed["start"] < 0, 0, temp_bed["start"])
        # For each chromosome, set end for all rows where end > chromosome_ends to chromosome_ends
        for each_chromosome in chromosome_ends["chromosome"].unique():
            chromosome_end = chromosome_ends[chromosome_ends["chromosome"].__eq__(each_chromosome)]["end"].values[0]
            temp_bed["end"] = np.where((temp_bed["end"] > chromosome_end) & (temp_bed["chromosome"] == each_chromosome), chromosome_end, temp_bed["end"])


        # If select_opp_strand is True, copy dataframe, reverse strand column and append to original df
        #if select_opp_strand == True:
        #    # set temp_bed_opp_min to copy of temp_bed, but only where strand == "-"
        #    temp_bed_opp_min = temp_bed[temp_bed["strand"].str.contains("-")]
        #    temp_bed_opp_min["strand"] = temp_bed_opp_min["strand"].replace({"-":"+"})

            # copy only strand "+" from temp_bed_opp and reverse strand
        #    temp_bed_opp_plus = temp_bed[temp_bed["strand"].str.contains("\+")]
        #    temp_bed_opp_plus["strand"] = temp_bed_opp_plus["strand"].replace({"+":"-"})
        #    temp_bed = pd.concat([temp_bed,temp_bed_opp_min,temp_bed_opp_plus],ignore_index=True)

        #    temp_bed.sort_values(by=["chromosome","start"],ascending=True,inplace=True)
        #    temp_bed.reset_index(drop=True, inplace=True)

        #display_sample_rows(temp_bed, 5)
        # select only path from bed_file
        bed_file_path = os.path.dirname(bed_file)

        # Save bed file
        temp_bedfile = bed_file_path+"/"+each_type+".bed"
        temp_bedfile_gz = bed_file_path + "/"+each_type+".bed.gz"
        temp_bedfile_tbi = bed_file_path + "/" + each_type + ".bed.gz.tbi"
        temp_bed.to_csv(temp_bedfile, sep="\t",header=False,index=False)

        # Compress bed file
        with open(temp_bedfile_gz, 'wb') as output_file:
            result = subprocess.run(["bgzip", "-c", temp_bedfile],stdout=output_file, stderr=subprocess.PIPE)
        if result.returncode != 0:
            print(f"Error executing command: {result.stderr.decode('utf-8')}")
        else:
            print(f"{temp_bedfile} has been compressed successfully to {temp_bedfile_gz}")

        # Index bed file
        with open(temp_bedfile_tbi, 'wb') as output_file:
            result2 = subprocess.run(["tabix", "-f", "-p", "bed", temp_bedfile_gz],stdout=output_file,stderr=subprocess.PIPE)
        if result2.returncode != 0:
            print(f"Error executing command: {result2.stderr.decode('utf-8')}")
        else:
            print(f"Index created successfully for {temp_bedfile_tbi}")

        # For first iteration
        if bed == []:
            bed = [temp_bedfile_gz]

        # Otherwise append region to temporary bed file.
        else:
            bed.append(temp_bedfile_gz)
    print("Saved the following bedfiles:",bed)
    return bed


### Extract only reads in bam file that overlap selected regions, and subselect down using fraction
### NOT NECESSARY IF RUNNING WHOLE CHROMOSOMES
# Define function to generate redux .bam file

def subsample_bam(bam_file, condition, bam_frac, sample_index, output_stem):
    """
    Subsamples a .bam file based on the specified fraction.

    Parameters:
    - bam_file: Path to the .bam file.
    - condition: Condition name.
    - bam_frac: Fraction of reads to keep.
    - sample_index: Sample index.
    - output_stem: Path to the output .bam file.

    Returns:
    - Output path of the subsampled .bam file.
    """
    if bam_frac == 1:
        # Second command: Index the output BAM file
        # Index bam_file if bam_file + ".bai" does not exist
        if not os.path.exists(bam_file + ".bai"):
            result = subprocess.run(['samtools', 'index', bam_file], stderr=subprocess.PIPE)
            if result.returncode == 0:
                print(f"BAM file indexed successfully: {bam_file}.bai")
            else:
                print(f"Error indexing BAM file: {result.stderr.decode('utf-8')}")
        return bam_file
    else:
        lock1 = multiprocessing.Lock()
        bam_list = []
        with lock1:
            # find root of output_stem path
            output_bamfile = os.path.dirname(bam_file) + "/mod_mappings_" + condition + str(sample_index)+"_" + str(bam_frac) + ".frac_sorted.bam"
            # Add output bam file to list
            # If output_bamfile exists, do not run command
            if os.path.exists(output_bamfile):
                print(f"Output file already exists: {output_bamfile}")

            else:
                print("starting on: ", bam_file)
                with open(output_bamfile, 'wb') as output_file:
                    process1 = subprocess.Popen(['samtools', 'view', '-h', '-s', str(bam_frac), bam_file],
                                                stdout=subprocess.PIPE)
                    process2 = subprocess.Popen(['samtools', 'view', '-h', '-b', '-'], stdin=process1.stdout,
                                                stdout=output_file)
                    process1.stdout.close()
                    process2.communicate()

                if process2.returncode != 0:
                    print("Error executing command")
                else:
                    print(f"Command executed successfully. Output written to {output_bamfile}")

                # Second command: Index the output BAM file
                result = subprocess.run(['samtools', 'index', output_bamfile], stderr=subprocess.PIPE)
                if result.returncode == 0:
                    print(f"BAM file indexed successfully: {output_bamfile}.bai")
                else:
                    print(f"Error indexing BAM file: {result.stderr.decode('utf-8')}")
            # return bam_list
            return output_bamfile


def parallel_subsample_bam(bam_files, conditions, bam_fracs, sample_indices, output_stems):
    """
    Parallel processing wrapper for subsampling BAM files. Constructs args_list internally.

    Parameters:
    - bam_files: List of paths to the .bam files.
    - conditions: List of condition names.
    - bam_fracs: List of fractions of reads to keep.
    - sample_indices: List of sample indices.
    - output_stems: List of paths to the output .bam files.
    """
    # Construct args_list here
    args_list = [(bam_file, condition, bam_frac, sample_index, output_stem)
                 for bam_file, condition, bam_frac, sample_index, output_stem
                 in zip(bam_files, conditions, bam_fracs, sample_indices, output_stems)]

    new_bam_files = []
    with multiprocessing.Pool() as pool:  # Adjust the number of processes if necessary
        new_bam_files = pool.starmap(subsample_bam, args_list)
    return new_bam_files

def extract_m6A_per_region(bam_file, bed_file, threshold, condition):
    """
    Extracts the number of m6A sites per region from a .bam file and a .bed file based on a threshold.

    Parameters:
    - bam_file: Path to the .bam file.
    - bed_file: Path to the .bed file.
    - threshold: m6A threshold for inclusion.
    - condition: Condition name for labeling output.

    Returns:
    - DataFrame: Contains the number of m6A sites and additional region information.
    """
    # Initialize a counter for the number of regions with at least one m6A site
    min_m6a_for_incl = 0
    null_counter = 0
    insuf_m6a_counter = 0
    global bed_window

    # Load the BAM file
    bam_ext = pysam.AlignmentFile(bam_file, "rb")

    # Load the BED file
    #regions = pysam.BedFile(bed_file)
    regions = pysam.TabixFile(bed_file)

    # Initialize a list to store the results
    results = []

    # Iterate over the regions in the BED file
    for region in regions.fetch(multiple_iterators=True):
        # Split the region string into the chromosome, start, and end positions
        chromosome, start, end, strand, region_type, chr_type = region.split()
        start = int(start)
        end = int(end)

        # Initialize counters for the total number of bases and the total number of m6A
        total_bases = 0
        total_m6A = 0
        read_counter = 0

        # Iterate over the reads that overlap the region
        # Note the need for multiple iterators: https://pysam.readthedocs.io/en/latest/faq.html?highlight=.fetch#pysam-coordinates-are-wrong
        for read in bam_ext.fetch(chromosome, start, end, multiple_iterators=True):
            # if read seq is not null:
            if (read.query_alignment_sequence != None):
                # print first 10 characters and length of read sequence
                # print("Read name:",read.query_name," | Read seq [0:10]",read.query_alignment_sequence[:10]," | length:",len(read.query_alignment_sequence))
                # print("starting on read:", read)

                # Count the total number of "A" bases in the read that overlap the region
                if read.is_forward == True:
                    # print(read.query_alignment_sequence)
                    if ('A', 0, 'Y') in read.modified_bases.keys() or ('A', 0, 'a') in read.modified_bases.keys():
                        try:
                            # try to get the modified base values from the forward strand with a capital Y
                            m6A_dict_values = [x[1] for x in read.modified_bases_forward[('A', 0, 'Y')]]
                        except:
                            # if that doesn't work, try to get the modified base values from the forward strand with a lowercase a
                            m6A_dict_values = [x[1] for x in read.modified_bases_forward[('A', 0, 'a')]]

                        read_bases = read.query_alignment_sequence.count("A")
                        read_mod_bases = sum(i > threshold for i in m6A_dict_values)

                        if read_bases > 0:
                            if read_mod_bases/read_bases > min_m6a_for_incl:
                                total_bases = total_bases + read_bases
                                total_m6A = read_mod_bases + total_m6A
                                read_counter = read_counter + 1
                            else:
                                insuf_m6a_counter = insuf_m6a_counter + 1
                        else:
                            insuf_m6a_counter = insuf_m6a_counter + 1

                    else:
                        # Unable to find modified base dictionary in read
                        print("No mod A dict!")
                        print(read)
                        print(read.modified_bases.keys())
                        print(read.modified_bases)
                else:
                    # BAM file stores the reverse complement of the actual read sequence for reverse strand reads.
                    # Therefore we need to count Ts for reverse strands
                    if ('A', 1, 'Y') in read.modified_bases.keys() or ('A', 1, 'a') in read.modified_bases.keys():
                        try:
                            m6A_dict_values = [x[1] for x in read.modified_bases[('A', 1, 'Y')]]
                        except:
                            m6A_dict_values = [x[1] for x in read.modified_bases[('A', 1, 'a')]]
                        read_bases = read.query_alignment_sequence.count("T")
                        read_mod_bases = sum(i > threshold for i in m6A_dict_values)

                        if read_bases > 0:
                            if read_mod_bases/read_bases > min_m6a_for_incl:
                                total_bases = total_bases + read_bases
                                total_m6A = read_mod_bases + total_m6A
                                read_counter = read_counter + 1
                            else:
                                insuf_m6a_counter = insuf_m6a_counter + 1
                        else:
                            insuf_m6a_counter = insuf_m6a_counter + 1
                    else:
                        # Unable to find modified base dictionary in read
                        print("No mod A dict!")
                        print(read)
                        print(read.modified_bases.keys())
                        print(read.modified_bases)

                # Reset the m6A_dict_values list for the next read.
                m6A_dict_values = []
            else:
                # If read has no sequence
                # increment null counter and print total
                null_counter = null_counter + 1
                if (null_counter % 1000 == 0):
                    print("Null counter:", null_counter)


        # Add the region information to the results list
        results.append([chromosome, start, end, region_type, chr_type, total_bases, total_m6A, read_counter])
        # print("Appending:",results)

    print("Insufficient m6A reads = ", insuf_m6a_counter)
    # Close the BAM file
    bam_ext.close()
    del bam_ext

    # Convert the results list to a pandas dataframe
    df = pd.DataFrame(results,
                      columns=["chromosome", "start", "end", "region_type", "chr_type", "total_bases", "total_m6A",
                               "overlapping_reads"])
    df["m6A_frac"] = df["total_m6A"] / df["total_bases"]
    df["norm_m6A_frac"] = df["m6A_frac"]/((df['total_m6A'].sum()/df['total_bases'].sum()))
    df["condition"] = condition
    print(condition," results df:",df.head())
    # return the dataframe
    return df

def extract_m6A_per_region_parallelized(bam_file, condition, bam_frac, file_prefix, selection, m6A_thresh, output_stem, new_bed_files):
    """
    Parallelized extraction of m6A sites per region from a .bam file across multiple bed files.

    Parameters:
    - bam_file: Path to the .bam file.
    - condition: Condition name.
    - bam_frac: Fraction of reads to keep.
    - file_prefix: Prefix for output files.
    - selection: List of types to select.
    - m6A_thresh: m6A threshold for inclusion.
    - output_stem: Path to output directory.
    - new_bed_files: List of bed file paths.

    Returns:
    - DataFrame: Combined data from all bed files.
    """
    combined_df = pd.DataFrame()
    for each_bedfile, each_type in zip(new_bed_files,selection):
        lock1 = multiprocessing.Lock()
        with lock1:
            print("starting with m6A_thresh =", m6A_thresh, "on ", each_type, " :", bam_file," with bedfile:",each_bedfile)
        output_df = extract_m6A_per_region(bam_file, each_bedfile, m6A_thresh, condition)
        output_df.to_csv(
            output_stem + file_prefix + "m6A_frac_" + condition + "_" + str(m6A_thresh) + "_" + each_type + ".csv",
            index=False, mode='w')
        combined_df = pd.concat([combined_df, output_df])
    print(combined_df)
    return combined_df

def frequency_table(df, n_dist=None):
    """
    Generates a frequency table of nucleosome repeat lengths from a DataFrame.

    Parameters:
    - df: DataFrame containing nucleosome positions.
    - n_dist: Optional; number of distances to consider.

    Returns:
    - DataFrame: Frequency table of nucleosome repeat lengths.
    """
    # Default n_dist to 10 if not provided
    if n_dist == None:
        n_dist = 10
    frequency_df = pd.DataFrame(columns=["chr", "dist", "n-plus", 'smt_pos'])
    column_labels = []

    for i in range(1, n_dist + 1):
        # create list of column labels and assign them to a new df, based on the number of n+ to consider.
        column_label = (str(i))
        column_labels.append(column_label)
        df = df.assign(**{column_label: df['smt_pos'].diff(periods=i)})

    # Create an empty list to store the tuples for the new dataframe
    data = []

    # Iterate through the columns of the original dataframe
    for col in df[column_labels]:
        # Get the values for the current column
        values = df[col].values
        chr_values = df["chr"].values
        smt_pos = df["smt_pos"].values
        # Create a tuple for each value with the column name as the first element and the value as the second element
        for chr_val, val, smt_pos in zip(chr_values, values, smt_pos):
            data.append((chr_val, col, val, smt_pos))

    # Create the new dataframe using the list of tuples
    frequency_df = pd.DataFrame(data, columns=['chr', 'n-plus', 'dist', 'smt_pos'])
    frequency_df = frequency_df.dropna()
    frequency_df = frequency_df[frequency_df['dist'] > 0]

    for i in range(1, 100):
        max_val = frequency_df['dist'].max()
        frequency_df = frequency_df[frequency_df['dist'] != max_val]

    frequency_df['n-plus'] = frequency_df['n-plus'].astype(int)
    # frequency_df["nrl"] = frequency_df["dist"]/frequency_df["n-plus"]

    # get the values of the 'pos' column as a numpy array
    return frequency_df

def plot_NRL_dist(data, chromosome_name, color, output_prefix, smoothing_val=None):
    """
    Plots the distribution of nucleosome-nucleosome distances (NRL).

    Parameters:
    - data: DataFrame containing distance data.
    - chromosome_name: Name of the chromosome.
    - color: Color for the plot.
    - output_prefix: Prefix for saving the plot.
    - smoothing_val: Optional; smoothing value for KDE.

    Returns:
    - peaks: DataFrame containing peak positions.
    - fig: The plot figure object.
    """
    # Clear .plt
    #plt.clf()
    sns.set_style("white")
    fig = sns.displot(data, x="dist", color=color,label=chromosome_name, kind="kde", height=6, aspect=1.5, linewidth=3,bw_adjust=smoothing_val,clip=(-1000,1000))
    plt.title(f'{output_prefix}\n{chromosome_name} Nucleosome-Nucleosome Distance', fontsize=15)
    ax1 = fig.facet_axis(0, 0)
    #ax1.set_xlim(0, 1200)
    #ax1.set_ylim(0, 0.002)

    line_x = ax1.lines[0].get_xdata() # Get the x data of the distribution
    line_y = ax1.lines[0].get_ydata() # Get the y data of the distribution
    peaks_index = find_peaks(line_y)
    peak_number = np.array(range(1, len(peaks_index[0]) + 1))
    peaks_x = [line_x[i] for i in peaks_index[0]]
    peaks = pd.DataFrame({'x': peak_number, f'{chromosome_name}-peak-position': peaks_x})
    #peaks.drop(peaks.tail(len(peak_number)-5).index, inplace=True)
    # pos_peaks = number of peaks > 0
    neg_peaks = len(peaks[peaks[f'{chromosome_name}-peak-position'] < 0])
    inc = -neg_peaks
    for i, _x in enumerate(peaks[f'{chromosome_name}-peak-position']):
        # add a vertical line for the value
        ax1.axvline(x=_x, lw=1, ls=":", color='grey')
        # add a label below the line
        ax1.text(_x, 0.00001, str(int(_x)), ha='center', fontsize=10)
        if inc == 0:
            inc+= 1
            #ax1.text(_x, 0.0002, "N", ha='center', fontsize=10)
        if inc >0:
            ax1.text(_x, 0.0002, "N+" + str(inc), ha='center', fontsize=10)
        else:
            ax1.text(_x, 0.0002, "N" + str(inc), ha='center', fontsize=10)
        inc += 1
    ax1.legend(loc="best")
    plt.xlabel('nucleosome-nucleosome distance (bp)', fontsize=12)
    fig.savefig(f'images/{output_prefix}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'images/{output_prefix}.svg', bbox_inches='tight')

    plt.show()
    # return fig and peaks
    return peaks,fig

def plot_NRL_dist_compare(data, chromosome_name, output_prefix, smoothing_val=None, norm_bool=True, hue='condition', window=1000):
    """
    Compares NRL distributions across conditions.

    Parameters:
    - data: DataFrame with NRL data.
    - chromosome_name: Chromosome name for labeling.
    - output_prefix: Output prefix for saved plots.
    - smoothing_val: Optional smoothing value for KDE plots.
    - norm_bool: Normalize distributions if True.
    - hue: DataFrame column name to color-code data.
    - window: Range window for plotting.

    Returns:
    - plt: Plot object.
    """
    sns.set_style("white")
    fig = sns.displot(data, x="dist", hue = hue,label=chromosome_name, kind="kde", height=6, aspect=1.5, linewidth=3,common_norm=norm_bool,clip=(-window,window),bw_adjust=smoothing_val)
    plt.title(f'{output_prefix}\n{chromosome_name} Nucleosome-Nucleosome Distance', fontsize=15)
    ax1 = fig.facet_axis(0, 0)
    #ax1.set_xlim(0, 1200)
    #ax1.set_ylim(0.00045, )
    #ax1.set_ylim(0, 0.002)

    ax1.legend(loc="best")
    plt.xlabel('nucleosome-nucleosome distance (bp)', fontsize=12)
    # set theme 
    #fig.savefig(f'images/{output_prefix}.png', dpi=300, bbox_inches='tight')
    #fig.savefig(f'images/{output_prefix}.svg', bbox_inches='tight')
    #plt.show()
    # return fig and peaks
    return plt

def plot_NRL_regression(dataframe, title, output_prefix, equation_table):
    """
    Plots a regression analysis of nucleosome-nucleosome distances.

    Parameters:
    - dataframe: DataFrame containing the data.
    - title: Title for the plot.
    - output_prefix: Prefix for output files.
    - equation_table: Table of regression equations.

    Displays the plot and saves output files.
    """
    # set default sns background to white
    sns.set_style("white")
    # Clear any previous plots
    plt.clf()

    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))

    # fetch text that preceeds '-n-plus' in column names in dataframe, ignoring all other column names
    unique_chromosomes = set([col.split('-n-plus')[0] for col in dataframe.columns if '-n-plus' in col])
    
    num_chrom = len(unique_chromosomes)

    # Create a dictionary "chromosomes" that assigns a single color to each element in unique_chromosomes
    chromosomes= dict(zip(unique_chromosomes, sns.color_palette("husl", num_chrom)))

    for chrom, color in chromosomes.items():
        sns.regplot(x=f'{chrom}-n-plus', y=f'{chrom}-peak-position', data=dataframe, ax=ax, label=chrom, color=color,ci=None)

    # remove transparency from plot

    ax.legend(loc="best")
    ax.set_title(title)
    ax.set_ylabel('nucleosome-nucleosome distance (bp)', fontsize=12)
    ax.set_xlabel('peak number (N+)', fontsize=12)
    # Add each element in equations table to lower right hand of the plot. Making sure that the last element is at the bottom.
    for i, equation in enumerate(equation_table[::-1]):
        ax.text(0.95, 0.05 + i * 0.05, equation, horizontalalignment='right', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)

    fig.savefig(f'jupyter_output/{output_prefix}_NRL_trend.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'jupyter_output/{output_prefix}_NRL_trend.svg', bbox_inches='tight')
    # set plt background to white
    plt.style.use('default')
    plt.show()

### Define cuntion to nucleosome positions by genomic features
def filter_nucs_by_features(data_df, bed_tss, region_cutoff):
    """
    Filters nucleosome positions by proximity to features in a bed file.

    Parameters:
    - data_df: DataFrame with nucleosome positions.
    - bed_tss: DataFrame with TSS positions from a bed file.
    - region_cutoff: Distance cutoff for filtering.

    Returns:
    - DataFrame: Filtered nucleosome positions.
    """
    all_results = []

    for chrm in data_df.chr.unique():
        bed_tss_chrm = bed_tss.loc[bed_tss['chromosome'] == chrm]
        bed_tss_chrm = bed_tss_chrm.sort_values(by="start")
        data_df_chrm = data_df.loc[data_df['chr'] == chrm]
        data_df_chrm = data_df_chrm.sort_values(by="smt_pos")

        # merge_asof
        data_df_chrm = pd.merge_asof(
            data_df_chrm,
            bed_tss_chrm,
            left_on='smt_pos',
            right_on='start'
        )

        # Filtering by the region cutoff
        data_df_chrm['sub'] = abs(data_df_chrm['smt_pos'] - data_df_chrm['start'])
        data_df_chrm = data_df_chrm.loc[data_df_chrm['sub'] < region_cutoff]
        data_df_chrm = data_df_chrm.iloc[:, :-len(bed_tss.columns)-1]

        all_results.append(data_df_chrm)

    result_df = pd.concat(all_results, ignore_index=True)

    return result_df

### Calculate bam summary statistics
def get_summary_from_bam(sampling_frac, a_threshold, modkit_path, bam_path, each_condition, each_exp_id, thread_ct=4, chromosome=None, start=None, end=None, bed_file=None):
    """
    Fetches summary statistics from a BAM file using modkit.

    Parameters:
    - sampling_frac: Fraction of reads to sample.
    - a_threshold: Threshold for 'a' modifications.
    - modkit_path: Path to the modkit executable.
    - bam_path: Path to the BAM file.
    - each_condition: Condition label.
    - each_exp_id: Experiment ID.
    - thread_ct: Number of threads to use.
    - chromosome: Optional; specific chromosome to analyze.
    - start: Optional; start position for analysis.
    - end: Optional; end position for analysis.
    - bed_file: Optional; path to a bed file for targeted analysis.

    Returns:
    - DataFrame: Summary statistics.
    """

    #print("Starting on ", each_condition, " with exp_id: ", each_exp_id)

    command = [
        modkit_path,
        "summary",
        #"--ignore",
        #"m",
        "--threads",
        f"{thread_ct}",
        "--sampling-frac",
        f"{sampling_frac}",
        "--seed",
        "10",
        #"--mod-thresholds",
        #f"a:{a_threshold}",
        #"--no-filtering",
        "--filter-threshold",
        f"A:{1-a_threshold}",
        "--filter-threshold",
        f"C:{1 - a_threshold}",
        "--mod-thresholds",
        f"a:{a_threshold}",
        "--mod-thresholds",
        f"m:{a_threshold}",
        "--log-filepath",
        f"temp_files/modkit_summary_{each_condition}_{each_exp_id}.log"
    ]

    # if bed_file is not None, add the bed file argument to the end of the command
    if bed_file is not None:
        command.extend(["--include-bed",
                        f"{bed_file}", "--only-mapped"])

    # if chromosome, start and end are not None, add the region argument
    elif chromosome is not None and start is not None and end is not None:
        region_argument = f"{chromosome}:{start}-{end}"
        command.extend(["--region", region_argument])

    command.extend([bam_path])

    # Run the command and capture the stdout
    #print("Running command:", " ".join(command))
    result = subprocess.run(command,  text=True, capture_output=True)
    if chromosome is None:
        print('stdout:')
        print(result.stdout)

    # Extract only the table data starting from the column names
    # Split the stdout using the substring
    split_data = result.stdout.split(" base  ")

    # Check if we have at least two items after the split
    if len(split_data) > 1:
        table_data = "base " + split_data[1]  # Adding back the column name
    else:
        # Here you can decide on what to do when the substring doesn't exist
        # For this example, I'll return an empty DataFrame
        print("Warning: Expected substring 'base' not found in the output.")
        return pd.DataFrame(columns=['condition', 'exp_id'])

    # Convert the extracted data to a DataFrame
    df = pd.read_csv(StringIO(table_data), sep="\s+", engine='python')

    # Add condition column
    df["condition"] = each_condition
    df["exp_id"] = each_exp_id

    return df

def display_sample_rows(df, n=5):
    """
    Displays sample rows from a DataFrame.

    Parameters:
    - df: DataFrame to display samples from.
    - n: Number of rows to display from the head, random, and tail of the DataFrame.

    No return value.
    """
    frame = inspect.currentframe().f_back
    varname = [k for k, v in frame.f_locals.items() if v is df][0]

    # if length of df is less than n then display all rows
    if len(df) < n:
        n = len(df)
        print(f"| {varname} | first {n} out of total {len(df)} rows.")
        sampled_df = df.head(n)
    else:
        print(f"| {varname} | first {n}, random {n} and last {n} out of total {len(df)} rows.")
        sampled_df = pd.concat([df.head(n), df.sample(n=n, random_state=1), df.tail(n)], ignore_index=True)
    display(sampled_df)

def generate_modkit_bed(new_bed_files, down_sample_autosome, select_opp_strand, output_name):
    """
    Generates a bed file formatted for modkit input from a list of bed files,
    with options for downsampling and strand selection.

    Parameters:
    - new_bed_files: List of bed file paths to combine.
    - down_sample_autosome: Boolean indicating whether to downsample autosomal regions.
    - select_opp_strand: Boolean indicating whether to select opposite strands for regions.
    - output_name: Path for the output bed file.

    Returns:
    - DataFrame: Combined bed data ready for modkit input.
    """
    # Initialize an empty DataFrame to store the combined data from all bed files
    combined_bed_df = pd.DataFrame()

    # Read each bed file and append it to the combined DataFrame
    for each_bed in new_bed_files:
        bed_path = each_bed[:-3]
        temp_df = pd.read_csv(bed_path, sep="\t", header=None)
        """# Downsample autosome genes if specified
        if down_sample_autosome and temp_df[0].str.contains("X").sum() < len(temp_df) - temp_df[0].str.contains("X").sum():
            x_genes = temp_df[temp_df[0].str.contains("X")]
            autosome_genes = temp_df[~temp_df[0].str.contains("X")].sample(n=len(x_genes),random_state=1)
            temp_df = pd.concat([x_genes, autosome_genes], ignore_index=True)"""
        combined_bed_df = combined_bed_df.append(temp_df)

    # Drop the last column
    combined_bed_df.drop(combined_bed_df.columns[len(combined_bed_df.columns)-1], axis=1, inplace=True)
    combined_bed_df.drop(combined_bed_df.columns[len(combined_bed_df.columns)-1], axis=1, inplace=True)

    # Insert two "." columns before the last column
    combined_bed_df.insert(len(combined_bed_df.columns) - 1, len(combined_bed_df.columns) - 1, ".", allow_duplicates=True)
    combined_bed_df.insert(len(combined_bed_df.columns) - 1, len(combined_bed_df.columns) - 1, ".", allow_duplicates=True)


    # Reset the column labels
    combined_bed_df.columns = range(len(combined_bed_df.columns))

    # Downsample autosome genes if specified
    if down_sample_autosome and combined_bed_df[0].str.contains("X").sum() < len(combined_bed_df) - combined_bed_df[0].str.contains("X").sum():
        x_genes = combined_bed_df[combined_bed_df[0].str.contains("X")]
        autosome_genes = combined_bed_df[~combined_bed_df[0].str.contains("X")].sample(n=len(x_genes), random_state=1)
        combined_bed_df = pd.concat([x_genes, autosome_genes], ignore_index=True)

    # Select opposite strands if specified
    if select_opp_strand:
        minus_strand_df = combined_bed_df[combined_bed_df[5] == "-"].copy()
        minus_strand_df[5] = "+"

        plus_strand_df = combined_bed_df[combined_bed_df[5] == "+"].copy()
        plus_strand_df[5] = "-"

        combined_bed_df = pd.concat([combined_bed_df, minus_strand_df, plus_strand_df], ignore_index=True)

    # Sort DataFrame and reset index
    combined_bed_df.sort_values([0, 1], inplace=True)
    combined_bed_df.reset_index(drop=True, inplace=True)

    # Save DataFrame to file
    combined_bed_df.to_csv(output_name, sep="\t", header=False, index=False)

    # Renaming combined_bed_df to modkit_bed_df to fulfill your requirement
    modkit_bed_df = combined_bed_df

    return modkit_bed_df

def extract_n50_from_bam(bam_file, fraction_to_sample=0.01):
    """
    Extracts sample read lengths from a BAM file to calculate N50 statistics.

    Parameters:
    - bam_file: Path to the BAM file.
    - fraction_to_sample: Fraction of reads to sample for N50 calculation.

    Returns:
    - DataFrame: Sampled read lengths.
    """
    # Get the total number of reads in the BAM file
    total_reads = int(subprocess.getoutput(f"samtools view -c {bam_file}"))

    # Calculate the fraction needed to get at least 10,000 reads
    min_reads = 5000
    min_fraction = min_reads / total_reads
    final_fraction = max(fraction_to_sample, min_fraction)

    # Seed can be any integer (here it's 42), and final_fraction is the fraction of reads you want to sample
    cmd = ["samtools", "view", "-s", f"42.{int(final_fraction * 100)}", bam_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"Samtools error: {result.stderr}")
        return None

    output = result.stdout
    read_data = []

    for line in StringIO(output):
        fields = line.strip().split("\t")
        query_name = fields[0]
        flag = int(fields[1])
        query_length = len(fields[9])

        read_data.append([query_name, query_length])

    df = pd.DataFrame(read_data, columns=['query_name', 'query_length'])

    return df

def process_bam_file_for_n50(args):
    """
    Processes a BAM file to extract read lengths for N50 calculation,
    designed to be used with multiprocessing.

    Parameters:
    - args: Tuple containing (bam_file, condition, exp_id).

    Returns:
    - DataFrame: Sampled read lengths with condition and experiment ID.
    """
    bam_file, condition, exp_id = args
    print(f"Starting on: {bam_file}, {condition}")
    temp_df = extract_n50_from_bam(bam_file, 0.1)
    temp_df['condition'] = condition
    temp_df['exp_id'] = exp_id
    return temp_df

def calculate_n50(lengths):
    """
    Calculates the N50 value from a list of lengths.

    Parameters:
    - lengths: List of sequence lengths.

    Returns:
    - int: The N50 value.
    """
    sorted_lengths = sorted(lengths, reverse=True)
    half_total_length = sum(sorted_lengths) / 2.0
    current_length = 0
    for length in sorted_lengths:
        current_length += length
        if current_length >= half_total_length:
            return length

def calculate_and_plot_n50(new_bam_files, conditions, exp_ids):
    """
    Calculates and plots N50 statistics for a set of BAM files.

    Parameters:
    - new_bam_files: List of paths to BAM files.
    - conditions: List of condition names corresponding to BAM files.
    - exp_ids: List of experiment IDs corresponding to BAM files.

    Returns:
    - Plotly graph object: Bar chart of N50 values by condition.
    """

    # Create a list of arguments
    args_list = [(bam_file, condition, exp_id) for bam_file, condition, exp_id in zip(new_bam_files, conditions, exp_ids)]

    # Initialize an empty DataFrame to hold the results
    combined_n50_df = pd.DataFrame()

    # Use a process pool to process the BAM files in parallel
    with Pool() as pool:  # Using all available CPUs
        results = pool.map(process_bam_file_for_n50, args_list)

    # Append the results to the combined DataFrame
    for temp_df in results:
        combined_n50_df = combined_n50_df.append(temp_df, ignore_index=True)

    # Reset the index
    combined_n50_df.reset_index(drop=True, inplace=True)

    # Calculate N50 for each 'condition' and create a new DataFrame
    n50_data = []
    for name, group in combined_n50_df.groupby('condition'):
        n50_value = calculate_n50(group['query_length'])
        n50_data.append([n50_value, name])

    n50_df = pd.DataFrame(n50_data, columns=['N50', 'condition'])

    # Create the bar chart using Plotly
    fig = px.bar(n50_df, x='condition', y='N50', color='condition', title='N50 by Condition', template='plotly_white', text='N50')
    fig.update_layout(autosize=False, width=400, height=400)
    return fig

#Source: https://stackoverflow.com/questions/67505252/plotly-box-p-value-significant-annotation
def add_p_value_annotation(fig, array_columns, subplot=None, _format=None):
    """
    Adds p-value annotations between box plots in a figure.

    Parameters:
    - fig: Plotly figure object containing box plots.
    - array_columns: Array specifying which columns to compare for p-values.
    - subplot: Specifies subplot to annotate (if applicable).
    - _format: Dictionary specifying formatting options for annotations.

    Returns:
    - fig: Modified Plotly figure object with p-value annotations.
    """
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01 + i * _format['interline'], 1.02 + i * _format['interline']]

    # Get values from figure
    fig_dict = fig.to_dict()

    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ''
        else:
            subplot_str = str(subplot)
        indices = []  # Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict['data']):
            # print(index, data['xaxis'], 'x' + subplot_str)
            if data['xaxis'] == 'x' + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ''

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair

        # Mare sure it is selecting the data and subplot you want
        #print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
        #print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

        # Get the p-value
        pvalue = stats.ttest_ind(
            fig_dict['data'][data_pair[0]]['y'],
            fig_dict['data'][data_pair[1]]['y'],
            equal_var=False,
        )[1]
        if pvalue >= 0.05:
            symbol = 'ns'
        elif pvalue >= 0.01:
            symbol = '*'
        elif pvalue >= 0.001:
            symbol = '**'
        else:
            symbol = '***'

        y_adjust = -0.1
        # Vertical line
        fig.add_shape(type="line",
                      xref="x" + subplot_str, yref="y" + subplot_str + " domain",
                      x0=column_pair[0], y0=y_range[index][0]+y_adjust,
                      x1=column_pair[0], y1=y_range[index][1]+y_adjust,
                      line=dict(color=_format['color'], width=2, )
                      )
        # Horizontal line
        fig.add_shape(type="line",
                      xref="x" + subplot_str, yref="y" + subplot_str + " domain",
                      x0=column_pair[0], y0=y_range[index][1]+y_adjust,
                      x1=column_pair[1], y1=y_range[index][1]+y_adjust,
                      line=dict(color=_format['color'], width=2, )
                      )
        # Vertical line
        fig.add_shape(type="line",
                      xref="x" + subplot_str, yref="y" + subplot_str + " domain",
                      x0=column_pair[1], y0=y_range[index][0]+y_adjust,
                      x1=column_pair[1], y1=y_range[index][1]+y_adjust,
                      line=dict(color=_format['color'], width=2, )
                      )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(dict(font=dict(color=_format['color'], size=14),
                                x=(column_pair[0] + column_pair[1]) / 2,
                                y=y_range[index][1] * _format['text_height']+y_adjust, #+ 0.075,
                                showarrow=False,
                                text=symbol,
                                textangle=0,
                                xref="x" + subplot_str,
                                yref="y" + subplot_str + " domain"
                                ))
    return fig

def create_bigwig_trace(filepath, plot_df):
    """
    Creates a trace from a BigWig file for plotting alongside genomic data.

    Parameters:
    - filepath: Path to the BigWig file.
    - plot_df: DataFrame containing genomic regions for plotting.

    Returns:
    - List of Plotly trace objects: Trace for the BigWig data and optional peak annotations.
    """
    if not os.path.isfile(filepath):
        print("Error: File not found")
        return None
    min_pos = plot_df['ref_position'].min()
    min_rel_pos= plot_df['rel_pos'].min()
    max_pos = plot_df['ref_position'].max()
    max_rel_pos = plot_df['rel_pos'].max()
    chromosome = plot_df['chrom'].iloc[0]  # Assuming all rows have the same chromosome
    # replace "CHROMOSOME" with "chr" in chromosome name
    chromosome = chromosome.replace("CHROMOSOME_", "chr")

    # Open the bigWig file
    bw = pyBigWig.open(filepath)

    # Extract the region of interest
    values = bw.values(chromosome, min_pos, max_pos)
    #print("values:",values)
    bw.close()

    # Find peaks
    peaks, _ = find_peaks(values)

    # Create the subplot for the bigWig data
    trace = go.Scatter(
        x=list(range(min_rel_pos, max_rel_pos)),
        y=values,
        mode='lines',
        name='BigWig Data',
        # set color to grey
        line=dict(color='grey', width=2)
    )

    traces = [trace]

    # Adding peak lines to the trace
    peak_traces = []
    for peak in peaks:
        peak_line = go.Scatter(
            x=[peak + min_rel_pos, peak + min_rel_pos],
            y=[min(values), max(values)],  # Assuming you want the line to span the full range of the y-axis
            mode='lines',
            name='Peak',
            showlegend=False,
            line=dict(color='grey', width=0.5, dash='dash')
        )
        traces.append(peak_line)

    return traces

def random_alpha_numeric(length=6):
    """
    Generates a random alphanumeric string.

    Parameters:
    - length: Desired length of the string.

    Returns:
    - str: Random alphanumeric string.
    """
    letters_and_digits = string.ascii_letters + string.digits
    return ''.join((random.choice(letters_and_digits) for i in range(length)))

### Plotting mod_base by chromosome
def prepare_chr_plotting_data(coverage_df):
    """
    Prepares data for plotting by grouping by condition and chromosome and computing m6A_frac.
    """
    # Group by condition and chromosome, and aggregate the data
    grouped_df = coverage_df.groupby(['condition', 'chromosome', 'start']).agg({
        'total_mod_base': 'sum',
        'total_canonical_base': 'sum'
    }).reset_index()

    # Compute the m6A/A ratio for the aggregated data
    grouped_df['m6A_frac'] = grouped_df['total_mod_base'] / grouped_df['total_canonical_base']

    # Create the "CHROMOSOME_X" and "Others" categories
    grouped_df['chromosome_group'] = grouped_df['chromosome'].apply(lambda x: "X" if x == 'CHROMOSOME_X' else 'Autosome')

    return grouped_df

def create_chr_type_box_plots_px(grouped_df, title_prefix):
    # Determine the number of unique conditions to set the number of columns for subplots
    unique_conditions = grouped_df['condition'].unique()
    num_conditions = len(unique_conditions)

    # Define color mapping for chromosome groups
    color_mapping = {'Autosome': 'blue', 'X': 'red'}

    # Create subplots with one column for each condition
    fig = make_subplots(rows=1, cols=num_conditions, subplot_titles=unique_conditions)

    # Iterate over each condition to create a separate plot
    for i, condition in enumerate(unique_conditions, start=1):
        # Filter the DataFrame for the current condition
        condition_data = grouped_df[grouped_df['condition'] == condition]

        # Create a box plot using plotly express
        box_fig = px.box(
            condition_data,
            x='chromosome_group',
            y='m6A_frac',
            color='chromosome_group',
            color_discrete_map=color_mapping,
            title=f"{title_prefix} - {condition}",
            category_orders={"chromosome_group": ["Autosome", "X"]}, # Ensure consistent order
            template='plotly_white'
        )

        # Extract the traces from the box plot and add them to the corresponding subplot
        for trace in box_fig['data']:
            fig.add_trace(trace, row=1, col=i)

    # Update layout for the entire figure
    fig.update_layout(
        title=f"{title_prefix} m6A/A Ratio by Chromosome Group for Each Condition",
        yaxis_title="m6A/A Ratio",
        template="plotly_white",
        height=600,
        width=400 * num_conditions,  # Adjust width based on the number of conditions
        showlegend=False
    )

    # Update y-axis format for all plots
    fig.update_yaxes(tickformat='.1%')  # Ensure y-axis formatting is consistent across all plots
    fig.show()

def create_condition_wise_chromosome_box_plots(grouped_df, title_prefix):
    # Determine the number of unique conditions and chromosomes
    unique_conditions = grouped_df['condition'].unique()
    unique_chromosomes = grouped_df['chromosome'].unique()

    # Create a subplot figure with 1 row and a number of columns equal to the number of conditions
    fig = make_subplots(rows=1, cols=len(unique_conditions),
                        subplot_titles=[f"{title_prefix} - {condition}" for condition in unique_conditions],
                        horizontal_spacing=0.05)  # Adjust spacing as needed

    # Generate a color sequence for the chromosomes using Plotly Express
    color_sequence = px.colors.qualitative.Alphabet

    # Iterate over each condition to create a box plot for all chromosomes within that condition
    for i, condition in enumerate(unique_conditions, start=1):
        # Filter the DataFrame for the current condition
        condition_data = grouped_df[grouped_df['condition'] == condition]

        # For each condition, add a box plot for each chromosome
        for j, chrom in enumerate(unique_chromosomes, start=0):
            # Filter the DataFrame for the current chromosome within the condition
            chrom_data = condition_data[condition_data['chromosome'] == chrom]

            fig.add_trace(
                go.Box(
                    y=chrom_data['m6A_frac'],
                    name=chrom,
                    boxpoints='all',  # Show all points
                    jitter=0.4,  # Slight jitter so points don't overlap
                    notched=True,  # Notch to show median confidence interval
                    marker=dict(color=color_sequence[j % len(color_sequence)]),
                ),
                row=1, col=i
            )

    # Update layout for the entire figure
    fig.update_layout(
        title=f"{title_prefix} m6A/A Ratio by Chromosome for Each Condition",
        yaxis_title="m6A/A Ratio",
        template="plotly_white",
        height=600,  # Height of the entire figure
        width=400 * len(unique_conditions),  # Width of the entire figure
        showlegend=True
    )

    # Update y-axis format for all subplots
    fig.update_yaxes(tickformat='.1%', row=1)

    # Display the figure
    fig.show()

def process_region_for_summary(args):
    chromosome, start, end, _, _, strand, each_bam, each_condition, each_thresh, each_exp_id, sampling_frac = args
    #print("starting on chromosome:", chromosome," | start:", start,"| end:",end, "with thresh:", each_thresh, "and condition:", each_condition, "and exp_id:", each_exp_id)
    temp_df = get_summary_from_bam(sampling_frac, each_thresh, "/Data1/software/modkit/modkit", each_bam, each_condition, each_exp_id, 12, chromosome, start, end)
    temp_df["chromosome"] = chromosome
    temp_df["start"] = start
    temp_df["end"] = end
    return temp_df


def process_and_export_summary_by_region(sampling_frac, type_selected,new_bam_files, conditions, thresh_list, exp_ids, modkit_bed_df, force_replace=False, output_directory="temp_files"):
    """
    Processes BAM files by region, summarizes them, and exports the summary to a CSV file.

    Parameters:
    - new_bam_files: List of paths to BAM files.
    - conditions: List of conditions corresponding to each BAM file.
    - thresh_list: List of threshold values for processing.
    - exp_ids: List of experiment IDs.
    - modkit_bed_df: DataFrame containing regions to process (chromosome, start, end, etc.).
    - force_replace: Whether to force regeneration of the summary file even if it exists.
    - output_directory: Directory where the output CSV file will be saved.
    """

    # Define filename for summary table based on selected conditions
    summary_table_name = os.path.join(output_directory, "_".join([conditions[0], conditions[-1], str(sampling_frac), "thresh"+str(thresh_list[0]), type_selected[0]]) + "_summary_table.csv")

    if not force_replace and os.path.exists(summary_table_name):
        print("Summary table exists, importing...")
        summary_bam_df = pd.read_csv(summary_table_name, sep="\t", header=0)
    else:
        # Prepare the arguments for multiprocessing
        args_list = [(row[0], row[1], row[2], row[3], row[4], row[5], each_bam, each_condition, each_thresh, each_exp_id, sampling_frac)
                     for _, row in modkit_bed_df.iterrows()
                     for each_bam, each_condition, each_thresh, each_exp_id in zip(new_bam_files, conditions, thresh_list, exp_ids)]

        # Use multiprocessing to process regions in parallel
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(process_region_for_summary, tqdm(args_list, total=len(args_list)))

        # Concatenate the results
        summary_bam_df = pd.concat(results, ignore_index=True)

        # Export the concatenated DataFrame to CSV
        summary_bam_df.to_csv(summary_table_name, sep="\t", header=True, index=False)
        print(f"Exported summary to {summary_table_name}")

    return summary_bam_df

def create_coverage_df(code, summary_bam_df, coverage_df_name):
    """
    Creates a coverage dataframe for a specified code ('a' or 'm') and saves it.
    """
    # Filter summary_bam_df based on code and canonical base (A for 'a', C for 'm')
    base_filter = "A" if code == 'a' else "C"
    filtered_df_canonical = summary_bam_df[summary_bam_df['base'] == base_filter]
    # filter based on code
    filtered_df_mod = summary_bam_df[summary_bam_df['code'] == code]

    # Calculate total_canonical_base
    total_canonical_base = filtered_df_canonical.groupby(['condition', 'chromosome', 'start'])['pass_count'].sum().reset_index()
    total_canonical_base.rename(columns={'pass_count': 'total_canonical_base'}, inplace=True)

    # Calculate total_mod_base
    total_mod_base = filtered_df_mod.groupby(['condition', 'chromosome', 'start'])['pass_count'].sum().reset_index()
    total_mod_base.rename(columns={'pass_count': 'total_mod_base'}, inplace=True)


    # Merge total_mod_base and total_canonical_base DataFrames
    coverage_df = pd.merge(total_mod_base, total_canonical_base, on=['condition', 'chromosome', 'start'], how='outer').fillna(0)

    # Calculate m6A_frac
    coverage_df['mod_frac'] = coverage_df['total_mod_base'] / (coverage_df['total_mod_base'] + coverage_df['total_canonical_base'])

    # Save coverage df
    coverage_df.to_csv(coverage_df_name, sep="\t", header=True, index=False)
    return coverage_df

