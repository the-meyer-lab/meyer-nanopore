### NANOTOOLS
# File description: Collection of python functions for analysis of nanopore datasets
#
# Created: Jul 28th 2023
# By: Yuri Malina
# ymalina@berkeley.edu

# --------------------------------------------------------------------
### IMPORTS
import pandas as pd
import numpy as np
import subprocess
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from io import StringIO
import inspect
import random
import pysam
from scipy import stats
import multiprocessing # used for parallel processing
from multiprocessing import Pool # used for parallel processing
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
from typing import Optional
import pyBigWig
from scipy.signal import find_peaks
import string


def filter_bed_file(bed_file,sample_source,selection,chromosome_selected,chr_type_selected,type_selected,strand_selected,max_regions,bed_window):
    """ Function takes in a bed file and filters it based on the following criteria. Saves 1 bed file per sample source type.
    to the path of the original bedfile with the folllowing name: "/Data1/reference/temp_do_not_use_"+each_type+".bed.gz"
    # Required fields:
    * bed_file = path to complete bedfile, e.g. "/Data1/reference/tss_tes_rex_combined_v2.bed"
    ** Must have the following columns: "chromosome","start","end","strand","type","chr-type"
    * sample_source = "type" "chr_type" or "chromosome" (which column based on which to filter bed file.)
    * selection = list of types to select, e.g. ["tss","tes"]
    * chromosome_selected = list of chromosomes to select, e.g. ["chrI","chrII"]
    * chr_type_selected = list of chr-types to select, e.g. ["chrI_telo","chrII_telo"]
    * type_selected = list of types to select, e.g. ["tss","tes"]
    * strand_selected = list of strands to select, e.g. ["+","-"]
    * max_regions = maximum number of regions to select per type, e.g. 1000
    * bed_window = number of basepairs to add to each side of the bed file, e.g. 1000
    """
    ### Select bed file
    full_bed = pd.read_csv(bed_file,sep='\t')
    bed=[]

    # Select rows where type == "whole_chr" and store end value to chromosome_ends
    chromosome_ends = full_bed[full_bed["type"] == "whole_chr"]

    for each_type in selection:
    # REGION CONFIGURATION
        if sample_source == "type":
            temp_bed = full_bed[full_bed["chromosome"].isin(chromosome_selected) &
                                full_bed["chr-type"].isin(chr_type_selected) &
                                full_bed["type"].str.contains(each_type) &
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

        # Adjust by bed_window, generally used for TSS or other elements with 0 width.
        temp_bed["start"]=temp_bed["start"] - bed_window
        temp_bed["end"]=temp_bed["end"] + bed_window
        temp_bed.loc[temp_bed["start"]<0,"start"]=0

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

        display_sample_rows(temp_bed, 5)
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
    """ Function takes in a .bam file and a .bed file and extracts only the reads that overlap the regions in the .bed file.
    # Required fields:
    * bam_file = path to .bam file
    * condition = condition name
    * bam_frac = fraction of reads to keep
    * selection = list of types to select
    * m6A_thresh = m6A threshold
    * output_stem = path to output .bam file"""
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


def extract_m6A_per_region(bam_file, bed_file, threshold, condition):
    """ Function takes in a .bam file, a .bed file, and a threshold and extracts the number of m6A sites per region.
    Required fields:
    * bam_file = path to .bam file
    * bed_file = path to .bed file
    * threshold = m6A threshold
    * condition = condition name"""
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

def extract_m6A_per_region_parellized(bam_file, condition,bam_frac,file_prefix, selection, m6A_thresh, output_stem,new_bed_files):
    # for each_bedfile and each_selection in zip(new_bed_files,selection):
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

# Define function that takes in a danpos nucleosome position data frame and returns
# Frequency table of nucleosome repeat length.
# Def frequency_table function with input df and optional n_dist parameter:
def frequency_table(df, n_dist=None):
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

# Define function that takes in a nucleosome RPL frequency tables, plots the distribution
# and returns the peak values as a dataframe.
def plot_NRL_dist(data, chromosome_name, color, output_prefix,smoothing_val = None):
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

def plot_NRL_dist_compare(data, chromosome_name, output_prefix,smoothing_val = None,norm_bool=True,hue='condition',window=1000):
    # Clear .plt
    #plt.clf()
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

def plot_NRL_regression(dataframe, title, output_prefix,equation_table):
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
def get_summary_from_bam(sampling_frac: float, a_threshold: float, modkit_path: str,bam_path: str,each_condition: str,each_exp_id: str, chromosome: Optional[str] = None, start: Optional[int] = None, end: Optional[int] = None) -> pd.DataFrame:
    """
    Fetches data from the given bam file using modkit and returns it as a pandas DataFrame.

    Parameters:
    - m_threshold (float): Threshold for 'm'.
    - a_threshold (float): Threshold for 'a'.
    - bam_path (str): Path to the bam file.

    Returns:
    - DataFrame: Pandas DataFrame containing the fetched data.

    "--mod-thresholds",
        f"m:{m_threshold}",
        "--mod-thresholds",
        f"a:{a_threshold}",

    """

    command = [
        modkit_path,
        "summary",
        "--ignore",
        "m",
        "--threads",
        "6",
        "--sampling-frac",
        f"{sampling_frac}",
        "--seed",
        "10",
        #"--mod-thresholds",
        #f"a:{a_threshold}",
        #"--no-filtering",
        "--filter-threshold",
        f"A:{1-a_threshold}",
        "--mod-thresholds",
        f"a:{a_threshold}",
        "--log-filepath",
        f"temp_files/modkit_summary_{each_condition}_{each_exp_id}.log"
    ]

    # if chromosome, start and end are not None, add the region argument
    if chromosome is not None and start is not None and end is not None:
        region_argument = f"{chromosome}:{start}-{end}"
        command.extend(["--region", region_argument])

    command.extend([bam_path])

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

#take in a list of bed file paths and return a single bed file with all regions formated for modkit
# also take a boolean on whether to down sample autosomes, and whether to select opp strands for regions
# saves output to temp file
def generate_modkit_bed(new_bed_files, down_sample_autosome, select_opp_strand, output_name):
    # Initialize an empty DataFrame to store the combined data from all bed files
    combined_bed_df = pd.DataFrame()

    # Read each bed file and append it to the combined DataFrame
    for each_bed in new_bed_files:
        bed_path = each_bed[:-3]
        temp_df = pd.read_csv(bed_path, sep="\t", header=None)
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

### Functions to return sample of read lengths to a dataframe, used for calculating N50 statistics
def extract_n50_from_bam(bam_file, fraction_to_sample=0.01):
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
    bam_file, condition, exp_id = args
    print(f"Starting on: {bam_file}, {condition}")
    temp_df = extract_n50_from_bam(bam_file, 0.1)
    temp_df['condition'] = condition
    temp_df['exp_id'] = exp_id
    return temp_df

# Function to calculate N50
def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    half_total_length = sum(sorted_lengths) / 2.0
    current_length = 0
    for length in sorted_lengths:
        current_length += length
        if current_length >= half_total_length:
            return length

def calculate_and_plot_n50(new_bam_files, conditions, exp_ids):
    # Your original process_bam_file function

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
def add_p_value_annotation(fig, array_columns, subplot=None,
                           _format=dict(interline=0.07, text_height=1.07, color='black')):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)

    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
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
    # check that file_path exists and is a bigWig file otherwise return error
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

# generate random alphanumeric code with default 6 characters
def random_alpha_numeric(length=6):
    letters_and_digits = string.ascii_letters + string.digits
    return ''.join((random.choice(letters_and_digits) for i in range(length)))