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
#import pyranges as pr

# Set seaborn style for all plots
sns.set(style="whitegrid")

def _autocorr_vec(x: np.ndarray, max_lag: int) -> np.ndarray:
    """
    True Pearson-style autocorrelation up to max_lag:
     - pairwise-complete (ignores NaNs)
     - per-lag separate means & variances
    """
    n = len(x)
    ac = np.full(max_lag+1, np.nan)
    if n < 2:
        return ac

    for k in range(0, min(max_lag, n-1) + 1):
        x1 = x[:n-k]
        x2 = x[k:]
        valid = ~np.isnan(x1) & ~np.isnan(x2)
        if valid.sum() < 2:
            continue

        x1v = x1[valid]
        x2v = x2[valid]
        m1 = x1v.mean()
        m2 = x2v.mean()

        num = np.sum((x1v - m1) * (x2v - m2))
        v1  = np.sum((x1v - m1)**2)
        v2  = np.sum((x2v - m2)**2)
        denom = np.sqrt(v1 * v2)

        ac[k] = num/denom if denom > 0 else np.nan

    return ac

def get_color(key):
    """
    Return a hex colour for *key*.

    Priority:
      1) Condition keywords (dpy27/sdc2/n2*/sdc3/dpy21)  → fixed colors
      2) Types: MEX_motif / MEXII_motif / univ_nuc       → fixed colors
      3) Degron keys and other fallbacks                 → existing logic
    """
    from plotly.colors import qualitative, sequential
    import re

    # ── helpers ───────────────────────────────────────────────────────────────
    def _rgb_to_hex(s):
        m = re.fullmatch(r"rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", s)
        if m:
            r, g, b = map(int, m.groups())
            return f"#{r:02X}{g:02X}{b:02X}"
        return s

    def lighten_color(hexstr, amount=0.7):
        hexstr = hexstr.lstrip('#')
        r, g, b = [int(hexstr[i:i+2], 16) for i in (0, 2, 4)]
        r = int(r + (255 - r) * amount)
        g = int(g + (255 - g) * amount)
        b = int(b + (255 - b) * amount)
        return f"#{r:02X}{g:02X}{b:02X}"

    PALETTE = qualitative.Plotly
    k_lower = str(key).lower()

    # ── 1) condition overrides (unchanged) ────────────────────────────────────
    if "dpy27" in k_lower:   return "#47B562"
    if "sdc2"  in k_lower:   return "#b12537"
    if "n2_old" in k_lower:  return "#4974a5"
    if "n2_mid" in k_lower:  return "#61baad"
    if "n2_young" in k_lower:return "#dea01e"
    if "sdc3"  in k_lower:   return "#9d55ac"
    if "dpy21" in k_lower:   return "#808080"

    # ── 2) type colours ───────────────────────────────────────────────────────
    TYPE_COLOR_MAP = {
        "mex_motif":   "#1f77b4",  # blue
        "mexii_motif": "#ff7f0e",  # orange
        "univ_nuc":    "#2ca02c",  # green
    }
    for token, hexcol in TYPE_COLOR_MAP.items():
        if token in k_lower:
            return hexcol

    # ── 3) degron time-course palette (as before) ─────────────────────────────
    degron_keys = [
        "T0_rep2", "T0_rep3", "T0p5_rep3", "T1_rep3",
        "T1p5_rep3", "T2_rep2", "T3_rep2", "T4_rep3",
    ]
    ylorrd_hex = [_rgb_to_hex(c) for c in sequential.YlOrRd[::-1][:8]]
    DEGRON_MAP = dict(zip(degron_keys, ylorrd_hex))

    base = None
    if k_lower.startswith("sdc2_degron_mid_t"):
        parts = str(key).split("_")
        if len(parts) >= 5:
            short_key = f"{parts[3]}_{parts[4]}"
            base = DEGRON_MAP.get(short_key)

    # ── keyword-based palette for all others (unchanged) ──────────────────────
    if base is None:
        if "mexiiscramble" in k_lower:
            base = PALETTE[8]       # pink
        elif "4thctog" in k_lower:
            base = "#7F7F7F"        # grey
        else:
            base = PALETTE[6]       # default fallback

    if "intergenic" in k_lower:
        base = lighten_color(_rgb_to_hex(base), amount=0.25)

    return base


def get_colors(keys):
    """
    Accept a single key or an iterable of keys and return
    a single color or list of colors.
    """
    if hasattr(keys, "__iter__") and not isinstance(keys, str):
        return [get_color(k) for k in keys]
    else:
        return get_color(keys)




# Function to impute and smooth scores for a bedgraph dataframe:
def process_chromosome(chr_data):
    chrom, scores, weights, impute_window, smooth_window, fill_value = chr_data
    #print("Imputing chromosome:", chrom)

    # Imputation
    if impute_window > 0:
        # Calculate weighted means using convolution
        weights_conv = np.convolve(weights.values, np.ones(impute_window), mode='same')
        weighted_mean = np.convolve(scores.values * weights.values, np.ones(impute_window), mode='same') / weights_conv
        imputed_scores = scores.combine_first(pd.Series(weighted_mean, index=scores.index)).fillna(fill_value)
        imputed_weights = weights.combine_first(
            weights.rolling(window=impute_window, center=True, min_periods=1).mean()).fillna(fill_value)
    else:
        imputed_scores = scores.fillna(fill_value)
        imputed_weights = weights.fillna(fill_value)

    #print("Smoothing chromosome:", chrom)
    # Smoothing
    if smooth_window > 0:
        # Calculate weighted sums and counts using convolution
        weights_conv_smooth = np.convolve(imputed_weights.values, np.ones(smooth_window), mode='same')
        weighted_sum = np.convolve(imputed_scores.values * imputed_weights.values, np.ones(smooth_window), mode='same')
        weighted_count = np.convolve(imputed_weights.values, np.ones(smooth_window), mode='same')

        # Avoid division by zero
        weighted_count[weighted_count == 0] = np.nan

        # Compute smoothed scores
        smoothed_scores = pd.Series(weighted_sum / weighted_count, index=scores.index)
        smoothed_scores = smoothed_scores.combine_first(
            imputed_scores.rolling(window=smooth_window, center=True, min_periods=1).mean()).fillna(fill_value)

        # Smoothed weights (just mean of weights over window)
        smoothed_weights = pd.Series(weights_conv_smooth / smooth_window, index=weights.index)
    else:
        smoothed_scores = imputed_scores
        smoothed_weights = imputed_weights

    return chrom, imputed_scores, imputed_weights, smoothed_scores, smoothed_weights

def parallel_impute_and_smooth(data, impute_window, smooth_window, fill_value=0):
    #print("DataFrame columns:", data.columns)
    
    # Prepare data for parallel processing
    chr_data = [
        (chrom, group['score'], group['coverage'], impute_window, smooth_window, fill_value)
        for chrom, group in data.groupby('chromosome')
    ]
    print("Starting parallel processing of chromosomes...")
    # Use multiprocessing.Pool for parallel processing
    results = [process_chromosome(ch_data) for ch_data in chr_data]

    # Combine results
    imputed_scores = pd.concat([result[1] for result in results])
    imputed_weights = pd.concat([result[2] for result in results])
    smoothed_scores = pd.concat([result[3] for result in results])
    smoothed_weights = pd.concat([result[4] for result in results])
    # print head of imputed_scores and smoothed_scores
    #print("Smoothed scores head:", smoothed_scores.head())
    # print head of imputed_weights and smoothed_weights
    #print("Smoothed weights head:", smoothed_weights.head())
    return imputed_scores, imputed_weights, smoothed_scores, smoothed_weights

def get_chromosome_sizes(genome_fasta_path):
    # Ensure the .fai index file exists
    if not os.path.exists(genome_fasta_path + '.fai'):
        subprocess.run(['samtools', 'faidx', genome_fasta_path], check=True)

    # Read the .fai file to get chromosome sizes
    chromosome_sizes = {}
    with open(genome_fasta_path + '.fai', 'r') as fai_file:
        for line in fai_file:
            parts = line.strip().split('\t')
            chromosome_sizes[parts[0]] = int(parts[1])

    return chromosome_sizes


# Assuming get_chromosome_sizes is defined elsewhere
# from your_module import get_chromosome_sizes

# Function to load and process bedgraph data into a dataframe
def load_bedgraph_file(file_path, chr=None, start=None, end=None, nan_fill=True):
    # Load data
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chromosome', 'start', 'end', 'score', 'coverage'])

    # fill "coverage" column with 1s if it is empty or not present or filled with nan
    if 'coverage' not in df.columns or df['coverage'].isnull().all():
        df['coverage'] = 1

    # Get chromosome sizes using samtools (you may need to adjust this part)
    chromosome_sizes = get_chromosome_sizes("/Data1/reference/c_elegans.WS235.genomic.fa")

    if chr is not None and start is not None and end is not None:
        # Filter rows based on chr, start, end
        df = df[(df['chromosome'] == chr) & (df['start'] > start) & (df['start'] < end)]

    if nan_fill:
        # Create a DataFrame with all possible positions for the entire genome
        all_positions = pd.DataFrame({
            'chromosome': np.concatenate([[chr] * length for chr, length in chromosome_sizes.items()]),
            'start': np.concatenate([np.arange(0, length ) for length in chromosome_sizes.values()]),
            'end': np.concatenate([np.arange(1, length + 1) for length in chromosome_sizes.values()])
        })

        # Merge with original DataFrame
        df = pd.merge(all_positions, df, on=['chromosome', 'start', 'end'], how='left')

        # Fill NA values
        df['score'].fillna(np.nan, inplace=True)
        df['coverage'].fillna(np.nan, inplace=True)

    return df

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
    Converts standard bed files into a BED formatted for downstream
    look‑ups.

    Parameters
    ----------
    new_bed_files : list[str]
        Paths to <*.bed.gz> files produced by `filter_bed_file`.

    Returns
    -------
    pandas.DataFrame
        Combined BED (cols: chrom, bed_start, bed_end, bed_strand,
        type, chr_type).
    """
    import os
    import pandas as pd

    combined_bed_df = pd.DataFrame()

    # ==================================================================
    #  Read each BED (un‑gzipped), skipping the deliberately empty ones
    # ==================================================================
    for each_bed in new_bed_files:
        bed_path = each_bed[:-3]       # strip ".gz"

        # ─── PATCH START – skip empty BEDs ───────────────────────────
        try:
            if os.path.getsize(bed_path) == 0:
                print(f"skipping empty file … {bed_path}")
                continue
            bed_df = pd.read_csv(bed_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            print(f"skipping empty file … {bed_path}")
            continue
        # ─── PATCH END ───────────────────────────────────────────────

        combined_bed_df = combined_bed_df.append(bed_df)

    # Ensure we actually read something
    if combined_bed_df.empty:
        raise RuntimeError("No non‑empty BED files found – "
                           "combined lookup BED would be empty.")

    combined_bed_df.columns = [
        "chrom", "bed_start", "bed_end",
        "bed_strand", "type", "chr_type"
    ]
    combined_bed_df.sort_values(
        by=["chrom", "bed_start"],
        ascending=True,
        inplace=True
    )
    combined_bed_df.dropna(inplace=True)
    combined_bed_df.drop_duplicates(inplace=True)
    combined_bed_df.reset_index(drop=True, inplace=True)

    return combined_bed_df



def filter_bed_file(bed_file, sample_source, selection,
                    chromosome_selected, chr_type_selected, type_selected,
                    strand_selected, max_regions, bed_window,
                    intergenic_window):
    """
    Filters a bed file based on various criteria and saves 1 bed file per
    sample source type.  Applies max_regions limit per type within each
    selection.

    Parameters
    ----------
    bed_file : str
        Path to the complete bed file.
    sample_source : {"type", "chr_type", "chromosome"}
        Column to iterate over.
    selection : list[str]
        Values drawn from *sample_source* to keep (processed one‑by‑one).
    chromosome_selected, chr_type_selected, type_selected : list[str]
        Additional filtering criteria.
    strand_selected : list[str]
        Strands to retain.
    max_regions : int
        Maximum regions per *type* (0 ⇒ no limit).
    bed_window, intergenic_window : int
        Windows added to bed coordinates.

    Returns
    -------
    list[str]
        Paths to the bgzipped bed files generated.
    """
    print("Filtering bed file...")
    print("Configs:", sample_source, selection, chromosome_selected,
          chr_type_selected, type_selected, strand_selected,
          max_regions, bed_window, intergenic_window)

    # ------------------------------------------------------------------
    #  Load the full BED; extract chromosome endpoints
    # ------------------------------------------------------------------
    full_bed = pd.read_csv(bed_file, sep="\t")
    chromosome_ends = full_bed[full_bed["type"] == "whole_chr"]
    print("Chromosome ends:", chromosome_ends)

    bed = []          # collect <*.bed.gz> paths here

    # ==================================================================
    #  Iterate through “X”, “Autosome”, etc.
    # ==================================================================
    for each_type in selection:
        # ───── REGION CONFIGURATION ──────────────────────────────────
        if sample_source == "type":
            temp_bed = full_bed[
                full_bed["chromosome"].isin(chromosome_selected) &
                full_bed["chr-type"].isin(chr_type_selected) &
                full_bed["type"].eq(each_type) &
                full_bed["strand"].isin(strand_selected)
            ]
        elif sample_source == "chr_type":
            temp_bed = full_bed[
                full_bed["chromosome"].isin(chromosome_selected) &
                full_bed["chr-type"].str.contains(each_type) &
                full_bed["type"].isin(type_selected) &
                full_bed["strand"].isin(strand_selected)
            ]
        elif sample_source == "chromosome":
            temp_bed = full_bed[
                full_bed["chromosome"].eq(each_type) &
                full_bed["chr-type"].isin(chr_type_selected) &
                full_bed["type"].isin(type_selected) &
                full_bed["strand"].isin(strand_selected)
            ]

        temp_bed = temp_bed.copy()   # avoid SettingWithCopyWarning

        # ─── PATCH 1 START – handle *totally empty* temp_bed ─────────
        if temp_bed.empty:
            print(f"No regions found for “{each_type}”.  Creating empty BED.")
            bed_file_path   = os.path.dirname(bed_file)
            temp_bedfile    = f"{bed_file_path}/{each_type}.bed"
            temp_bedfile_gz = f"{bed_file_path}/{each_type}.bed.gz"
            temp_bedfile_tbi = f"{bed_file_path}/{each_type}.bed.gz.tbi"

            pd.DataFrame(columns=full_bed.columns).to_csv(
                temp_bedfile, sep="\t", header=False, index=False
            )

            with open(temp_bedfile_gz, "wb") as out_fh:
                res = subprocess.run(
                    ["bgzip", "-c", temp_bedfile],
                    stdout=out_fh, stderr=subprocess.PIPE
                )
            if res.returncode != 0:
                print("bgzip error:", res.stderr.decode())
            else:
                print(f"{temp_bedfile} → {temp_bedfile_gz}")

            with open(temp_bedfile_tbi, "wb") as out_fh:
                res2 = subprocess.run(
                    ["tabix", "-f", "-p", "bed", temp_bedfile_gz],
                    stdout=out_fh, stderr=subprocess.PIPE
                )
            if res2.returncode != 0:
                print("tabix warning:", res2.stderr.decode())
            else:
                print(f"Index created: {temp_bedfile_tbi}")

            bed.append(temp_bedfile_gz)
            continue          # move to next each_type
        # ─── PATCH 1 END ─────────────────────────────────────────────

        # ───────────────────── apply max_regions ─────────────────────
        if max_regions > 0:
            type_groups = []
            for type_name in temp_bed["type"].unique():
                subset = temp_bed[temp_bed["type"] == type_name]
                excess = len(subset) - max_regions
                if excess > 0:
                    keep_idx = np.random.choice(
                        subset.index, size=max_regions, replace=False
                    )
                    subset = subset.loc[keep_idx]
                type_groups.append(subset)

            # ─── PATCH 2 START – guard against empty *type_groups* ───
            if not type_groups:          # nothing survived filtering
                print(f"All regions dropped for “{each_type}” "
                      f"after applying max_regions={max_regions}. "
                      "Creating empty BED.")
                bed_file_path   = os.path.dirname(bed_file)
                temp_bedfile    = f"{bed_file_path}/{each_type}.bed"
                temp_bedfile_gz = f"{bed_file_path}/{each_type}.bed.gz"
                temp_bedfile_tbi = f"{bed_file_path}/{each_type}.bed.gz.tbi"

                pd.DataFrame(columns=full_bed.columns).to_csv(
                    temp_bedfile, sep="\t", header=False, index=False
                )

                with open(temp_bedfile_gz, "wb") as out_fh:
                    subprocess.run(
                        ["bgzip", "-c", temp_bedfile],
                        stdout=out_fh, stderr=subprocess.PIPE, check=False
                    )
                with open(temp_bedfile_tbi, "wb") as out_fh:
                    subprocess.run(
                        ["tabix", "-f", "-p", "bed", temp_bedfile_gz],
                        stdout=out_fh, stderr=subprocess.PIPE, check=False
                    )
                bed.append(temp_bedfile_gz)
                continue
            # ─── PATCH 2 END ─────────────────────────────────────────

            temp_bed = pd.concat(type_groups, ignore_index=True)

        # ─────────────────── adjust windows, clip ends ───────────────
        temp_bed.sort_values(by=["chromosome", "start"],
                             inplace=True, ascending=True)
        temp_bed.reset_index(drop=True, inplace=True)

        if "intergenic_control" in temp_bed["type"].unique():
            temp_bed["start"] = np.where(
                temp_bed["type"] != "intergenic_control",
                temp_bed["start"] - bed_window,
                temp_bed["start"] - intergenic_window
            )
            temp_bed["end"] = np.where(
                temp_bed["type"] != "intergenic_control",
                temp_bed["end"] + bed_window,
                temp_bed["end"] + intergenic_window
            )
        else:
            temp_bed["start"] -= bed_window
            temp_bed["end"]   += bed_window

        temp_bed["start"] = temp_bed["start"].clip(lower=0)

        for chrom in chromosome_ends["chromosome"].unique():
            chrom_end = chromosome_ends.loc[
                chromosome_ends["chromosome"].eq(chrom), "end"
            ].values[0]
            mask = temp_bed["chromosome"].eq(chrom)
            temp_bed.loc[mask, "end"] = temp_bed.loc[mask, "end"].clip(
                upper=chrom_end
            )

        # ───────────────────── write / bgzip / tabix ─────────────────
        bed_file_path   = os.path.dirname(bed_file)
        temp_bedfile    = f"{bed_file_path}/{each_type}.bed"
        temp_bedfile_gz = f"{bed_file_path}/{each_type}.bed.gz"
        temp_bedfile_tbi = f"{bed_file_path}/{each_type}.bed.gz.tbi"

        temp_bed.to_csv(temp_bedfile, sep="\t", header=False, index=False)

        with open(temp_bedfile_gz, "wb") as out_fh:
            subprocess.run(
                ["bgzip", "-c", temp_bedfile],
                stdout=out_fh, stderr=subprocess.PIPE, check=False
            )
        with open(temp_bedfile_tbi, "wb") as out_fh:
            subprocess.run(
                ["tabix", "-f", "-p", "bed", temp_bedfile_gz],
                stdout=out_fh, stderr=subprocess.PIPE, check=False
            )

        bed.append(temp_bedfile_gz)

    print("Saved the following bedfiles:", bed)
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
def get_summary_from_bam(sampling_frac, a_threshold, modkit_path, bam_path, each_condition, each_exp_id,c_threshold = None, thread_ct=4, chromosome=None, start=None, end=None, bed_file=None):
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

    # if c_thresh is not provided, set it to a_threshold
    if c_threshold is None:
        c_threshold = a_threshold

    command = [
        modkit_path,
        "summary",
        #"--ignore",
        #"m",
        "--tsv",
        "--threads",
        f"{thread_ct}",
        # "--sampling-frac",
        # f"{sampling_frac}",
        "--seed",
        "10",
        #"--no-filtering",
        "--filter-threshold",
        f"A:{0.7}",
        "--filter-threshold",
        f"C:{0.7}",
        "--mod-thresholds",
        f"a:{a_threshold}",
        "--mod-thresholds",
        f"m:{c_threshold}",
        "--log-filepath",
        f"temp_files/modkit_summary_{each_condition}_{each_exp_id}.log"
    ]

    # sampling
    if sampling_frac is None or sampling_frac >= 1:
        command.append("--no-sampling")
    else:
        command += ["--sampling-frac", str(sampling_frac)]  # <-- correct flag

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
    print("Running command:", " ".join(command))
    result = subprocess.run(command,  text=True, capture_output=True)
    if chromosome is None:
        print('stdout:')
        print(result.stdout)

    # Extract only the table data starting from the column names
    # Split the stdout using the substring
    # ------------- replace everything from "split_data = ..." down with:
    from io import StringIO
    df = pd.read_csv(StringIO(result.stdout), sep="\t")

    if df.empty:
        print("Warning: modkit returned no rows; stdout was:")
        print(result.stdout[:300])  # first 300 chars for debug
        return pd.DataFrame(columns=['condition', 'exp_id'])

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

def generate_modkit_bed(new_bed_files,
                        down_sample_autosome: bool,
                        select_opp_strand: bool,
                        output_name: str):
    """
    Build the “include‑bed” for `modkit` by concatenating the filtered
    BEDs produced by *filter_bed_file*.

    Parameters
    ----------
    new_bed_files : list[str]
        Paths to *.bed.gz* files returned by `filter_bed_file`.
    down_sample_autosome : bool
        If True, randomly down‑sample autosome regions to match the
        number of X‑chromosome regions.
    select_opp_strand : bool
        If True, create a second entry for each region on the opposite
        strand.
    output_name : str
        Name of the final BED written to disk.

    Returns
    -------
    pandas.DataFrame
        The concatenated BED as a DataFrame.
    """
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    dfs = []

    # ==================================================================
    #  Iterate through each gzipped BED emitted by *filter_bed_file*
    # ==================================================================
    for each_bed in new_bed_files:
        # Strip “.gz” to get the plain‑text BED path
        bed_path = each_bed[:-3]

        # ──────────────────────────────────────────────────────────────
        #                PATCH  – Gracefully skip empty files
        # ──────────────────────────────────────────────────────────────
        try:
            if os.path.getsize(bed_path) == 0:
                print(f"skipping empty file … {bed_path}")
                continue
            temp_df = pd.read_csv(bed_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            print(f"skipping empty file … {bed_path}")
            continue
        # ─────────────── PATCH END ────────────────────────────────────

        # Ensure we have exactly 6 BED columns
        if temp_df.shape[1] < 6:
            raise ValueError(f"{bed_path} has <6 columns; "
                             "unexpected BED format")

        # Optionally down‑sample autosome regions so that the number of
        # autosome rows equals the number of X rows per chromosome type.
        if down_sample_autosome:
            is_x = temp_df[0].str.contains("X")
            num_x = is_x.sum()
            num_auto = len(temp_df) - num_x
            if num_auto > num_x:
                auto_idx = temp_df.loc[~is_x].index
                drop_n = num_auto - num_x
                drop_idx = np.random.choice(auto_idx, drop_n, replace=False)
                temp_df = temp_df.drop(drop_idx).reset_index(drop=True)

        # Optionally duplicate entries with the strand flipped
        if select_opp_strand:
            # Flip “+” ↔ “−” in column 5
            flipped = temp_df.copy()
            flipped[5] = flipped[5].replace({"+": "-", "-": "+"})
            temp_df = pd.concat([temp_df, flipped], ignore_index=True)

        dfs.append(temp_df)

    # Concatenate *all* surviving DataFrames
    if not dfs:
        raise RuntimeError("No non‑empty BED files found – nothing to do.")

    modkit_bed_df = pd.concat(dfs, ignore_index=True)

    # Write concatenated BED
    out_path = Path(output_name)
    modkit_bed_df.to_csv(out_path, sep="\t", header=False, index=False)
    print(f"Combined BED written to {out_path} "
          f"({len(modkit_bed_df):,} rows)")

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
def prepare_chr_plotting_data(coverage_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add convenience columns used in downstream plots.
    Keeps both ‘condition’ *and* ‘exp_id’.
    """
    df = coverage_df.rename(columns={"mod_frac": "m6A_frac"}).copy()

    # simple chromosome group annotation (autosomes vs X, etc.)
    def _chr_group(chr_name):
        return "Autosome" if chr_name not in {"CHROMOSOME_X", "CHROMOSOME_V"} else chr_name.lstrip("CHROMOSOME_")
    df["chromosome_group"] = df["chromosome"].map(_chr_group)

    return df

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

import os
import io
import contextlib
import warnings

def process_region_for_summary(args):
    """
    Wrapper around get_summary_from_bam that captures all stdout/stderr
    from the modkit call into a buffer, then only prints it if the
    returned DataFrame is empty (so you can inspect what went wrong).
    """
    (chromosome, start, end, _, _, strand,
     each_bam, each_condition, each_thresh,
     each_exp_id, sampling_frac, bed_file) = args

    # 1) Run the summary call, capturing all output into `buffer`
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        temp_df = get_summary_from_bam(
            sampling_frac=sampling_frac,
            a_threshold=each_thresh,
            modkit_path="/Data1/software/modkit_v0.3/modkit",
            bam_path=each_bam,
            each_condition=each_condition,
            each_exp_id=each_exp_id,
            c_threshold=each_thresh,  # keep thresholds matched
            thread_ct=12,
            chromosome=chromosome,
            start=start,
            end=end,
        )

    # 2) If we got nothing back, dump the captured text for inspection
    if temp_df.empty:
        debug_txt = buffer.getvalue()
        print(f"[DEBUG modkit output for {chromosome}:{start}-{end}  "
              f"(cond={each_condition}  exp_id={each_exp_id})]:\n{debug_txt}")
        warnings.warn(
            f"Empty result for {chromosome}:{start}-{end}  "
            f"(cond={each_condition}  exp_id={each_exp_id})",
            RuntimeWarning
        )
        return temp_df

    # 3) Otherwise attach the genomic coords and return as before
    temp_df["chromosome"] = chromosome
    temp_df["start"]      = start
    temp_df["end"]        = end
    return temp_df

import shutil   # add at top if not already imported
import os
import pandas as pd
import multiprocessing
from tqdm.auto import tqdm
from tempfile import NamedTemporaryFile

def _process_bam_with_bed(args):
    """Wrapper when a single BED is given; no per‑region looping."""
    (bam, cond, thr, exp, sampling_frac, bed_path) = args
    return get_summary_from_bam(
        sampling_frac=sampling_frac,
        a_threshold=thr,
        modkit_path="/Data1/software/modkit_v0.3/modkit",
        bam_path=bam,
        each_condition=cond,
        each_exp_id=exp,
        c_threshold=thr,
        thread_ct=2,
        bed_file=bed_path,          # ← key line
    )

from io import StringIO
import pandas as pd

def extract_call_table(path_or_str):
    """
    Given the text of modkit --tsv output (or a file path), return just
    the modification‑calls table (columns: base code pass_count …).
    """
    # read whole text
    if os.path.exists(path_or_str):
        with open(path_or_str, "r") as fh:
            lines = fh.readlines()
    else:  # already a string
        lines = path_or_str.splitlines()

    # find the header row that starts the calls table
    for i, ln in enumerate(lines):
        if ln.startswith("base\t"):
            start = i
            break
    else:
        raise ValueError("Could not find 'base' header in the TSV text")

    table_text = "".join(lines[start:])          # join from header onward
    df_calls   = pd.read_csv(StringIO(table_text), sep="\t")

    return df_calls

def process_and_export_summary_by_region(
    sampling_frac,
    type_selected,
    new_bam_files,
    conditions,
    thresh_list,
    exp_ids,
    modkit_bed_df,
    force_replace=False,
    output_directory="temp_files",
    bed_file=None,                 # <── pass your BED path here
):
    """Streams summary TSVs to disk; now BED‑aware."""
    summary_table_name = os.path.join(
        output_directory,
        "_".join([conditions[0], conditions[-1],
                  str(sampling_frac), f"thresh{thresh_list[0]}",
                  type_selected[0]]) + "_summary_table.csv")

    if not force_replace and os.path.exists(summary_table_name):
        print("Summary table exists → loading from disk …")
        return pd.read_csv(summary_table_name, sep="\t")

    # ----------------------------------------------------------------
    # Build job list
    # ----------------------------------------------------------------
    if bed_file:  # ── only ONE job per BAM
        args_list = [
            (bam, cond, thr, exp, sampling_frac, bed_file)
            for bam, cond, thr, exp in zip(
                new_bam_files, conditions, thresh_list, exp_ids)
        ]
        worker   = _process_bam_with_bed
    else:        # ── original per‑region jobs
        args_list = [
            (row[0], row[1], row[2], row[3], row[4], row[5],
             bam, cond, thr, exp, sampling_frac, None)
            for _, row in modkit_bed_df.iterrows()
            for bam, cond, thr, exp in zip(
                new_bam_files, conditions, thresh_list, exp_ids)
        ]
        worker   = process_region_for_summary

    total_jobs = len(args_list)
    os.makedirs(output_directory, exist_ok=True)

    header_written = False
    with NamedTemporaryFile(mode="w", delete=False) as tmp, \
         multiprocessing.Pool(processes=32) as pool:

        for temp_df in tqdm(
                pool.imap_unordered(worker, args_list, chunksize=1),
                total=total_jobs, desc="Jobs", unit="job"):
            if temp_df.empty:
                continue
            temp_df.to_csv(tmp, sep="\t", header=not header_written, index=False)
            header_written = True

    if not header_written:
        raise ValueError("All summaries empty—nothing written.")

    shutil.move(tmp.name, summary_table_name)

    print(f"Exported summary → {summary_table_name}")
    return pd.read_csv(summary_table_name, sep="\t")

import pandas as pd
import numpy as np

def ensure_tidy_summary(df_in: pd.DataFrame) -> pd.DataFrame:
    """
    Convert Modkit 'summary' outputs that use a wide 'mod_bases' layout into
    a tidy schema with at least:
        ['condition','exp_id','base','code','pass_count', ... (region cols if present)]
    Passes through untouched if already tidy.
    """
    df = df_in.copy()

    # Already tidy?
    if {'base','code','pass_count'}.issubset(df.columns):
        return df

    # Expect wide layout with 'mod_bases' and a single numeric value column
    if 'mod_bases' not in df.columns:
        raise KeyError(
            "Expected either tidy columns {'base','code','pass_count'} "
            "or wide layout with 'mod_bases'. Got: "
            f"{df.columns.tolist()}"
        )

    # Prefer common numeric column names; else pick the first numeric column
    preferred = [c for c in ['C,A','value','count','pass_count','n'] if c in df.columns]
    if preferred:
        vcol = preferred[0]
    else:
        id_like = {'mod_bases','condition','exp_id','chromosome','start','end','region','name','strand'}
        numeric_candidates = [c for c in df.columns
                              if c not in id_like and pd.api.types.is_numeric_dtype(df[c])]
        if not numeric_candidates:
            raise ValueError("No numeric value column found in wide summary.")
        vcol = numeric_candidates[0]

    id_cols = [c for c in ['condition','exp_id','chromosome','start','end','region','name','strand']
               if c in df.columns]

    rows = []
    # Group by identifiers (condition/exp_id/region coords, if present)
    for _, sub in df.groupby(id_cols, dropna=False) if id_cols else [(None, df)]:
        def get(name: str) -> float:
            exact = sub.loc[sub['mod_bases'] == name, vcol]
            if not exact.empty:
                return float(exact.iloc[0])
            pref = sub.loc[sub['mod_bases'].str.startswith(name, na=False), vcol]
            return float(pref.iloc[0]) if not pref.empty else 0.0

        a_mod = get('A_pass_calls_modified_a')
        a_un  = get('A_pass_calls_unmodified')
        c_mod = get('C_pass_calls_modified_m')
        c_un  = get('C_pass_calls_unmodified')

        base_rows = [
            {'base':'A','code':'a','pass_count':a_mod},
            {'base':'A','code':'-','pass_count':a_un},
            {'base':'C','code':'m','pass_count':c_mod},
            {'base':'C','code':'-','pass_count':c_un},
        ]

        common = {k: (sub.iloc[0][k] if id_cols else None) for k in id_cols}
        for r in base_rows:
            rows.append({**common, **r})

    out = pd.DataFrame(rows)
    # Keep region coords if they exist; create placeholders if not
    # (create_coverage_df expects 'chromosome' and 'start' for per-region plots)
    for col in ['chromosome','start']:
        if col not in out.columns:
            out[col] = np.nan

    # Ensure types
    if out['pass_count'].dtype != float and out['pass_count'].dtype != int:
        out['pass_count'] = pd.to_numeric(out['pass_count'], errors='coerce').fillna(0)

    return out


import pandas as pd
from pathlib import Path
import pandas as pd
import numpy as np
from pathlib import Path

def create_coverage_df(code: str,
                       summary_bam_df: pd.DataFrame,
                       out_csv: str,
                       *,
                       save: bool = True) -> pd.DataFrame:
    """
    Build a per-region (or genome-wide if region coords are missing) coverage table.

    Returns columns:
      condition · exp_id · [chromosome · start] · mod_pass · canon_pass · mod_frac
    """
    # Normalize if still wide
    if not {'base','code','pass_count'}.issubset(summary_bam_df.columns):
        summary_bam_df = ensure_tidy_summary(summary_bam_df)

    base_filter = "A" if code == "a" else "C"
    mod_code, canon_code = code, "-"

    required = {'condition','exp_id','code','pass_count','base'}
    missing = required - set(summary_bam_df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # Decide grouping keys: include region keys only if they exist AND have any non-NaN
    region_keys = []
    for k in ("chromosome", "start"):
        if k in summary_bam_df.columns and summary_bam_df[k].notna().any():
            region_keys.append(k)
    group_keys = ["condition", "exp_id"] + region_keys

    subset = summary_bam_df.loc[
        (summary_bam_df["base"] == base_filter)
        & (summary_bam_df["code"].isin([mod_code, canon_code])),
        group_keys + ["code", "pass_count"]
    ]
    if subset.empty:
        raise ValueError(
            f"No rows for base={base_filter} with codes {mod_code!r} or '-' "
            "after filtering; check upstream summary."
        )

    # Sum pass_count per group × code; keep NaN groups if any
    agg = (
        subset.groupby(group_keys + ["code"], dropna=False, as_index=False)["pass_count"]
              .sum()
    )

    # Pivot and ensure both code columns exist
    pivot = agg.pivot_table(
        index=group_keys, columns="code", values="pass_count", fill_value=0
    )
    # If only one code is present, add the missing one as zeros
    for c in (mod_code, canon_code):
        if c not in pivot.columns:
            pivot[c] = 0

    pivot = (
        pivot.rename(columns={mod_code: "mod_pass", canon_code: "canon_pass"})
             .reset_index()
    )
    pivot.columns.name = None  # drop the 'code' axis name if present

    denom = (pivot["mod_pass"] + pivot["canon_pass"]).replace(0, np.nan)
    pivot["mod_frac"] = pivot["mod_pass"] / denom

    if save:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        pivot.to_csv(out_csv, sep="\t", index=False)

    return pivot

def prepare_chr_plotting_data(coverage_df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare per-region coverage for plotting.
    Safe if 'chromosome' is missing by checking alternatives.
    """
    df = coverage_df.rename(columns={"mod_frac": "m6A_frac"}).copy()

    # tolerate alt column names
    if "chromosome" not in df.columns:
        if "chrom" in df.columns:
            df = df.rename(columns={"chrom": "chromosome"})
        else:
            raise KeyError("Expected 'chromosome' in coverage_df.")

    # annotate X vs autosomes
    def _grp(c):
        return "Autosome" if c not in {"CHROMOSOME_X", "CHROMOSOME_V"} else c.replace("CHROMOSOME_","")
    df["chromosome_group"] = df["chromosome"].map(_grp)

    return df
