import pysam
from joblib import Parallel, delayed

from typing import List, Tuple, Union
import pandas as pd
import multiprocessing
import numpy as np
import os.path

from dimelo.utils import  create_sql_table, execute_sql_command


class Region(object):
    def __init__(self, region: Union[str, pd.Series]):
        """Represents a region of genetic data.
        Attributes:
                - chromosome: string name of the chromosome to which the region applies
                - begin: integer start position of region
                - end: integer end position of region
                - size: length of region
                - string: string representation of region
                - strand: string specifying forward or reverse strand; either "+" or "-" (default +)
        TODO:
                - There should be no reason the second datatype can't be something like typing.Sequence, but technically a pd.Series isn't a valid Sequence. Not sure what the best way to annotate this would be.
                - Can I re-enable errors being raised from this method?
                - Consider changing size to be a property rather than an attribute
                - Change string to be the magic string method
                - Determine why the string representation exists, and why it is important
        """
        self.chromosome = None
        self.begin = None
        self.end = None
        self.size = None
        self.string = None
        self.strand = "+"

        if isinstance(region, str):  # ":" in region:
            # String of format "{CHROMOSOME}:{START}-{END}"
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                try:
                    # see if just integer chromosomes are used
                    self.chromosome = int(self.chromosome)
                except ValueError:
                    pass
                self.begin, self.end = [int(i) for i in interval.split("-")]
            except ValueError:
                raise TypeError(
                    "Invalid region string. Example of accepted format: 'chr5:150200605-150423790'"
                )
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        elif isinstance(region, pd.Series):
            # Ordered sequence containing [CHROMOSOME, START, END] and optionally [STRAND], where STRAND can be either "+" or "-"
            self.chromosome = region[0]
            self.begin = region[1]
            self.end = region[2]
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
            # strand of motif to orient single molecules
            # if not passed just keep as all +
            if len(region) >= 4:
                if (region[3] == "+") or (region[3] == "-"):
                    self.strand = region[3]
                # handle case of bed file with additional field that isn't strand +/-
                else:
                    self.strand = "+"
            else:
                self.strand = "+"
        else:
            raise TypeError(
                "Unknown datatype passed for Region initialization"
            )

def reads_window(
    fileName: str,
    basemod: str,
    window: Region,
    threshA: int = 129,
    extractAllBases: bool = False,
) -> None:
    window_mods=0
    window_bps=0
    bam = pysam.AlignmentFile(fileName, "rb")
    for read in bam.fetch(
        reference=window.chromosome, start=window.begin, end=window.end
        ):
            read_mods, read_bps = get_modified_reference_positions(
                read,
                basemod,
                window,
                threshA,
                fileName,
                extractAllBases,
            )

            window_mods += read_mods
            window_bps += read_bps

    if window_bps != 0:
        fracs = float(window_mods)/float(window_bps)
    else:
        fracs = 0


    #print(df["Start_bp"].to_list())

    table_name = "methylationAggregate_XChrom_500bin"

    command = (
            """INSERT OR IGNORE INTO """
            + table_name
            + """ VALUES(?,?,?,?);"""
    )

    data_fill = (window.begin, window_mods, window_bps, fracs)
    execute_sql_command(command, "/Data2/seq_data/500bp-bins.db", data_fill)

    print("The fraction of methylated bases in the window beggining at bp " + str(window.begin) + " is " + str(fracs))

#    if len(df["Start_bp"].to_list()) % 500 == 0:
#        print("saving")
#        if os.path.exists('/Data2/seq_data/Dataframe_frac-500bp-bins_complete.pkl'):
#            df_complete = pd.read_pickle('/Data2/seq_data/Dataframe_frac-500bp-bins_complete.pkl')
#            df_complete = pd.concat([df_complete, df], ignore_index=True)
#            df_complete.to_pickle('/Data2/seq_data/Dataframe_frac-500bp-bins_complete.pkl')

#        else:
#            df.to_pickle('/Data2/seq_data/Dataframe_frac-500bp-bins_complete.pkl')
#        df = pd.DataFrame(columns=['Start_bp', 'methylated_bases', 'total_bases', 'frac'])
#        df.to_pickle('/Data2/seq_data/Dataframe_frac-500bp-bins.pkl')
#
#    else:
#        df.to_pickle('/Data2/seq_data/Dataframe_frac-500bp-bins.pkl')






def get_modified_reference_positions(
    read: pysam.AlignedSegment,
    basemod: str,
    window: Region,
    threshA: int,
    fileName: str,
    extractAllBases: bool,
):
    """Extract mA and mC pos & prob information for the read
    Args:
            :param read: single read from bam file
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param window: window from bed file
            :param threshA: threshold above which to call an A base methylated
    Return:
        For each mod, you get the positions where those mods are and the probabilities for those mods (parallel vectors)
    TODO:
            - This is referenced in plot_joint_enrichment
            - Oh boy, what does this return? At minimum it should have a type annotation. At maximum, it should have a class.
              - See get_mod_reference_positions_by_mod() for details
    """
    read_num_mod = 0
    read_num_bp = 0
    if (read.has_tag("Mm")) & (";" in read.get_tag("Mm")):
        mod1 = read.get_tag("Mm").split(";")[0].split(",", 1)[0]
        mod2 = read.get_tag("Mm").split(";")[1].split(",", 1)[0]
        # mod1_list = read.get_tag("Mm").split(";")[0].split(",", 1)
        # mod2_list = read.get_tag("Mm").split(";")[1].split(",", 1)
        base = basemod[0]  # this will be A, C, or A
        if basemod == "A+CG":
            base2 = basemod[2]  # this will be C for A+CG case
        else:  # in the case of a single mod will just be checking that single base
            base2 = base
        # if len(mod1_list) > 1 and (base in mod1 or base2 in mod1):
        if base in mod1 or base2 in mod1:
            mod1, bp1 = get_mod_reference_positions_by_mod(
                read,
                mod1,
                0,
                window,
                threshA,
                fileName,
                extractAllBases,
            )
            read_num_mod += mod1
            read_num_bp += bp1

        # if len(mod2_list) > 1 and (base in mod2 or base2 in mod2):
        if base in mod2 or base2 in mod2:
            mod2, bp2  = get_mod_reference_positions_by_mod(
                read,
                mod2,
                1,
                window,
                threshA,

                fileName,
                extractAllBases,
            )
            read_num_mod += mod2
            read_num_bp += bp2
        return (read_num_mod, read_num_bp)



def get_mod_reference_positions_by_mod(
    read: pysam.AlignedSegment,
    basemod: str,
    index: int,
    window: Region,
    threshA: int,
    fileName: str,
    extractAllBases: bool,
):
    """Get positions and probabilities of modified bases for a single read
    Args:
            :param read: one read in bam file
            :param mod: which basemod, reported as base+x/y/m
            :param window: window from bed file
            :param center: report positions with respect to reference center (+/- window size) if True or in original reference space if False
            :param threshA: threshold above which to call an A base methylated
            :param threshC: threshold above which to call a C base methylated
            :param windowSize: window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only applicable for center=True
            :param index: 0 or 1

    TODO:
            - What is index? What is its type? What does it represent?
            - Oh boy, what does this return? At minimum it should have a type annotation. At maximum, it should have a class.
                -  -> Tuple(str, List[int], List[int])
    """
    modded = 0
    allbp = 0
    modsPresent = True
    base, mod = basemod.split("+")
    num_base = len(read.get_tag("Mm").split(";")[index].split(",")) - 1
    # get base_index
    base_index = np.array(
        [
            i
            for i, letter in enumerate(read.get_forward_sequence())
            if letter == base
        ]
    )
    # get reference positons
    refpos = np.array(read.get_reference_positions(full_length=True))
    if read.is_reverse:
        refpos = np.flipud(refpos)
    modified_bases = []
    if num_base == 0:
        modsPresent = False
    if modsPresent:
        deltas = [
            int(i) for i in read.get_tag("Mm").split(";")[index].split(",")[1:]
        ]
        Ml = read.get_tag("Ml")
        if index == 0:
            probabilities = np.array(Ml[0:num_base], dtype=int)
        if index == 1:
            probabilities = np.array(Ml[0 - num_base :], dtype=int)
        # determine locations of the modified bases, where index_adj is the adjustment of the base_index
        # based on the cumulative sum of the deltas
        locations = np.cumsum(deltas)
        # loop through locations and increment index_adj by the difference between the next location and current one + 1
        # if the difference is zero, therefore, the index adjustment just gets incremented by one because no base should be skipped
        index_adj = []
        index_adj.append(locations[0])
        i = 0
        for i in range(len(locations) - 1):
            diff = locations[i + 1] - locations[i]
            index_adj.append(index_adj[i] + diff + 1)
        # get the indices of the modified bases
        modified_bases = base_index[index_adj]

    # extract CpG sites only rather than all mC
    keep = []
    prob_keep = []
    all_bases_index = []
    probs = []
    i = 0
    seq = read.get_forward_sequence()
    # deal with None for refpos from soft clipped / unaligned bases

    for b in base_index:
        if refpos[b] is not None:
            all_bases_index.append(
                b
            )  # add to all_bases_index whether or not modified
            if b in modified_bases:
                if probabilities[i] >= threshA:
                    keep.append(b)
                    prob_keep.append(i)
            if extractAllBases:
                if b in modified_bases:
                    probs.append(probabilities[i])
                else:
                    probs.append(0)
        # increment for each instance of modified base
        if b in modified_bases:
            i = i + 1

    for pos in refpos[all_bases_index]:
        if pos >= window.begin and pos <= window.end:
            if pos in refpos[keep]:
                modded += 1
                allbp += 1
            else:
                allbp += 1

    return (modded, allbp)


cores = 16
bedFile = "/Data2/reference/bins_X_500bp.bed"


cores_avail = multiprocessing.cpu_count()
if cores is None:
    num_cores = cores_avail
else:
    # if more than available cores is specified, process with available cores
    if cores > cores_avail:
        num_cores = cores_avail
    else:
        num_cores = cores

bed = pd.read_csv(bedFile, sep="\t", header=None)
windows = []
for _, row in bed.iterrows():
    windows.append(Region(row))

table_name = "methylationAggregate_XChrom_500bin"
cols = ['Start_bp', 'methylated_bases', 'total_bases', 'frac']
dtypes = ["INT", "INT", "INT", "FlOAT"]
create_sql_table("/Data2/seq_data/500bp-bins.db", table_name, cols, dtypes)


Parallel(n_jobs=num_cores)(
        delayed(reads_window )(
            "/Data2/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam",
            "A",
            window,
             )
        for window in windows
    )

print("finished")
