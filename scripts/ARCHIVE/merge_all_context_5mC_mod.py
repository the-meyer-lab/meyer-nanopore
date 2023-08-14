import pysam

# Open the .bam files using pysam
all_context_file = pysam.AlignmentFile("/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/basecalls/m6A/mod_mappings.sorted.bam",
                                "rb")
mC_file = pysam.AlignmentFile("/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/basecalls/5mC/mod_mappings.sorted.bam",
                                "rb")

count=0

for all_context_read, mC_read in zip(all_context_file.fetch(until_eof=True),mC_file.fetch(until_eof=True)):
    print("ALL CONTEXT BAM")
    print(pysam.all_context_read.query_name(all_context_read))
    print(all_context_read.get_tags())
    print("5mC BAM")
    print(mC_read.query_name())
    print(mC_read.get_tags())
    count += 1
    if count >= 1:
        break


# BAM FILE FORMAT
# ReadID 0  chr 1   40  3184M   *   start   end seq

# Close the .bam files
all_context_file.close()
mC_file.close()

