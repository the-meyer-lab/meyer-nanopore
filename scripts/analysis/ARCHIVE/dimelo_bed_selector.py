# NOTE: Requires "conda activate dimelo"

import dimelo as dm
import pandas as pandas
import numpy as np

sample_source = "chromosome" # "chr_type" or "type" or "chromosome"
sampleName = ["CHROMOSOME_I", "CHROMOSOME_II", "CHROMOSOME_III", "CHROMOSOME_IV", "CHROMOSOME_V","CHROMOSOME_X"] # "TES_q1" "strong_rex" "weak_rex" "type", "X", "Autosome"; Must be same number of unique values in selected bed rows.
chr_type_selected = ["X","Autosome"] # 'X' or "Autosome"
type_selected = ["whole_chr"] #TES_q1-4 | #TSS_q1-4 | strong/weak rex | whole_chr
max_regions = 0 # max regions to consider.0 = full set;
chromosome_selected = ["CHROMOSOME_I", "CHROMOSOME_II", "CHROMOSOME_III", "CHROMOSOME_IV", "CHROMOSOME_V","CHROMOSOME_X"]
strand_selected = ["+","-"] #+ and/or -
bed = "/Data1/reference/tss_tes_rex_combined.bed"
full_bed = pandas.read_csv(bed,sep='\t')
mods = "A" # {A,CG,A+CG}
### Tube 4
#bam1 = "/Data1/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/basecalls/mod_mappings.sorted.m6Aonly.bam"
#outDir = "/Data1/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_1121_22/viz/whole_chr"
### Tube D
bam1 = "/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/basecalls/m6A/mod_mappings.sorted.bam"
outDir = "/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/viz/whole_chr"
### Tube H
#bam1 = "/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/basecalls/m6A/mod_mappings.sorted.bam"
#outDir = "/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/viz/whole_chr"


m6A_thresh = 129 #default is 129
mC_thresh = 129 #default is 129
window = 1000 # size of window to use
coreNum = 96 # cores to use

bed=[]

if sample_source == "chr_type":
  selection = chr_type_selected
if sample_source == "type":
  selection = type_selected
if sample_source == "chromosome":
  selection = chromosome_selected

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
    temp_bed = full_bed[full_bed["chromosome"].str.contains(each_type) &
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

  drop_indices = np.random.choice(temp_bed.index, drop_count, replace=False)
  temp_bed.drop(drop_indices,inplace=True)
  temp_bed.reset_index(drop=True, inplace=True)
  temp_bed["start"]=temp_bed["start"]
  temp_bed["end"]=temp_bed["end"]
  #print(temp_bed)
  temp_bedfile = "/Data1/reference/temp_do_not_use_"+each_type+".bed"
  temp_bed.to_csv(temp_bedfile, sep="\t",header=False,index=False)

  # For first iteration
  if bed == []:
    bed = [temp_bedfile]

  # Otherwise append region to temporary bed file.
  else:
    bed.append(temp_bedfile)

print("BAM:",bam1)
print("BED: ",bed)

for eachBed, eachSample in zip(bed,sampleName):
  dm.parse_bam(bam1,eachSample,outDir,eachBed,mods,center=True,windowSize=window,threshA=m6A_thresh,cores=coreNum)

#parse_bam(
#          fileName: str,
#          sampleName: str,
#          outDir: str,
#          bedFile: str = None,
#          basemod: str = DEFAULT_BASEMOD,
#          center: bool = False,
#          windowSize: int = DEFAULT_WINDOW_SIZE,
#          region: str = None,
#          threshA: int = DEFAULT_THRESH_A,
#          threshC: int = DEFAULT_THRESH_C,
#          extractAllBases: bool = False,
#          cores: int = None,
#  )
#### Run dimelo plot enrichment file on subselected bed
#dm.plot_enrichment_profile(bam1, sampleName, bed, mods, outDir, windowSize=1000, dotsize=0.05, cores=coreNum,colors=["#053C5E","#BB4430"],threshA = m6A_thresh,threshC = mC_thresh)

#dm.plot_enrichment_profile(
#    "/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam",
#    'N2MixedStageEmbryos_tss_10k', "/Data1/reference/ce11-all-tss-10kb.bed", "A",
#    "/Data1/viz/barcode01_m6A/", windowSize=10000, dotsize=0.5, smooth=50, min_periods=20, cores=64)
#dm.parse_bam(
#    "/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam",
#    'AMA-1-Tube31D', '/Data1/seq_data/AMA1-31D-10-08-22/basecalls/', bedFile= "/Data1/reference/ce11-genebodies-2.bed", basemod="A",
#    center=True, extractAllBases=True, windowSize=1000
#)

#dm.plot_browser(bam, sampleName, "CHROMOSOME_X:10000000-10500000", mods, outDir, threshA=153, static=False, smooth=100, min_periods=10)

#dm.plot_enrichment_profile(bam, sampleName, bed, mods, outDir, windowSize=2000, dotsize=0.05)
#dm.plot_enrichment_profile(bam, "Rex", "/Data1/reference/rex1kb.bed", mods, outDir, windowSize=2000, dotsize=0.05, cores=96)

#sampleName = ["TES_Q1","TES_Q2","TES_Q3","TES_Q4"]
#bed = "/Data1/reference/Kreusi_allTSS_WS235_EMBRYO.bed"
#bed = ["/Data1/reference/dimelo_tes_q1.bed","/Data1/reference/dimelo_tes_q2.bed","/Data1/reference/dimelo_tes_q3.bed","/Data1/reference/dimelo_tes_q4.bed"]


#dm.plot_enrichment("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", "AMA1_Tube31D_0p1fix", "/Data1/reference/ce11-TSS.bed", "A", "/Data1/seq_data/AMA1-31D-10-08-22/viz")
# BE CAREFUL NO HYPHENS ALLOWED IN SAMPLE NAME!!!!!!!!!!

#dm.qc_report("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", 'AMA1-Tube31D-0p1fix', "/Data1/seq_data/AMA1-31D-10-08-22/viz/", cores=4)
#dm.plot_enrichment(bam, ["chrI","chrII","chrIII","chrIV","chrV","chrX"], "/Data1/reference/ce11-chrm-regions_labeled_dimelo.bed", "A", outDir)
#dm.plot_enrichment("/Data1/seq_data/N2_Yeast_Raja_06-14-21/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", "chrX", "/Data1/reference/ce11-10kb_test.bed", "A", "/Data1/viz/barcode01_m6A/", cores=8)
#dm.plot_enrichment("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", "chrX", "/Data1/reference/ce11-10kb_test.bed", "A", "/Data1/seq_data/AMA1-31D-10-08-22/viz", cores=8)

'''
### Commans to generate plot browser from Dimelo
# Remember to execute
# conda activate dimelo
import dimelo as dm
#bam = "/Data1/seq_data/N2_Yeast_Raja_06-14-21/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam"
sampleName = ["N2_Mixed_Stage_2uM_Hia5"]#,"TSS_Q2","TSS_Q3","TSS_Q4"]
#bed = "/Data1/reference/Kreusi_allTSS_WS235_EMBRYO.bed"
#bed = ["/Data2/reference/dimelo_tss_q1.bed"]#,"/Data1/reference/dimelo_tss_q2.bed","/Data1/reference/dimelo_tss_q3.bed","/Data1/reference/dimelo_tss_q4.bed"]
mods = "A"
#outDir = "/Data1/seq_data/N2_Yeast_Raja_06-14-21/megalodon/barcode01_m6A/viz"

bam1 = ["/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/basecalls/m6A/mod_mappings.commonsorted.bam"]
outDir = "/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/viz"
#bam1 = ["/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/basecalls/m6A/mod_mappings.sorted.bam"]
#outDir = "/Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/viz"

#sampleName = 'AMA1_31D_0p1fix_TSS'
#bed = "/Data1/reference/Kreusi_allTSS_WS235_EMBRYO.bed"
#mods = "A"
### Plot Browser
dm.plot_browser(
   bam1,
   sampleName,
   "CHROMOSOME_X:7180757-7182759",
   mods,
   outDir,
)

#Athresh = 204
#dm.plot_enrichment_profile(
#    "/Data2/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam",
#    'N2MixedStageEmbryos_tss_10k', "/Data2/reference/ce11-all-tss-10kb.bed", "A",
#    "/Data2/viz/barcode01_m6A/", windowSize=10000, dotsize=0.5, smooth=50, min_periods=20, cores=64)
#dm.parse_bam(
#    "/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam",
#    'AMA-1-Tube31D', '/Data1/seq_data/AMA1-31D-10-08-22/basecalls/', bedFile= "/Data1/reference/ce11-genebodies-2.bed", basemod="A",
#    center=True, extractAllBases=True, windowSize=1000
#)

#dm.plot_browser(bam, sampleName, "CHROMOSOME_X:10000000-10500000", mods, outDir, threshA=153, static=False, smooth=100, min_periods=10)

#dm.plot_enrichment_profile(bam, sampleName, bed, mods, outDir, windowSize=2000, dotsize=0.05)
#dm.plot_enrichment_profile(bam, sampleName, bed, mods, outDir, windowSize=2000, dotsize=0.05, cores=96)

#sampleName = ["TES_Q1"]#,"TES_Q2","TES_Q3","TES_Q4"]
#bed = "/Data1/reference/Kreusi_allTSS_WS235_EMBRYO.bed"
#bed = ["/Data2/reference/dimelo_tes_q1.bed","/Data2/reference/dimelo_tes_q2.bed","/Data2/reference/dimelo_tes_q3.bed","/Data2/reference/dimelo_tes_q4.bed"]
#dm.plot_enrichment_profile(bam, sampleName, bed, mods, outDir, windowSize=2000, dotsize=0.05,cores=96)

#dm.plot_enrichment("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", "AMA1_Tube31D_0p1fix", "/Data1/reference/ce11-TSS.bed", "A", "/Data1/seq_data/AMA1-31D-10-08-22/viz")
# BE CAREFUL NO HYPHENS ALLOWED IN SAMPLE NAME!!!!!!!!!!

#dm.qc_report("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", 'AMA1-Tube31D-0p1fix', "/Data1/seq_data/AMA1-31D-10-08-22/viz/", cores=4)
#dm.plot_enrichment("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", ["chrI","chrII","chrIII","chrIV","chrV","chrX"], "/Data1/reference/ce11-chrm-regions.bed", "A", "/Data1/seq_data/AMA1-31D-10-08-22/viz")
#dm.plot_enrichment("/Data1/seq_data/N2_Yeast_Raja_06-14-21/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", "chrX", "/Data1/reference/ce11-10kb_test.bed", "A", "/Data1/viz/barcode01_m6A/", cores=8)
#dm.plot_enrichment("/Data1/seq_data/AMA1-31D-10-08-22/basecalls/mod_mappings.sorted.m6Aonly.bam", "chrX", "/Data1/reference/ce11-10kb_test.bed", "A", "/Data1/seq_data/AMA1-31D-10-08-22/viz", cores=8)
'''

