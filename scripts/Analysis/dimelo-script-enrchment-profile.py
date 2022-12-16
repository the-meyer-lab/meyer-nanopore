import dimelo as dm
import pandas as pandas
import numpy as np


sampleName = ["X","Autosome"] #"TES_q1" "strong_rex" "weak_rex" from bed column "type"
type_selected = ["TES_q1","TES_q2","TES_q3","TES_q4"]
strand_selected = ["+"] #+ and/or -
bed = "/Data1/reference/tss_tes_rex_combined.bed"
full_bed = pandas.read_csv(bed,sep='\t')
mods = "A"
bam = "/Data1/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/basecalls/mod_mappings.sorted.m6Aonly.bam"
outDir = "/Data1/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/viz"
m6A_thresh = 129 #default is 129

bed=[]
for each_type in sampleName:
  # REGION CONFIGURATION
  max_regions = 200 #1 = full set; 0.1 = 90% dropped at random.
  chromosome_selected = ["CHROMOSOME_X","CHROMOSOME_I","CHROMOSOME_II","CHROMOSOME_III","CHROMOSOME_IV","CHROMOSOME_V"]
  temp_bed = full_bed[full_bed["chromosome"].isin(chromosome_selected) &
                      full_bed["chr-type"].str.contains(each_type) &
                      full_bed["type"].isin(type_selected) &
                      full_bed["strand"].isin(strand_selected)]

  drop_count = len(temp_bed)-max_regions
  if(drop_count<0):
    drop_count=0

  drop_indices = np.random.choice(temp_bed.index, drop_count, replace=False)
  temp_bed.drop(drop_indices,inplace=True)
  temp_bed.reset_index(drop=True, inplace=True)
  temp_bed["start"]=temp_bed["start"]-1
  temp_bed["end"]=temp_bed["end"]+1
  print(temp_bed)

  temp_bedfile = "/Data1/reference/temp_do_not_use_"+each_type+".bed"
  temp_bed.to_csv(temp_bedfile, sep="\t",header=False,index=False)
  if bed == []:
    bed = [temp_bedfile]
  else:
    bed.append(temp_bedfile)

dm.plot_enrichment_profile(bam, sampleName, bed, mods, outDir, windowSize=1000, dotsize=0.05, cores=96,colors=["#053C5E","#BB4430"],threshA = m6A_thresh)

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

