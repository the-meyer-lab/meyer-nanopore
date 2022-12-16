# Remember to execute
# conda activate dimelo
import dimelo as dm
#bam = "/Data1/seq_data/N2_Yeast_Raja_06-14-21/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam"
sampleName = ["N2_Mixed_Stage_2uM_Hia5"]#,"TSS_Q2","TSS_Q3","TSS_Q4"]
#bed = "/Data1/reference/Kreusi_allTSS_WS235_EMBRYO.bed"
#bed = ["/Data2/reference/dimelo_tss_q1.bed"]#,"/Data1/reference/dimelo_tss_q2.bed","/Data1/reference/dimelo_tss_q3.bed","/Data1/reference/dimelo_tss_q4.bed"]
mods = "A"
#outDir = "/Data1/seq_data/N2_Yeast_Raja_06-14-21/megalodon/barcode01_m6A/viz"

bam = "/Data2/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/basecalls/mod_mappings.sorted.m6Aonly.bam"
#sampleName = 'AMA1_31D_0p1fix_TSS'
#bed = "/Data1/reference/Kreusi_allTSS_WS235_EMBRYO.bed"
#mods = "A"
outDir = "/Data2/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/viz"

### Plot Browser
dm.plot_browser(
   bam,
   sampleName,
   "CHROMOSOME_X:998400-1000400",
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

