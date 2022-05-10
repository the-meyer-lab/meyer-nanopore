import dimelo as dm
dm.plot_enrichment_profile("/home/ubuntu/Data2/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", 'N2MixedStageEmbryos5', "/home/ubuntu/Data2/reference/ce11-all-rex-10kb.bed", "A", "/home/ubuntu/Data2/viz/barcode01_m6A/", windowSize=10000, dotsize=0.5, smooth=50, min_periods=20, cores=8)
#dm.qc_report("/home/ubuntu/Data2/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam", "Mixed Stage Embryos N2 XX", "/home/ubuntu/Data2/seq_data/210614_Raja/megalodon/barcode01_CpG/dimelo", cores=4)
#dm.plot_enrichment("/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", ["chrI","chrII","chrIII","chrIV","chrV","chrX"], "/Data1/reference/ce11-chrm-regions.bed", "A", "/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/dimelo")
#dm.plot_enrichment("/home/ubuntu/Data2/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", "chrX", "/home/ubuntu/Data2/reference/ce11-10kb_test.bed", "A", "/home/ubuntu/Data2/viz/barcode01_m6A/", cores=8)

