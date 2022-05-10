import dimelo as dm
#dm.qc_report("/Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam", "Mixed Stage Embryos N2 XX", "/Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/dimelo", cores=12)
#dm.plot_enrichment("/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", ["chrI","chrII","chrIII","chrIV","chrV","chrX"], "/Data1/reference/ce11-chrm-regions.bed", "A", "/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/dimelo")
dm.plot_enrichment("/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam", "chrX", "/Data1/reference/ce11-chrm-regions_Xonly.bed", "A", "/Data1/viz/barcode01_m6A/", cores=2)
