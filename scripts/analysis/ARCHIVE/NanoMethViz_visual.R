#######
# Author: Yuri Malina
# Contact: ymalina@gmail.com
# Purpose: R scripts for visualizing methylation levels using
# Nanomethviz to plot methyltion across aligned regions, e.g. TSS and TES
#######


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("NanoMethViz", quietly = TRUE))
  BiocManager::install("NanoMethViz")

if (!require("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")

library(NanoMethViz)
library(dplyr)
library(GenomicRanges)

### Import c. elegans gene annotations
#Source: vignette("ExonAnnotations", package = "NanoMethViz")
exon_tibble <- rtracklayer::import(system.file(package = "NanoMethViz", "c_elegans.gtf.gz"))
exon_tibble <- exon_tibble %>%
  as.data.frame() %>%
  dplyr::rename(
    chr = "seqnames",
    symbol = "gene_name"
  ) %>%
  dplyr::select("gene_id", "chr", "strand", "start", "end", "transcript_id", "symbol")
head(exon_tibble)
###

#vignette("ImportingExportingData", package = "NanoMethViz")
#vignette("Introduction", package = "NanoMethViz")

### CONFIGURATION
sample = c("N2")
group = c("N2 Mixed Stage Embryos + Hia5_group")
sample_anno = data.frame(sample,group,stringsAsFactors = FALSE)
methy = "~/Data1/seq_data/nanomethviz_files/Tube4_N2_2uM-Hia5_fiber-seq_11_21_22.bgz"
nano_result = NanoMethResult(methy,sample_anno,exons = NULL)
master_region <- read.table("~/Data1/reference/kreusi_tss_corrected_allq.bed",sep="\t",header=TRUE)
flank_set = 1000 #base pairs on either side of the plot
span_set = 0.01 #smoothing on plot
meth_thresh = 0.5 # 0.8 default for methylation cutoff

# REGION CONFIGURATION
max_regions = 100 #1 = full set; 0.1 = 90% dropped at random.
quartile_selected = c("q1","q2") # "q1","q2","q3","q4" are all options
chromosome_selected = c("CHROMOSOME_X","CHROMOSOME_I","CHROMOSOME_II","CHROMOSOME_III","CHROMOSOME_IV","CHROMOSOME_V")
region_selection = master_region[master_region$chrom %in% chromosome_selected &
                                   master_region$N2_RNAi_quartile %in% quartile_selected,]
drop_count = nrow(region_selection)-max_regions
if(drop_count<0){
  drop_count=0}
region_selection = region_selection[-sample(1:nrow(region_selection),drop_count),]
 

agg_plot <- plot_agg_regions(
  nano_result,
  region_selection,
  flank = flank_set,
  stranded = TRUE,
  span = span_set,
  group_col = "N2_chromosome_type",
  binary_threshold = meth_thresh
)

agg_plot + ylim(0.3,0.6) + theme(text = element_text(size = 20))  + ggtitle(paste(group,"Meth Threshold = ",meth_thresh,max_regions,"of",quartile_selected,"regions"))  

### REX 19
annotations_frame= data.frame(
  "chr"=c("CHROMOSOME_X"),
  "start"=c(1492414),
  "end"=c(1492434))

plot_grange(
  nano_result, 
  GenomicRanges::GRanges("CHROMOSOME_X:1491424-1493424"),
  binary_threshold = meth_thresh,
  heatmap=TRUE,
  spaghetti = TRUE,
  span = 0.075,
  anno_regions=annotations_frame)

### REX 32
annotations_frame= data.frame(
  "chr"=c("CHROMOSOME_X"),
  "start"=c(2997067),
  "end"=c(2997087))

plot_grange(
  nano_result, 
  GenomicRanges::GRanges("CHROMOSOME_X:2996077-2998077"),
  binary_threshold = meth_thresh,
  heatmap=TRUE,
  spaghetti = TRUE,
  span = 0.075,
  anno_regions=annotations_frame)

### REX 45
annotations_frame= data.frame(
  "chr"=c("CHROMOSOME_X"),
  "start"=c(8036353),
  "end"=c(8036373))

plot_grange(
  nano_result, 
  GenomicRanges::GRanges("CHROMOSOME_X:8035363-8037363"),
  binary_threshold = meth_thresh,
  heatmap=TRUE,
  spaghetti = TRUE,
  span = 0.075,
  anno_regions=annotations_frame)