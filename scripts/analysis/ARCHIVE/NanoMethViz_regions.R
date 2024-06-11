#######
# Author: Yuri Malina
# Contact: ymalina@gmail.com
# Purpose: R scripts for visualizing methylation levels using
# Nanomethviz to calculate % methylation across regions in a bed file
# and ggplot to output box plots.
#######

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("NanoMethViz", quietly = TRUE))
  BiocManager::install("NanoMethViz")

library(NanoMethViz)
#vignette("ImportingExportingData", package = "NanoMethViz")
#vignette("Introduction", package = "NanoMethViz")


### Function to load sample
load_sample<- function(x,y,z){
  sample = c(x)
  group = c(x)
  sample_anno = data.frame(sample,group,stringsAsFactors = FALSE)
  nano_result <<- NanoMethResult(z,sample_anno,exons = NULL)
  region <<- read.table(y,sep="\t",header=TRUE)
}

### Block, can be copied for multiple samples
print("Starting First Sample")
bgz_file = "~/Data1/seq_data/nanomethviz_files/Tube4_N2_2uM-Hia5_fiber-seq_11_21_22.bgz"
sample_name ="N2"
region_name ="~/Data1/reference/ce11-100kb_windows_labels_nano.bed"
output_name = "~/Data1/seq_data/nanomethviz_files/N2_Mixed_Stage_Hia5_2uM_t0p8_100k_chr.tsv"
meth_thresh = 0.8 #default of megalodon is 0.8

load_sample(sample_name,region_name,bgz_file)
query_results <- region_methy_stats(nano_result,region,threshold=meth_thresh)
#query_results <- query_methy(nano_result, region$chr, region$start, region$end, simplify = TRUE, force = FALSE)
write.table(query_results,output_name, sep="\t",quote = F)

sample_prev <- read.table(output_name, sep="\t", header=TRUE)

ggplot(data = sample_prev, aes(x=chromosome, y=prevalence)) + 
  geom_boxplot(aes(fill=chromosome)) + ggtitle(c("N2 Mixed Stage Embryos + Hia5 2uM")) + 
  theme(text = element_text(size = 20)) + 
  geom_text(data = sdc2_means, aes(label = prevalence, y = prevalence))

### END OF BLOCK ###

#ama1_n2_prev <- read.table("~/Data1/seq_data/nanomethviz_files/AMA1_N2_10kb_m6a.tsv", sep="\t", header=TRUE)
#
#ggplot(data = ama1_n2_prev, aes(x=chromosome, y=prevalence)) + 
#  geom_boxplot(aes(fill=target)) +
#  ggtitle("DiMeLo-seq target: AMA-1 versus N2") +
#  facet_wrap(~target) + 
#  theme(text = element_text(size = 20))  +
#  stat_summary(fun=mean, geom="point", shape=10, size=3,
#              position = position_dodge2(width = 0.1,   
#                                         preserve = "single")) 


