#######
# Author: Yuri Malina
# Contact: ymalina@gmail.com
# Purpose: R scripts for generating tabix file use by nanomethviz.
#######

### Instructions for launching rstudio
# Launching docker with Data1 bind mounted.
# docker run --rm -ti -v /Data1:/home/rstudio/Data1 -e PASSWORD=jupyter -e ROOT=true -p 8787:8787 rocker/rstudiov2
# Installed on ubuntu 18 server using https://rocker-project.org/images/versioned/rstudio.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("NanoMethViz", quietly = TRUE))
  BiocManager::install("NanoMethViz")

library(NanoMethViz)

### Nanomethviz Vignettes
#vignette("Introduction", package = "NanoMethViz")
#vignette("ImportingExportingData", package = "NanoMethViz")

### CONFIGURATIONS
samples_for_tabix <- c("N2")
basecalls_path <- ("~/Data1/seq_data/TubeD1a_N2_Fiberseq_Hia5_MSssI_12_22_22/basecalls/m6A/")
methy_tabix_combined <- "~/Data1/seq_data/nanomethviz_files/TubeD1a_N2_2uM-Hia5_fiber-seq_12_22_22.bgz"

### If tabix files does not exist, create it:
# you should see messages when running this yourself
if (file.exists(methy_tabix_combined) == FALSE){
  create_tabix_file(
    c(basecalls_path,"per_read_modified_base_calls.txt"),
      methy_tabix_combined, samples_for_tabix)
}




