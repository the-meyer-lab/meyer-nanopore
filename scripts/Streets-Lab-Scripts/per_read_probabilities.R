#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-04-06
#########################

#######################################################################################
# the per_read_base_calls.txt file itself is too large for RStudio's memory,
# so you'll need to use mawk on the command line to pick out the lines you want first
# An SBATCH jobscript and bash script have been written to do this on the computing cluster:
# extract_chr-js.sh and extract_chr.sh
# EXAMPLE:
# cat per_read_modified_base_calls.txt | mawk '$2 ~ /^III$/ {print $0}' > chrIII.txt
#######################################################################################

library(data.table)
library(tidyverse)
library(wesanderson)

mega_directory <- "/Volumes/brothersseq/210403_Hello/"
chr <- "XI" #which chromosome?
barcode <- "01"

probs <- fread(sprintf("%schr%s_%s.txt", mega_directory, chr, barcode),
                    select = c(1, 2, 3, 4, 5),
                    col.names = c("read_id", "chrm", "strand", "pos", "mod_log_prob"))

#1. create a binary column; if 10^mod_log_prob is >0.8, set as methylated (m6A). if < 0.8, set as unmethylated (A)
#2. add start and end positions for each read_id
#3. find the average methylation of each read (for ordering on plot)
#4. order by strand, avg methyl and read_id
probs_filtered <- probs[, methylation := ifelse(10^mod_log_prob > 0.8, TRUE, FALSE)][
  , list(read_id, pos, methylation, strand)][
    , start_pos := min(pos), by = read_id][
      , end_pos := max(pos), by = read_id][
        , avg_methyl := mean(methylation == TRUE), by = read_id][
          order(-avg_methyl, read_id)]

# INCLUDE NA's (methylation = NA if prob is between 0.2 and 0.8)
# probs <- probs[, methylation := ifelse(10^mod_log_prob > 0.8, TRUE,
#                                        ifelse(10^mod_log_prob < 0.2, FALSE, NA))][
#   , list(read_id, pos, methylation, strand)][
#     , start_pos := min(pos), by = read_id][
#       , end_pos := max(pos), by = read_id][
#         , avg_methyl := mean(methylation == TRUE, na.rm = TRUE), by = read_id][
#           order(strand, -avg_methyl, read_id)]

#extract each unique read_id to set the order (in same order as in the data table to start with) for single read plots
read_names <- unique(probs_filtered$read_id)
probs_filtered$read_id = factor(probs_filtered$read_id, levels = read_names)

##################################
####PLOT AVERAGE PROBABILITIES####
##################################

plot_prob <- function(input, title) {
  ggplot(input, aes(x = pos, y = avg_prob)) +
    theme_classic() +
    geom_point(position = "jitter") +
    labs(
      title = sprintf("%s region", title),
      x = sprintf("chr %s", chr),
      y = "average log_mod_prob") +
    ylim(0, 1)
}

control_avg <- probs[pos %between% c(100e3, 105e3),
                     .(avg_prob = mean(mod_prob)),
                     by = pos
]
plot_prob(control_avg, "control")

HML_avg <- probs[pos %between% c(11e3, 16e3),
                 .(avg_prob = mean(mod_prob)),
                 by = pos
]
plot_prob(HML_avg, "HML")

HMR_avg <- probs[pos %between% c(291e3, 296e3),
                 .(avg_prob = mean(mod_prob)),
                 by = pos
]
plot_prob(HMR_avg, "HMR")

#########################
####PLOT SINGLE READS####
#########################

#plot function
plot_binary <- function(data) {
  ggplot(data, aes(x = pos, y = read_id, color = methylation)) +
    geom_point(shape = 15) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 15, color = "black", family = "Arial"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_color_manual(values = c("gray90", "mediumpurple4"), na.value = "white")
}

#############################
#### unmethylated region ####
#############################

control_region <- c(95e3, 100e3)
control <- probs_filtered[start_pos <= control_region[1]][end_pos >= control_region[2]][
  pos %between% control_region]
plot_binary(control)

#############
#### HML ####
#############

HML_region <- c(10.5e3, 15.5e3)
HML_E <- c(11237, 11268)
HML_I <- c(14600, 14711)
HML_linkers = c(#9407, 9587.5, 9067, 9747, 9923, 10166, 10331, 
  10585, 10748,
  10944, 11118, 11418, 11645, 12021, 12251, 12436, 12649, 12842,
  13017, 13396, 13558, 13829, 14011, 14221, 14883, 15229, 15406
  #15573, 15984, 16244
)

#read in segments with high and low methylation from steady state data
HML_segments <- data.table(readRDS("~/sequencing/sir_spreading/data/HML_bins.rds"))
HML_highmeth <- HML_segments[meth_status == "high"]
HML_lowmeth <- HML_segments[meth_status == "high"]

#plus or minus strands only
HML_plus <- probs_filtered[strand == "+"][start_pos <= HML_region[1]][end_pos >= HML_region[2]][
  pos %between% HML_region]
HML_plus <- droplevels(HML_plus)
HML_minus <- probs_filtered[strand == "-"][start_pos <= HML_region[1]][end_pos >= HML_region[2]][
  pos %between% HML_region]
HML_minus <- droplevels(HML_minus)

#all reads
HML_all <- probs_filtered[start_pos <= HML_region[1]][end_pos >= HML_region[2]][
  pos %between% HML_region]
HML_all <- droplevels(HML_all)

#plot
hmlp <- plot_binary(HML_all)
hmlp
#annotate silencers (black) and high-meth regions (purple)
hmlp + annotate("rect", xmin = c(HML_highmeth$start), xmax = c(HML_highmeth$end),
                ymin = 0.5, ymax = nlevels(HML_all$read_id)+0.5, alpha = 0.3, fill = "mediumpurple4") +
  annotate("rect", xmin = c(HML_E[1], HML_I[1]), xmax = c(HML_E[2], HML_I[2]),
           ymin = 0.5, ymax = nlevels(HML_all$read_id)+0.5, alpha = 0.2, fill = "black")
hmlp + geom_vline(xintercept = c(mean(HML_E), mean(HML_I)))
hmlp + geom_vline(xintercept = HML_linkers)

hmlp + geom_vline(xintercept = c(mean(HML_E), mean(HML_I), 13282, 13809, 12386, 13018))
#############
#### HMR ####
#############

HMR_region <- c(291e3, 296e3)
HMR_E = c(292674, 292769)
HMR_I = c(294805, 294864)
# HMR_linkers = c(291100, 291312, 291644, 291863, 292129, 292322,
#                 292498, 292921, 293078, 293227, 293440, 293633, 293841, 294155,
#                 294515, 294699, 295239, 295555, 295743, 295906)

#read in segments with high and low methylation from steady state data
HMR_segments <- data.table(readRDS("~/sequencing/sir_spreading/data/HMR_bins.rds"))
HMR_highmeth <- HMR_segments[meth_status == "high"]
HMR_lowmeth <- HMR_segments[meth_status == "low"]

#plus or minus only reads
HMR_plus <- probs_filtered[strand == "+"][start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_plus <- droplevels(HMR_plus)
HMR_minus <- probs_filtered[strand == "-"][start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_minus <- droplevels(HMR_minus)

#all reads
HMR_all <- probs_filtered[start_pos <= HMR_region[1]][end_pos >= HMR_region[2]][
  pos %between% HMR_region]
HMR_all <- droplevels(HMR_all)

#plot
hmrp <- plot_binary(HMR_all)
hmrp
#annotate silencers (black) and high-meth regions (purple)
hmrp + annotate("rect", xmin = c(HMR_highmeth$start), xmax = c(HMR_highmeth$end),
                ymin = 0.5, ymax = nlevels(HMR_all$read_id)+0.5, alpha = 0.2, fill = "mediumpurple4") +
  annotate("rect", xmin = c(HMR_E[1], HMR_I[1]), xmax = c(HMR_E[2], HMR_I[2]),
           ymin = 0.5, ymax = nlevels(HMR_all$read_id)+0.5, alpha = 0.3, fill = "black")

hmrp + geom_vline(xintercept = c(mean(HMR_E), mean(HMR_I), 293835, 294321, 293179, 293538))

###################
#### telomeres ####
###################
#plot tel14L
tel14L_region <- c(0, 10e3)
tel14L <- probs_filtered[start_pos <= tel14L_region[1]+500][end_pos >= tel14L_region[2]][
  pos %between% tel14L_region]
plot_binary(tel14L)

#plot tel11R
tel11R_region <- c(662e3, 666500)
tel11R <- probs_filtered[start_pos <= tel11R_region[1]][end_pos >= tel11R_region[2]][
  pos %between% tel11R_region]
plot_binary(tel11R)

###############################
####HIERARCHICAL CLUSTERING####
###############################

#convert read_ids and positions to factors
HMR_plus$read_id <- factor(HMR_plus$read_id)
HMR_plus$pos <- factor(HMR_plus$pos)
HMRplus_readnames <- unique(HMR_plus$read_id)
HMRplus_positions <- unique(HMR_plus$pos)

#create a new dataframe with positions as column names and read_ids as row names
new_HMRplus <- data.frame()
for(i in HMRplus_positions){
  new_HMRplus[,i] <- NA
}
for(i in HMRplus_readnames){
  new_HMRplus[i,] <- logical()
}

#move the methylation data from the original dataframe to the new dataframe
for(x in HMRplus_readnames){
  for(y in HMRplus_positions){
    if(any(HMR_plus$read_id == x & HMR_plus$pos == y) == FALSE){
      next
    }
    new_HMRplus[as.character(x),as.character(y)] <- HMR_plus[which(HMR_plus$read_id == x & HMR_plus$pos == y),]$methylation
  }
}

#calculate distances and plot the resulting dendrogram
HMRplus_distances <- dist(new_HMRplus, method = "binary")
plot(hclust(HMRplus_distances, method = "average"))

# TO DO: figure out how to apply the clustering/distances to levels to plot the single reads by similarity/difference
