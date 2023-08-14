#!/bin/R

#########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2021-04-06
#########################

# With an input bedMethyl file:
# (https://www.encodeproject.org/data-standards/wgbs/),
# this script will plot the percentage of methylation
# at each position as a scatterplot

library(data.table)
library(tidyverse)
library(stringr)

columns <- c("chrom", "start", "end", "name", "score", 
             "strand", "startCodon", "stopCodon", "color", 
             "coverage", "percentage")

dt_1 <- fread("/Volumes/brothers_seq/210403_Hello/megalodon_output_01/modified_bases.aggregate01.6mA.bed")
colnames(dt_1) <- columns

dt_2 <- fread("/Volumes/brothers_seq/210403_Hello/megalodon_output_02/modified_bases.aggregate02.6mA.bed")
colnames(dt_2) <- columns

dt_3 <- fread("/Volumes/brothers_seq/210403_Hello/megalodon_output_03/modified_bases.aggregate03.6mA.bed")
colnames(dt_3) <- columns

# dt_4 <- fread("/Volumes/brothers_seq/210310_Russula/megalodon_output_04/modified_bases.aggregate04.6mA.bed")
# colnames(dt_4) <- columns
# 
# dt_5 <- fread("/Volumes/brothers_seq/210310_Russula/megalodon_output_05/modified_bases.aggregate05.6mA.bed")
# colnames(dt_5) <- columns
# 
# dt_6 <- fread("/Volumes/brothers_seq/210310_Russula/megalodon_output_06/modified_bases.aggregate06.6mA.bed")
# colnames(dt_6) <- columns

#get a new data.table only containing chrom, start, coverage, and percentage
#filter out mitochondrial DNA (MT) and coverage < 10

select_cols <- c("chrom", "start", "coverage", "percentage")
relevant_1 <- dt_1[chrom != "MT" & coverage > 10, ..select_cols]
relevant_2 <- dt_2[chrom != "MT" & coverage > 10, ..select_cols]
relevant_3 <- dt_3[chrom != "MT" & coverage > 10, ..select_cols]
# relevant_4 <- dt_4[chrom != "MT" & coverage > 10, ..select_cols]
# relevant_5 <- dt_5[chrom != "MT" & coverage > 10, ..select_cols]
# relevant_6 <- dt_6[chrom != "MT" & coverage > 10, ..select_cols]

#plot percentage of methylation in a particular region as a scatter plot
#opacity of dot corresponds to amount of coverage
plot_title = "sir3-8-EcoGII"

plot_methylation_dot <- function(data, chr) {
  ggplot(data, aes(x = start, y = percentage)) +
    geom_point(color = "mediumpurple4", alpha = 1/4) +
    scale_x_continuous(limits = c(min(data$start), max(data$start)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size=15, color = "black", family = "Arial")) +
    labs(title = plot_title,
         x = sprintf("position on chr %s", chr),
         y = "% methylated reads")
}

plot_methylation_bars <- function(data, data2, chr) {
  #break data into 200bp bins, find the mean and sd of % methylated reads in each bin
  numbreaks <- as.integer((max(data$start) - min(data$start)) / 100)
  data$window <- cut(data$start, breaks = numbreaks)
  data_mean <- aggregate(data$percentage, by = list(data$window), function(x){ mean(x)})
  data_sd <- aggregate(data$percentage, by = list(data$window), function(x){ sd(x)})
  
  numbreaks2 <- as.integer((max(data2$start) - min(data2$start)) / 100)
  data2$window <- cut(data2$start, breaks = numbreaks2)
  data2_mean <- aggregate(data2$percentage, by = list(data2$window), function(x){ mean(x)})
  data2_sd <- aggregate(data2$percentage, by = list(data2$window), function(x){ sd(x)})
  
  #make a data.frame of the stats found above and clean up for plotting
  data_stats <- data.frame(interval = data_mean[,1], mean = data_mean[,2], upper = data_mean[,2]+data_sd[, 2], lower = data_mean[,2]-data_sd[, 2])
  data_stats$interval <- gsub("[(]", "", data_stats$interval)
  data_stats$interval <- gsub("[,].*", "", data_stats$interval)
  data_stats$interval <- as.numeric(data_stats$interval)
  
  data2_stats <- data.frame(interval = data2_mean[,1], mean = data2_mean[,2], upper = data2_mean[,2]+data2_sd[, 2], lower = data2_mean[,2]-data2_sd[, 2])
  data2_stats$interval <- gsub("[(]", "", data2_stats$interval)
  data2_stats$interval <- gsub("[,].*", "", data2_stats$interval)
  data2_stats$interval <- as.numeric(data2_stats$interval)
  
  #get the points that have 100% reads methylated
  data_full <- data[data$percentage == 100, ]
  data2_full <- data2[data2$percentage == 100, ]
  
  #plot the bins and the points that are at 100%
  ggplot(mapping = aes(x = interval, y = mean)) +
    ylim(c(-20,100)) +
    geom_point(data = data_stats, color = "mediumpurple4") +
    geom_point(data = data2_stats, color = "cyan3") +
    geom_point(data = data_full, mapping = aes(x = start, y = percentage), color = "mediumpurple4", alpha = 0.3, inherit.aes = FALSE) +
    geom_point(data = data2_full, mapping = aes(x = start, y = percentage), color = "cyan3", alpha = 0.3, inherit.aes = FALSE) +
    geom_errorbar(data = data_stats, mapping = aes(ymin = lower, ymax = upper), color = "mediumpurple4", inherit.aes = TRUE) + 
    geom_errorbar(data = data2_stats, mapping = aes(ymin = lower, ymax = upper), color = "cyan3", inherit.aes = TRUE) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size = 15, color = "black", family = "Arial")) +
    labs(title = plot_title,
         x = sprintf("position on chr %s", chr),
         y = "average % methylated reads")
}

plot_methylation_bars_alone <- function(data, chr) {
  #break data into 200bp bins, find the mean and sd of % methylated reads in each bin
  numbreaks <- as.integer((max(data$start) - min(data$start)) / 100)
  data$window <- cut(data$start, breaks = numbreaks)
  data_mean <- aggregate(data$percentage, by = list(data$window), function(x){ mean(x)})
  data_sd <- aggregate(data$percentage, by = list(data$window), function(x){ sd(x)})
  
  #make a data.frame of the stats found above and clean up for plotting
  data_stats <- data.frame(interval = data_mean[,1], mean = data_mean[,2], upper = data_mean[,2]+data_sd[, 2], lower = data_mean[,2]-data_sd[, 2])
  data_stats$interval <- gsub("[(]", "", data_stats$interval)
  data_stats$interval <- gsub("[,].*", "", data_stats$interval)
  data_stats$interval <- as.numeric(data_stats$interval)
  
  #get the points that have 100% reads methylated
  data_full <- data[data$percentage == 100, ]
  
  #plot the bins and the points that are at 100%
  ggplot(mapping = aes(x = interval, y = mean)) +
    ylim(c(-15,100)) +
    geom_point(data = data_stats, color = "mediumpurple4") +
    geom_point(data = data_full, mapping = aes(x = start, y = percentage), color = "mediumpurple4", alpha = 0.3, inherit.aes = FALSE) +
    geom_errorbar(data = data_stats, mapping = aes(ymin = lower, ymax = upper), color = "mediumpurple4", inherit.aes = TRUE) + 
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size = 15, color = "black", family = "Arial")) +
    labs(title = plot_title,
         x = sprintf("position on chr %s", chr),
         y = "average % methylated reads")
}

#plot a negative control region
control_1 <- relevant_1[chrom == "III" & start > 80e3 & start < 105e3]
control_2 <- relevant_2[chrom == "III" & start > 80e3 & start < 105e3]
control_3 <- relevant_3[chrom == "III" & start > 80e3 & start < 105e3]
plot_methylation_bars_alone(control_3, "III")
plot_methylation_bars(control_2, control_1, "III")
plot_methylation_dot(control, "III")

#plot HML
HML_1 <- relevant_1[chrom == "III" & start > 0 & start < 25e3]
HML_2 <- relevant_2[chrom == "III" & start > 0 & start < 25e3]
HML_3 <- relevant_3[chrom == "III" & start > 0 & start < 25e3]
# HML_4 <- relevant_4[chrom == "III" & start > 0 & start < 25e3]
# HML_5 <- relevant_5[chrom == "III" & start > 0 & start < 25e3]
# HML_6 <- relevant_6[chrom == "III" & start > 0 & start < 25e3]
plot_methylation_bars_alone(HML_3, "III") + geom_vline(xintercept = c(11146, 14849))
plot_methylation_bars(HML_3, HML_2, "III")
plot_methylation_dot(HML_3, "III")

#plot HMR
HMR_1 <- relevant_1[chrom == "III" & start > 280e3 & start < 305e3]
HMR_2 <- relevant_2[chrom == "III" & start > 280e3 & start < 305e3]
HMR_3 <- relevant_3[chrom == "III" & start > 280e3 & start < 305e3]
# HMR_4 <- relevant_4[chrom == "III" & start > 280e3 & start < 305e3]
# HMR_5 <- relevant_5[chrom == "III" & start > 280e3 & start < 305e3]
# HMR_6 <- relevant_6[chrom == "III" & start > 280e3 & start < 305e3]
plot_methylation_bars_alone(HMR_3, "III") + geom_vline(xintercept = c(292388, 295034))
plot_methylation_bars(HMR_3, HMR_2, "III") + geom_vline(xintercept = c(292388, 295034))
plot_methylation_dot(HMR_37C, "III")

#plot telomere
tel6R_1 <- relevant_1[chrom == "VI" & start > 250e3 & start < 270e3]
tel6R_2 <- relevant_2[chrom == "VI" & start > 250e3 & start < 270e3]
tel6R_3 <- relevant_3[chrom == "VI" & start > 250e3 & start < 270e3]
plot_methylation_bars_alone(tel6R_2, "VI")
plot_methylation_bars(tel6R_2, tel6R_3, "VI")
plot_methylation_dot(tel6R, "VI")

#plot rDNA
rdna_1 <- relevant_1[chrom == "XII" & start > 445e3 & start < 475e3]
rdna_2 <- relevant_2[chrom == "XII" & start > 445e3 & start < 475e3]
plot_methylation_bars(rdna_1, rdna_2, "XII")

#plot chromosome
chrXV <- relevant[chrom == "XV"]
plot_methylation_dot(chrXV, "XV")
plot_methylation_bars(chrXV, "XV")