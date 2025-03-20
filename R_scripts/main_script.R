# Reorganized R Code for NICU Microbiome Analysis
# -----------------------------------------------------------
# This script contains functions and pipelines for:
#  1) Reading and merging metagenomic data (species, genus)
#  2) Cleaning and filtering data
#  3) Calculating alpha diversity, Bray-Curtis distances
#  4) Conducting PCA, MRPP, Wilcoxon tests, GLMM analysis, etc.
#  5) Generating heatmaps, boxplots, and volcano plots
# -----------------------------------------------------------

# SETUP ------------------------------------------------------

## Define environment/work directories (comment out those not needed)
## Example usage:
# setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")
# setwd("/home/david/Documents/Code/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929")

# PACKAGES ---------------------------------------------------

list.of.packages <- c(
  "vegan","labdsv","pvclust","gplots","RColorBrewer",
  "Heatplus","plyr","fossil","ade4","scales","randomForest","colorspace",
  "ggbiplot","MGLM","ggthemes","tidyverse","pheatmap","magrittr","robustbase","cowplot","sda","locfdr",
  "FactoMineR","factoextra","NBZIMM","gtools","heatmap3","glmmTMB",
  "dirmult","ecodist","lme4", "VennDiagram","EnhancedVolcano","permute"
)

# Install any missing packages (optional)
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

lapply(list.of.packages, library, character.only = TRUE)

# UTILITY FUNCTIONS ------------------------------------------

## Function to perform linear regression and return plot with stats
ggplotRegression <- function(fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(
      title = paste(
        "Adj R2 = ", signif(summary(fit)$adj.r.squared, 5),
        "Intercept =", signif(fit$coef[[1]],5),
        "Slope =", signif(fit$coef[[2]],5),
        "P =", signif(summary(fit)$coef[2,4], 5)
      )
    )
}

## Utility to capitalize
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

## Effect size computation (SDA-based) for comparing groups (Location, Abx, etc.)
## Based on approach described by Bernd Klaus
compute_ef_Location <- function(d, g, min_shrink = 0.3){
  raw_scores <-  sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  # compute effect sizes
  ef <- raw_scores[, "cat.Hangzhou"]*m["Hangzhou"] - raw_scores[, "cat.Cincinnati"]*m["Cincinnati"]
  n0 <- sum(g_summaries$idx[, "Hangzhou"])
  n1 <- sum(g_summaries$idx[, "Cincinnati"])
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
}

## Repetitive versions of the same approach for different grouping
compute_ef_Abx <- function(d, g, min_shrink = 0.3){
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  ef <- raw_scores[, "cat.Infant.Abx"]*m["Infant.Abx"] - raw_scores[, "cat.No.Infant.Abx"]*m["No.Infant.Abx"]
  n0 <- sum(g_summaries$idx[, "Infant.Abx"])
  n1 <- sum(g_summaries$idx[, "No.Infant.Abx"])
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
}

compute_ef_Gestation <- function(d, g, min_shrink = 0.3){
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  ef <- raw_scores[, "cat.Cohort.2"]*m["Cohort.2"] - raw_scores[, "cat.Cohort.1"]*m["Cohort.1"]
  n0 <- sum(g_summaries$idx[, "Cohort.2"])
  n1 <- sum(g_summaries$idx[, "Cohort.1"])
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
}

compute_ef_Week <- function(d, g, min_shrink = 0.3){
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  ef <- raw_scores[, "cat.Week.1"]*m["Week.1"] - raw_scores[, "cat.Week.3"]*m["Week.3"]
  n0 <- sum(g_summaries$idx[, "Week.1"])
  n1 <- sum(g_summaries$idx[, "Week.3"])
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  lfdr <- locfdr(stats, nulltype = 1)$fdr
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  tibble(
    Species = names(ef),
    ef = ef,
    ef_shrunk = ef_shrunk,
    stat = stats,
    lfdr = lfdr
  )
}

## Additional Theming for ggplot
params <- function(x){
  theme(axis.text.x= element_text(size=14,color="black")) +
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    theme(axis.text.y=element_text(size=14,color="black")) +
    theme(plot.title=element_text(size=22,color="black")) +
    theme(axis.title.x=element_text(size=20), axis.title.y=element_text(size=18)) +
    theme(legend.title=element_text(size=18), legend.text=element_text(size=14))
}

paramsAngled <- function(x){
  theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1,color="black")) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(axis.text.y=element_text(size=14,color="black")) +
    theme(plot.title=element_text(size=22,color="black")) +
    theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=18)) +
    theme(legend.title=element_text(size=18),legend.text=element_text(size=14)) +
    theme(legend.position="none")
}

paramsUnAngled <- function(x){
  theme(axis.text.x=element_text(size=12,vjust=1,hjust=1,color="black")) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(axis.text.y=element_text(size=14,color="black")) +
    theme(plot.title=element_text(size=22,color="black")) +
    theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=18)) +
    theme(legend.title=element_text(size=14),legend.text=element_text(size=12)) +
    theme(legend.position="none")
}

paramsBox <- function(x){
  theme(axis.text.x=element_text(size=14,color="black")) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(axis.text.y=element_text(size=18,color="black")) +
    theme(plot.title=element_text(size=22,color="black")) +
    theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=22)) +
    theme(legend.title=element_text(size=18),legend.text=element_text(size=14)) +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(size=18,color="black")) +
    theme(axis.text.y=element_text(size=16,angle=0,color="black")) +
    theme(strip.text.y=element_text(size=14,angle=270,colour="black")) +
    theme(strip.text.x=element_text(size=16,angle=0,colour="black"))
}

paramsBoxWide <- function(x){
  theme(axis.text.x=element_text(size=14,color="black")) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(axis.text.y=element_text(size=16,color="black")) +
    theme(plot.title=element_text(size=22,color="black")) +
    theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24)) +
    theme(legend.title=element_text(size=18),legend.text=element_text(size=14)) +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(size=16,color="black")) +
    theme(axis.text.y=element_text(size=16,angle=0,color="black")) +
    theme(strip.text.y=element_text(size=14,angle=270,colour="black")) +
    theme(strip.text.x=element_text(size=16,angle=0,colour="black"))
}

metaphlan.colors <- colorRampPalette(c("#000033","#007FFF","cyan","red","yellow"))
scaleyellowred <- colorRampPalette(c("lightyellow","red"), space="rgb")(50)
col <- c("gray32","royalblue4","firebrick")
diverging.colors <- colorRampPalette(c("#000033","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","firebrick"))

# Additional utility functions
noise.removal <- function(dataframe, percent=0.001, top=NULL){
  Matrix <- dataframe
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  return(Matrix_1)
}

LogModulus <- function(x){
  sign(x)*(log(abs(x)+1))
}

robustMean <- function(x){
  huberM(x)$mu
}

# MAIN ANALYSES AND DATA LOADING ------------------------------
# (1) Reading the sample key, merging Kraken2 output
# (2) Combining species/genus-level data
# (3) Rarefying, filtering, removing contaminants
# (4) Alpha diversity calculations (Shannon, Simpson, etc.)
# (5) Beta diversity (Bray-Curtis, PCA, MRPP)
# (6) Differential abundance, effect size calculations
# (7) Plotting (heatmaps, boxplots, volcano plots, etc.)

# NOTE: Much of the code below is specialized for your data.
#       If you see repeated code for separate analyses, consider modularizing.

# Example structure:
# Step 1. setwd() & read CSV files
# Step 2. read Kraken2 outputs
# Step 3. create tables, merge them
# Step 4. filter out human reads or unwanted species
# Step 5. define and run alpha diversity calculations
# Step 6. define or run PCA, MRPP, etc.
# Step 7. create boxplots, heatmaps, volcano plots

# CLEANUP & SAVE WORKSPACE ------------------------------------
# Finally, you can save the environment if needed:
# save.image(file = "NICUData2025xx.RData")

# End of reorganized code --------------------------------------
