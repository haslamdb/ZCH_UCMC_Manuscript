#!/usr/bin/env Rscript
# NICU Microbiome Analysis Utilities
# -----------------------------------------------------------
# This script contains utility functions for microbiome analysis:
#  1) Plotting functions and parameters
#  2) Effect size calculation utilities 
#  3) Data filtering and transformation functions
#  4) Statistical helper functions
# -----------------------------------------------------------

# Load required packages
load_packages <- function() {
  list.of.packages <- c(
    "vegan", "labdsv", "pvclust", "gplots", "RColorBrewer",
    "Heatplus", "plyr", "fossil", "ade4", "scales", 
    "randomForest", "colorspace", "ggbiplot", "MGLM", 
    "ggthemes", "tidyverse", "pheatmap", "magrittr", 
    "robustbase", "cowplot", "sda", "locfdr",
    "FactoMineR", "factoextra", "NBZIMM", "gtools", 
    "heatmap3", "glmmTMB", "dirmult", "ecodist", 
    "lme4", "VennDiagram", "EnhancedVolcano", "permute"
  )
  
  # Check and install missing packages (uncomment if needed)
  # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  # if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)
  
  # Load packages
  suppressMessages(invisible(lapply(list.of.packages, library, character.only = TRUE)))
  
  cat("Packages loaded successfully!\n")
}

# Color Palettes and Visualization Themes
# ---------------------------------------

# Tableau-inspired color palette
tableau.colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", 
                   "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7")

# Function to use tableau colors in ggplot
scale_colour_tableau <- function() {
  scale_colour_manual(values = tableau.colors)
}

# Standard color schemes
metaphlan.colors <- colorRampPalette(c("#000033", "#007FFF", "cyan", "red", "yellow"))
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
col <- c("gray32", "royalblue4", "firebrick")
diverging.colors <- colorRampPalette(c("#000033", "#4575b4", "#74add1", "#abd9e9", 
                                     "#e0f3f8", "#fee090", "#fdae61", "#f46d43", 
                                     "#d73027", "firebrick"))

# ggplot Themes for Different Plot Types
# ---------------------------------------

# Standard parameters
params <- function() {
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 14)
  )
}

# Parameters with angled x-axis labels
paramsAngled <- function() {
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.position = "none"
  )
}

# Parameters without angled labels
paramsUnAngled <- function() {
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none"
  )
}

# Parameters for box plots
paramsBox <- function() {
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.position = "none",
    strip.text.y = element_text(size = 14, angle = 270, colour = "black"),
    strip.text.x = element_text(size = 16, angle = 0, colour = "black")
  )
}

# Parameters for wide box plots
paramsBoxWide <- function() {
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, angle = 0, color = "black"),
    plot.title = element_text(size = 22, color = "black"),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.position = "none",
    strip.text.y = element_text(size = 14, angle = 270, colour = "black"),
    strip.text.x = element_text(size = 16, angle = 0, colour = "black")
  )
}

# Data Processing Functions
# ---------------------------------------

# Function to remove low-abundance features
noise.removal <- function(dataframe, percent = 0.001, top = NULL) {
  Matrix <- dataframe
  bigones <- rowSums(Matrix) * 100 / (sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones, ]
  return(Matrix_1)
}

# Log modulus transformation for data with positive and negative values
LogModulus <- function(x) {
  sign(x) * (log(abs(x) + 1))
}

# Robust mean calculation
robustMean <- function(x) {
  huberM(x)$mu
}

# Capitalize first letter of each word
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

# Fold change calculation
foldchange <- function(x, y) {
  # Returns fold change between x and y
  # If x > y, returns (x/y)
  # If y > x, returns -(y/x)
  # Handles zero values gracefully
  
  result <- numeric(length(x))
  
  for (i in 1:length(x)) {
    if (x[i] == 0 && y[i] == 0) {
      result[i] <- 0
    } else if (x[i] == 0) {
      result[i] <- -Inf
    } else if (y[i] == 0) {
      result[i] <- Inf
    } else if (x[i] > y[i]) {
      result[i] <- x[i] / y[i]
    } else {
      result[i] <- -y[i] / x[i]
    }
  }
  
  return(result)
}

# Convert fold change to log ratio
foldchange2logratio <- function(fc, base = 2) {
  # Converts fold change to log ratio
  # Handles positive and negative fold changes
  
  result <- numeric(length(fc))
  
  for (i in 1:length(fc)) {
    if (fc[i] > 0) {
      result[i] <- log(fc[i], base)
    } else if (fc[i] < 0) {
      result[i] <- -log(abs(fc[i]), base)
    } else {
      result[i] <- 0
    }
  }
  
  return(result)
}

# Visualization Functions
# ---------------------------------------

# Function to create regression plot with statistics
ggplotRegression <- function(fit) {
  require(ggplot2)
  
  # Get model summary data
  coef <- coefficients(fit)
  formula <- sprintf("y = %.3fx + %.3f", coef[2], coef[1])
  r2 <- sprintf("R^2 = %.3f", summary(fit)$r.squared)
  p <- sprintf("p = %.3f", summary(fit)$coefficients[2,4])
  
  # Create the plot
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste0(formula, "\n", r2, "\n", p))
}

# Statistical Analysis Functions
# ---------------------------------------

# Effect size computation for Location comparisons
compute_ef_Location <- function(d, g, min_shrink = 0.3) {
  # Computes effect sizes for Location comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  require(tibble)
  
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  # Compute effect sizes
  ef <- raw_scores[, "cat.Hangzhou"] * m["Hangzhou"] - raw_scores[, "cat.Cincinnati"] * m["Cincinnati"]
  n0 <- sum(g_summaries$idx[, "Hangzhou"])
  n1 <- sum(g_summaries$idx[, "Cincinnati"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
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

# Effect size computation for Antibiotic comparisons
compute_ef_Abx <- function(d, g, min_shrink = 0.3) {
  # Computes effect sizes for Antibiotic comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  require(tibble)
  
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  # Compute effect sizes
  ef <- raw_scores[, "cat.Infant.Abx"] * m["Infant.Abx"] - raw_scores[, "cat.No.Infant.Abx"] * m["No.Infant.Abx"]
  n0 <- sum(g_summaries$idx[, "Infant.Abx"])
  n1 <- sum(g_summaries$idx[, "No.Infant.Abx"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
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

# Effect size computation for Gestational age comparisons
compute_ef_Gestation <- function(d, g, min_shrink = 0.3) {
  # Computes effect sizes for Gestational age comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  require(tibble)
  
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  # Compute effect sizes
  ef <- raw_scores[, "cat.Cohort.2"] * m["Cohort.2"] - raw_scores[, "cat.Cohort.1"] * m["Cohort.1"]
  n0 <- sum(g_summaries$idx[, "Cohort.2"])
  n1 <- sum(g_summaries$idx[, "Cohort.1"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
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

# Effect size computation for Week comparisons
compute_ef_Week <- function(d, g, min_shrink = 0.3) {
  # Computes effect sizes for Week comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  require(tibble)
  
  raw_scores <- sda.ranking(d, g, diagonal = FALSE, verbose = TRUE, fdr = FALSE)
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs) / freqs / length(g))
  
  # Compute effect sizes
  ef <- raw_scores[, "cat.Week.1"] * m["Week.1"] - raw_scores[, "cat.Week.3"] * m["Week.3"]
  n0 <- sum(g_summaries$idx[, "Week.1"])
  n1 <- sum(g_summaries$idx[, "Week.3"])
  m_stats <- 1 / sqrt(1 / n0 + 1 / n1)
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

# Export all utility functions
utils_export <- function() {
  # This function prints the list of exported functions
  cat("NICU Microbiome Analysis Utilities\n")
  cat("----------------------------------\n")
  cat("Color palettes: tableau.colors, metaphlan.colors, diverging.colors\n")
  cat("Theme functions: params(), paramsAngled(), paramsUnAngled(), paramsBox(), paramsBoxWide()\n")
  cat("Data processing: noise.removal(), LogModulus(), foldchange(), foldchange2logratio()\n")
  cat("Visualization: ggplotRegression()\n")
  cat("Statistical: compute_ef_Location(), compute_ef_Abx(), compute_ef_Gestation(), compute_ef_Week()\n")
}

# Demo function
utils_demo <- function() {
  # Create sample data
  set.seed(123)
  x <- runif(100)
  y <- 2*x + rnorm(100, 0, 0.5)
  df <- data.frame(x = x, y = y)
  
  # Demo regression plot
  fit <- lm(y ~ x, data = df)
  plot <- ggplotRegression(fit)
  print(plot)
  
  # Demo color palettes
  par(mfrow = c(3, 1))
  image(matrix(1:100, ncol=1), col = metaphlan.colors(100), main = "Metaphlan Colors")
  image(matrix(1:100, ncol=1), col = diverging.colors(100), main = "Diverging Colors")
  barplot(rep(1, length(tableau.colors)), col = tableau.colors, main = "Tableau Colors")
  
  # Return to original plot settings
  par(mfrow = c(1, 1))
}
