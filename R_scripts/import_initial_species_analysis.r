#!/usr/bin/env Rscript
# NICU Microbiome Data Import and Analysis
# 
# This script performs the following operations:
# 1. Imports and processes sample metadata
# 2. Imports Kraken2 taxonomic alignments
# 3. Filters and normalizes species count data
# 4. Analyzes diversity metrics and associations
# 5. Performs statistical tests between groups
# 6. Creates visualizations of the results

# Load required libraries
library(vegan)      # For diversity analyses and ordination
library(ggplot2)    # For plotting
library(reshape2)   # For reshaping data
library(FactoMineR) # For PCA
library(factoextra) # For PCA visualization
library(dplyr)      # For data manipulation
library(sda)        # For effect size calculation

# Set your working directory in script or through setwd() call
# Project paths (user can modify this section)
PROJECT_DIR <- "."  # Set this to your project directory
DATA_DIR <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
KRAKEN_DIR <- file.path(PROJECT_DIR, "KrakenAlignments/Kraken2")

# Create directories if they don't exist
dir.create(DATA_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)

# Helper functions
# Function to calculate fold change
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

# Function to compute effect size for Location comparisons
compute_ef_Location <- function(d, g, min_shrink = 0.1) {
  # Computes effect sizes for Location comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  
  g1 <- which(g == "Hangzhou")
  g2 <- which(g == "Cincinnati")
  
  d1 <- d[g1, ]
  d2 <- d[g2, ]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  xbar1 <- colMeans(d1)
  xbar2 <- colMeans(d2)
  
  npool <- n1 + n2 - 2
  s1 <- apply(d1, 2, var)
  s2 <- apply(d2, 2, var)
  
  # Pool and get std error
  spool <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / npool)
  
  # Calculate effect size
  ef <- (xbar1 - xbar2) / spool
  
  # Shrink (regularize) the effect sizes
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - apply(d, 2, var))
  
  data.frame(
    Species = colnames(d),
    ef = ef,
    ef_shrunk = ef_shrunk
  )
}

# Function to compute effect size for Antibiotic comparisons
compute_ef_Abx <- function(d, g, min_shrink = 0.1) {
  # Computes effect sizes for Antibiotic comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  
  g1 <- which(g == "Infant.Abx")
  g2 <- which(g == "No.Infant.Abx")
  
  d1 <- d[g1, ]
  d2 <- d[g2, ]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  xbar1 <- colMeans(d1)
  xbar2 <- colMeans(d2)
  
  npool <- n1 + n2 - 2
  s1 <- apply(d1, 2, var)
  s2 <- apply(d2, 2, var)
  
  # Pool and get std error
  spool <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / npool)
  
  # Calculate effect size
  ef <- (xbar1 - xbar2) / spool
  
  # Shrink (regularize) the effect sizes
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - apply(d, 2, var))
  
  data.frame(
    Species = colnames(d),
    ef = ef,
    ef_shrunk = ef_shrunk
  )
}

# Function to compute effect size for Week comparisons
compute_ef_Week <- function(d, g, min_shrink = 0.1) {
  # Computes effect sizes for Week comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  
  g1 <- which(g == "Week.1")
  g2 <- which(g == "Week.3")
  
  d1 <- d[g1, ]
  d2 <- d[g2, ]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  xbar1 <- colMeans(d1)
  xbar2 <- colMeans(d2)
  
  npool <- n1 + n2 - 2
  s1 <- apply(d1, 2, var)
  s2 <- apply(d2, 2, var)
  
  # Pool and get std error
  spool <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / npool)
  
  # Calculate effect size
  ef <- (xbar1 - xbar2) / spool
  
  # Shrink (regularize) the effect sizes
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - apply(d, 2, var))
  
  data.frame(
    Species = colnames(d),
    ef = ef,
    ef_shrunk = ef_shrunk
  )
}

# Function to compute effect size for Gestation comparisons
compute_ef_Gestation <- function(d, g, min_shrink = 0.1) {
  # Computes effect sizes for Gestation comparisons
  # d: data matrix
  # g: grouping factor
  # min_shrink: minimum shrinkage parameter
  
  require(sda)
  
  g1 <- which(g == "Cohort.1")
  g2 <- which(g == "Cohort.2")
  
  d1 <- d[g1, ]
  d2 <- d[g2, ]
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  xbar1 <- colMeans(d1)
  xbar2 <- colMeans(d2)
  
  npool <- n1 + n2 - 2
  s1 <- apply(d1, 2, var)
  s2 <- apply(d2, 2, var)
  
  # Pool and get std error
  spool <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / npool)
  
  # Calculate effect size
  ef <- (xbar1 - xbar2) / spool
  
  # Shrink (regularize) the effect sizes
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - apply(d, 2, var))
  
  data.frame(
    Species = colnames(d),
    ef = ef,
    ef_shrunk = ef_shrunk
  )
}

# Function for regression plot
ggplotRegression <- function(fit) {
  # Create a regression plot
  # fit: lm model object
  
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

# Plot parameters
paramsBox <- function() {
  # Standard box plot parameters
  theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

paramsBoxWide <- function() {
  # Wide box plot parameters
  theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right")
}

paramsAngled <- function() {
  # Angled text parameters
  theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# Color schemes
tableau.colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7")
scale_colour_tableau <- function() {
  scale_colour_manual(values = tableau.colors)
}

# Load and process sample metadata
print("Loading and processing sample metadata...")

# Import the sample key
sample_key_file <- file.path(DATA_DIR, "AllNICUSampleKeyRevised20250206_for_HangzhouCincinnatiSamples.csv")
NewSampleKey <- read.csv(sample_key_file, header = TRUE, stringsAsFactors = FALSE)
NewSampleKey$PostNatalAntibiotics <- factor(NewSampleKey$PostNatalAntibiotics, levels = c("No.Infant.Abx", "Infant.Abx"))
NewSampleKey$PostNatalAntibioticsNew <- factor(NewSampleKey$PostNatalAntibioticsNew, 
                                              levels = c("No.Infant.Abx", "Low.Infant.Abx", "High.Infant.Abx"))
NewSampleKey$Sample <- NewSampleKey$SampleID


# Combine sample keys
AllSampleKey <- NewSampleKey

# Set factor levels for categorical variables
AllSampleKey$GestationCohort <- factor(AllSampleKey$GestationCohort, 
                                      levels = c("23-27 Weeks", "28-32 Weeks", "33-36 Weeks", "37-42 Weeks"))
AllSampleKey$SampleCollectionWeek <- factor(AllSampleKey$SampleCollectionWeek, levels = c("Week.1", "Week.3"))
AllSampleKey$SampleType <- factor(AllSampleKey$SampleType, levels = c("Nares", "Axilla", "Groin", "Stool"))
AllSampleKey$MaternalAntibiotics <- factor(AllSampleKey$MaternalAntibiotics, levels = c("No.Mat.Abx", "Mat.Abx"))

# Create combined factors for analysis
AllSampleKey$TypeWeek <- paste(AllSampleKey$SampleType, AllSampleKey$SampleCollectionWeek, sep = "-")
AllSampleKey$TypeAbx <- paste(AllSampleKey$SampleType, AllSampleKey$PostNatalAntibiotics, sep = "-")
AllSampleKey$TypeAbx <- factor(AllSampleKey$TypeAbx, 
                              levels = c("Axilla-No.Infant.Abx", "Axilla-Infant.Abx", 
                                         "Groin-No.Infant.Abx", "Groin-Infant.Abx", 
                                         "Stool-No.Infant.Abx", "Stool-Infant.Abx", 
                                         "Nares-No.Infant.Abx", "Nares-Infant.Abx"))

# Filter for unique samples and relevant sample types
AllSampleKey <- subset(AllSampleKey, !duplicated(AllSampleKey$SampleID) & 
                      AllSampleKey$SampleType %in% c("Stool", "Groin", "Axilla"))

# Import human reactive species data
human_reactive_file <- file.path(DATA_DIR, "HumanReactiveKraken2.csv")
HumanReactiveSpecies <- read.csv(human_reactive_file)

# Clean species names for matching
HumanReactiveSpecies$SubSpecies <- gsub("_", ".", HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies <- gsub("-", ".", HumanReactiveSpecies$SubSpecies)
HumanReactiveSpecies$SubSpecies <- gsub("/", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("X,", "", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("[", "", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("]", "", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("..", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies <- gsub("..", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)
HumanReactiveSpecies$SubSpecies <- gsub(" ", ".", HumanReactiveSpecies$SubSpecies, fixed = TRUE)

# Load Kraken2 alignments
print("Loading Kraken2 alignments...")
setwd(KRAKEN_DIR)

# Get list of Kraken2 species abundance files
AllK2Files <- list.files()
SpeciesFileList3 <- grep("_species_abundance.txt", AllK2Files)
K2SpeciesFiles <- AllK2Files[SpeciesFileList3]
K2FileList <- gsub("_species_abundance.txt", "", K2SpeciesFiles)

# Filter for NICU samples
NICUFiles <- as.character(AllSampleKey$SampleID)
K2NICUFiles <- K2FileList[K2FileList %in% NICUFiles]
ALLNICUFiles <- K2NICUFiles

# Identify missing NICU files
MissingNICU <- subset(NICUFiles, !NICUFiles %in% c(K2NICUFiles))
# write.csv(MissingNICU, file = file.path(RESULTS_DIR, "MissingNICUFiles.txt"))

# Load all species abundance files
for (f in 1:length(K2NICUFiles)) {
  fnr <- paste0(K2NICUFiles[f], "_species_abundance.txt")
  filename <- K2NICUFiles[f]
  
  assign(filename, read.csv(fnr, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
}

# Process and merge species data
print("Processing species abundance data...")

# Start with first file to create the base table
NewSpeciesTable <- read.csv(paste0(K2NICUFiles[1], "_species_abundance.txt"), sep = "\t", header = TRUE)
names(NewSpeciesTable)[1] <- "Species"

# Clean species names
NewSpeciesTable$Species <- gsub(" ", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("..", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("_", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("-", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("/", ".", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("X,", "", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("[", "", NewSpeciesTable$Species, fixed = TRUE)
NewSpeciesTable$Species <- gsub("]", "", NewSpeciesTable$Species, fixed = TRUE)

# Extract species and counts
NewSpeciesTable <- as.data.frame(NewSpeciesTable[,c(1,6)])
names(NewSpeciesTable) <- c("Species", paste(K2NICUFiles[1]))
NewSpeciesTable <- subset(NewSpeciesTable, !duplicated(NewSpeciesTable$Species))

# Merge data from all files
for (i in 2:length(K2NICUFiles)) {
  x <- get(K2NICUFiles[i])
  x <- as.data.frame(x[,c(1,6)])
  names(x) <- c("Species", paste(K2NICUFiles[i]))
  
  # Clean species names
  x$Species <- gsub(" ", ".", x$Species, fixed = TRUE)
  x$Species <- gsub("..", ".", x$Species, fixed = TRUE)
  x$Species <- gsub("_", ".", x$Species, fixed = TRUE)
  x$Species <- gsub("-", ".", x$Species, fixed = TRUE)
  x$Species <- gsub("/", ".", x$Species, fixed = TRUE)
  x$Species <- gsub("X,", "", x$Species, fixed = TRUE)
  x$Species <- gsub("[", "", x$Species, fixed = TRUE)
  x$Species <- gsub("]", "", x$Species, fixed = TRUE)
  
  x <- subset(x, !duplicated(x$Species))
  
  # Merge with existing table
  NewSpeciesTable <- merge(x, NewSpeciesTable, by = "Species", all = TRUE)
  NewSpeciesTable[is.na(NewSpeciesTable)] <- 0
}

# Convert to matrix format
row.names(NewSpeciesTable) <- NewSpeciesTable$Species
NewSpeciesTable$Species <- NULL
NewSpeciesTable[is.na(NewSpeciesTable)] <- 0

# Clean up temporary objects
rm(list = K2NICUFiles)

# Transpose species table for analysis
AllSpeciesRaw <- t(NewSpeciesTable)
Speciesdf <- as.data.frame(AllSpeciesRaw)
Speciesdf$SampleID <- row.names(AllSpeciesRaw)
SpeciesTable <- as.data.frame(t(NewSpeciesTable))

# Set working directory back to project root
setwd(PROJECT_DIR)

# Filter and normalize species data
print("Filtering and normalizing species data...")

# Create a filtered species table removing human DNA and low-abundance species
NICUSpecies <- as.data.frame(AllSpeciesRaw)

# Remove unwanted species
SaliniCol <- grep("Salinibacter.ruber", names(NICUSpecies))
if (length(SaliniCol) > 0) NICUSpecies <- NICUSpecies[,-SaliniCol]

HumanCol <- grep("Homo.sapiens", names(NICUSpecies))
if (length(HumanCol) > 0) NICUSpecies <- NICUSpecies[,-HumanCol]

# Remove human reactive species
HumanReactiveCols <- which(colnames(NICUSpecies) %in% HumanReactiveSpecies$SubSpecies)
if (length(HumanReactiveCols) > 0) NICUSpecies <- NICUSpecies[,-HumanReactiveCols]

# Filter for species present in at least 5% of samples
TenPercentCutoff <- floor(nrow(NICUSpecies) / 20)

NonZeroCounts <- list()
for (i in 1:ncol(NICUSpecies)) {
  NonZeroCounts[i] <- length(which(NICUSpecies[,i] > 0))
}

TenPercentNotZero <- which(NonZeroCounts >= TenPercentCutoff)
NICUSpeciesReduced <- NICUSpecies[,TenPercentNotZero]

# Remove samples with low read counts
LowSamples <- which(rowSums(NICUSpeciesReduced) <= 100000)
if (length(LowSamples) > 0) NICUSpeciesReduced <- NICUSpeciesReduced[-LowSamples,]

# Noise removal and rarefaction
NICUSpeciesNR <- as.data.frame(t(noise.removal(t(NICUSpeciesReduced), 0.001)))
minCount <- min(rowSums(NICUSpeciesNR))
NICUSpeciesNR <- data.frame(rrarefy(NICUSpeciesReduced, minCount))
NICUSpeciesNR$SampleID <- row.names(NICUSpeciesNR)

# Merge metadata with species data
NICUSpeciesNR <- merge(AllSampleKey, NICUSpeciesNR, by = "SampleID", all.x = TRUE)
NICUSpeciesNR <- subset(NICUSpeciesNR, !duplicated(NICUSpeciesNR$SampleID))

# Create long format species table for plotting
LongSpeciesTable <- melt(NICUSpeciesNR, id.vars = c("SampleID", "PatientID", "Sample", "SampleType", 
                                                  "Location", "GestationalAge", "GestationCohort", 
                                                  "GestationTime", "SampleCollectionWeek", 
                                                  "PostNatalAntibiotics", "MaternalAntibiotics", 
                                                  "Week1Abx", "Week3Abx", "PostNatalAntibioticsNew", 
                                                  "TypeWeek", "TypeAbx", "AnyMilk"))
names(LongSpeciesTable)[c(18,19)] <- c("Species", "Count")
LongSpeciesTable <- subset(LongSpeciesTable, !LongSpeciesTable$Count == 0)
write.csv(LongSpeciesTable, file = file.path(RESULTS_DIR, "LongSpeciesTableNoHuman.csv"))

# Process non-rarefied species data for diversity analysis
NICUSpeciesReduced$SampleID <- row.names(NICUSpeciesReduced) 
NICUSpeciesReduced <- subset(NICUSpeciesReduced, !duplicated(NICUSpeciesReduced$SampleID))

NonRarefiedSpecies <- merge(AllSampleKey, NICUSpeciesReduced, by = "SampleID", all.y = TRUE)
NonRarefiedSpecies <- subset(NonRarefiedSpecies, !duplicated(NonRarefiedSpecies$SampleID))

# Prepare for diversity analysis
NonRarefiedSpecies <- subset(NonRarefiedSpecies, !is.na(NonRarefiedSpecies$SampleType))
row.names(NonRarefiedSpecies) <- NonRarefiedSpecies$SampleID
SampleData <- NonRarefiedSpecies[, 1:17]
NonRarefiedSpecies <- NonRarefiedSpecies[, 18:ncol(NonRarefiedSpecies)]

DiversitySpecies <- NonRarefiedSpecies

# Calculate diversity metrics
print("Calculating diversity metrics...")

# Simpson diversity
simpson <- data.frame(diversity(DiversitySpecies, "simpson"))

# Shannon diversity
shannon <- data.frame(diversity(DiversitySpecies, "shannon"))

# Inverse Simpson
invsimp <- data.frame(diversity(DiversitySpecies, "inv"))

# Fisher's alpha
alpha <- data.frame(fisher.alpha(DiversitySpecies))

# Species richness and evenness
S <- data.frame(specnumber(DiversitySpecies))
H <- data.frame(diversity(DiversitySpecies))
J <- data.frame(H/log(S))

# Combine diversity metrics
Diversity <- cbind(simpson, shannon, invsimp, alpha, S, J)
Diversity$SampleID <- row.names(Diversity)
names(Diversity) <- c("Simpson", "Shannon", "InvSimpson", "Alpha", "SpeciesNo", "Evenness", "SampleID")
Diversity <- Diversity[,c(7,1,2,3,4,5,6)]

# Check correlations between diversity metrics
# pairs(cbind(simpson, shannon, alpha, J, S), pch="+", col="blue")

# Merge diversity metrics with sample data
Diversity <- merge(SampleData, Diversity, by = "SampleID", all.y = TRUE)
Diversity <- subset(Diversity, !duplicated(Diversity$SampleID))

# Save diversity data
write.csv(Diversity, file = file.path(RESULTS_DIR, "DiversityMetrics.csv"), row.names = FALSE)

# Plot diversity by location for stool samples
col <- tableau.colors
ShannonLocation <- ggplot(subset(Diversity, Diversity$Location %in% c("Cincinnati", "Hangzhou") & 
                                 Diversity$SampleCollectionWeek == "Week.1" &
                                 Diversity$SampleType == "Stool"),
                        aes(x=Location, y=Shannon, fill=Location)) + 
  geom_boxplot(lwd=1, aes(color=factor(Location), fill=NA), outlier.size=3) +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + 
  scale_fill_manual(values=col) + 
  geom_point(size=4, aes(color=factor(Location))) + 
  ylab("Shannon Diversity Index : Stool\n") + 
  xlab(NULL) + 
  paramsBox() + 
  labs(caption="*p = 0.036") + 
  scale_color_tableau()

ggsave(filename=file.path(RESULTS_DIR, "ShannonStoolLocation.pdf"), plot=ShannonLocation, width=5, height=8, limitsize=FALSE)

# Statistical tests for diversity differences
StoolDiversity <- subset(Diversity, Diversity$Location %in% c("Cincinnati", "Hangzhou") & 
                         Diversity$SampleCollectionWeek == "Week.1" & 
                         Diversity$SampleType == "Stool")
stool_location_test <- pairwise.wilcox.test(StoolDiversity$Shannon, StoolDiversity$Location)
print("Stool diversity by location test:")
print(stool_location_test)

# Plot Shannon diversity by collection week
ShannonCollectionWeek <- ggplot(subset(Diversity, Diversity$SampleType %in% c("Axilla", "Groin", "Stool")),
                             aes(x=SampleCollectionWeek, y=Shannon, fill=TypeWeek)) + 
  geom_boxplot(lwd=1, aes(color=factor(TypeWeek), fill=NA), outlier.size=3) +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + 
  scale_fill_manual(values=col) + 
  scale_colour_tableau() +
  geom_point(size=4, aes(color=factor(TypeWeek))) + 
  ylab("Shannon Diversity Index\n") + 
  xlab(NULL) + 
  paramsBox() + 
  facet_grid(Location ~ SampleType) + 
  scale_fill_manual(values=col)

ggsave(filename=file.path(RESULTS_DIR, "ShannonCollectionWeek.pdf"), plot=ShannonCollectionWeek, width=14, height=12, limitsize=FALSE)

# Kruskal-Wallis and post-hoc tests for diversity by sample type
DiversityAxilla <- subset(Diversity, Diversity$SampleType == "Axilla")
DiversityGroin <- subset(Diversity, Diversity$SampleType == "Groin")
DiversityStool <- subset(Diversity, Diversity$SampleType == "Stool")

axilla_week_test <- pairwise.wilcox.test(DiversityAxilla$Shannon, DiversityAxilla$SampleCollectionWeek)
groin_week_test <- pairwise.wilcox.test(DiversityGroin$Shannon, DiversityGroin$SampleCollectionWeek)
Week1Stools$PostNatalAntibiotics`
SpeciesGroupMeans <- as.data.frame(t(SpeciesGroupMeans[,2:ncol(SpeciesGroupMeans)]))
SpeciesGroupMeans$Species <- row.names(SpeciesGroupMeans)

# Calculate overall means
OverallMeans <- as.data.frame(colMeans(Week1Stools[3:ncol(Week1Stools)]))
OverallMeans$Species <- row.names(OverallMeans)
names(OverallMeans) <- c("OverallMean", "Species")

# Merge group means and overall means
SpeciesGroupMeans <- merge(SpeciesGroupMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesGroupMeans$No.Infant.Abx <- as.numeric(SpeciesGroupMeans$No.Infant.Abx)
SpeciesGroupMeans$Infant.Abx <- as.numeric(SpeciesGroupMeans$Infant.Abx)

# Calculate fold changes and log ratios
SpeciesGroupMeans$FC <- foldchange(SpeciesGroupMeans$Infant.Abx, SpeciesGroupMeans$No.Infant.Abx)
SpeciesGroupMeans$logratio <- foldchange2logratio(SpeciesGroupMeans$FC)

# Combine all results
SpeciesGroupTable <- merge(SpeciesGroupMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesGroupTable <- merge(SpeciesGroupTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesGroupTable <- SpeciesGroupTable[order(SpeciesGroupTable$Unadjusted_p),]
SpeciesGroupTable <- SpeciesGroupTable[,c(1,8,4,2,3,5,9,10)]
names(SpeciesGroupTable) <- c("Species", "Abundance", "OverallMean", 
                            "No.Infant.Abx.Mean", "Infant.Abx.Mean", 
                            "Fold.Change", "p_Unadjusted", "FDR")
SpeciesGroupTable$Abundance <- SpeciesGroupTable$Abundance * 100
SpeciesGroupTable$Fold.Change[is.na(SpeciesGroupTable$Fold.Change)] <- 0
SpeciesGroupTable <- SpeciesGroupTable[order(SpeciesGroupTable$Fold.Change, decreasing = TRUE),]

# Filter for significant species
SigSpeciesGroupTable <- subset(SpeciesGroupTable, SpeciesGroupTable$p_Unadjusted < 0.2)

# Save results
write.csv(SigSpeciesGroupTable, file = file.path(RESULTS_DIR, "SignificantAntibioticsDifferences.csv"))
write.csv(SpeciesGroupTable, file = file.path(RESULTS_DIR, "AllAntibioticsSpeciesDifferences.csv"))

# Calculate effect sizes for antibiotics
metadata <- Week1Stools[, 1:2]
cts <- as.matrix(Week1Stools[, -(1:2)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- Week1Stools$PostNatalAntibiotics

# Compute effect sizes
abx_effects <- compute_ef_Abx(d = cts_l2, g = Week1Stools$PostNatalAntibiotics, min_shrink = 0.3)
abx_effects$Species <- abx_effects$Species  # Maintain column name

# Join with statistical results
abx_effects <- arrange(left_join(abx_effects, SigSpeciesGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

abx_effects <- subset(abx_effects, !is.na(abx_effects$FDR))
resWk1Abx <- abx_effects

# Add labels for which group has higher abundance
resWk1Abx$HigherIn <- ifelse(resWk1Abx$ef_shrunk > 0, "Infant.Abx", "No.Infant.Abx")
resWk1Abx$HigherIn <- factor(resWk1Abx$HigherIn, levels = c("Infant.Abx", "No.Infant.Abx"))

# Save effect size results
write.csv(resWk1Abx, file = file.path(RESULTS_DIR, "SignificantWk1AbxTableNew.csv"), 
        row.names = FALSE)

# Plot effect sizes
Week1Abx_Effect <- ggplot(filter(resWk1Abx, abs(resWk1Abx$ef_shrunk) > 0.5),
                       aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = HigherIn)) +
  geom_bar(colour="black", stat="identity") +
  guides(fill=guide_legend(override.aes=list(colour=NA))) + 
  coord_flip() +
  labs(x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More Abundant In:")) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) +
  scale_fill_manual(values = tableau.colors[c(1,2)])

ggsave(filename = file.path(RESULTS_DIR, "Week1Abx_Effect.pdf"), 
     plot = Week1Abx_Effect, width = 10, height = 8, limitsize = FALSE)

# Visualization of specific bacteria by antibiotics status
# Create dummy dataset for log plotting (replacing zeros)
DummySpeciesNR <- NICUSpeciesNR

# Bifidobacterium by antibiotics status
DummySpeciesNR$Bifidobacterium.longum[DummySpeciesNR$Bifidobacterium.longum == 0] <- 1
BifidoAbx <- ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                  aes(x=PostNatalAntibiotics, y=as.numeric(Bifidobacterium.longum), fill=PostNatalAntibiotics)) + 
  geom_boxplot(lwd=1, aes(color=factor(PostNatalAntibiotics), fill=NA), outlier.size=3) +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + 
  scale_fill_manual(values=col) + 
  scale_colour_manual(values=col) +
  geom_point(size=4, aes(color=factor(PostNatalAntibiotics))) + 
  xlab(NULL) + 
  facet_grid(rows=SampleType ~ SampleCollectionWeek + Location) +
  ylab("Bifidobacterium.longum \n") + 
  paramsBoxWide() + 
  scale_y_log10() + 
  xlab(NULL) + 
  scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))

ggsave(BifidoAbx, file=file.path(RESULTS_DIR, "BifidoAbx.pdf"), width=8, height=8, limitsize=FALSE)

# Staphylococcus aureus by antibiotics status
DummySpeciesNR$Staphylococcus.aureus[DummySpeciesNR$Staphylococcus.aureus == 0] <- 1
StaphAbx <- ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                aes(x=PostNatalAntibiotics, y=as.numeric(Staphylococcus.aureus), fill=PostNatalAntibiotics)) + 
  geom_boxplot(lwd=1, aes(color=factor(PostNatalAntibiotics), fill=NA), outlier.size=3) +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + 
  scale_fill_manual(values=col) + 
  scale_colour_manual(values=col) +
  geom_point(size=4, aes(color=factor(PostNatalAntibiotics))) + 
  xlab(NULL) + 
  facet_grid(rows=SampleType ~ SampleCollectionWeek) +
  ylab("Staphylococcus.aureus \n") + 
  paramsBoxWide() + 
  scale_y_log10() + 
  xlab(NULL) + 
  scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))

ggsave(StaphAbx, file=file.path(RESULTS_DIR, "StaphAbx.pdf"), width=12, height=10, limitsize=FALSE)

# Analyze location-specific species differences
# Staphylococcus aureus by location
StaphLocation <- ggplot(subset(DummySpeciesNR, DummySpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool")),
                     aes(x=Location, y=as.numeric(Staphylococcus.aureus), fill=Location)) + 
  geom_boxplot(lwd=1, aes(color=factor(Location), fill=NA), outlier.size=3) +
  stat_summary(fun=mean, geom="point", shape=5, size=8) + 
  scale_fill_manual(values=col) + 
  scale_colour_manual(values=col) +
  geom_point(size=4, aes(color=factor(Location))) + 
  xlab(NULL) + 
  ylab("Staphylococcus.aureus \n") + 
  paramsBoxWide() + 
  scale_y_log10() + 
  xlab(NULL) + 
  facet_grid(rows=. ~ SampleCollectionWeek) + 
  scale_color_tableau() +
  theme(axis.title.y = element_text(face = "italic"))

ggsave(StaphLocation, file=file.path(RESULTS_DIR, "StaphLocation.pdf"), width=8, height=6, limitsize=FALSE)

# PCA analysis for sample types and weeks
print("Performing PCA analysis for sample types and weeks...")

# Filter data for sample type analysis
AxillaGroinStoolSpecies <- subset(NICUSpeciesNR, NICUSpeciesNR$SampleType %in% c("Axilla", "Groin", "Stool"))
AxillaGroinStoolSpecies$TypeWeek <- paste(AxillaGroinStoolSpecies$SampleType, AxillaGroinStoolSpecies$SampleCollectionWeek, sep = "-")

# Prepare data for PCA
metadata <- AxillaGroinStoolSpecies[, 1:17]
cts <- as.matrix(AxillaGroinStoolSpecies[, -(1:17)])
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolSpecies$TypeWeek

# Perform PCA
col <- tableau.colors
sample_type_pca <- PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
sample_type_pca_plot <- fviz_pca_ind(sample_type_pca,
                                  geom.ind = "point", 
                                  pointsize = 1.5,
                                  point.alph = 0.1,
                                  title = "Unsupervised Principal Coordinate Analysis",
                                  subtitle = "Samples Colored by Sample Type and Infant Age",
                                  col.ind = grps, 
                                  addEllipses = TRUE, 
                                  ellipse.alpha = 0.2,
                                  ellipse.type = "confidence",
                                  ellipse.level = 0.95,
                                  legend.title = "Sample", 
                                  legend.size = 11,
                                  mean.point = TRUE,
                                  palette = col,
                                  axes.linetype = "blank"
)

# Save PCA plot
pdf(file.path(RESULTS_DIR, "SampleTypeWeekPCA_AllSamples.pdf")) 
print(sample_type_pca_plot)
dev.off() 

# Filter for no antibiotics samples
NoAbxSampleInd <- AxillaGroinStoolSpecies$Sample[AxillaGroinStoolSpecies$PostNatalAntibiotics == "No.Infant.Abx"]

# Create PCA plot for samples without antibiotics
noabx_pca_plot <- fviz_pca_ind(sample_type_pca,
                             geom.ind = "point", 
                             pointsize = 1.5,
                             point.alph = 0.1,
                             title = "Unsupervised Principal Coordinate Analysis",
                             subtitle = "Samples Colored by Sample Type and Infant Age : No Antibiotics",
                             col.ind = grps, 
                             addEllipses = TRUE, 
                             ellipse.alpha = 0.2,
                             ellipse.type = "confidence",
                             ellipse.level = 0.95,
                             legend.title = "Sample", 
                             legend.size = 11,
                             mean.point = TRUE,
                             palette = col,
                             axes.linetype = "blank",
                             select.ind = list(name = NoAbxSampleInd)
)

# Save PCA plot for samples without antibiotics
pdf(file.path(RESULTS_DIR, "SampleTypeWeekPCA_NOABX.pdf")) 
print(noabx_pca_plot)
dev.off() 

# Filter for antibiotic-exposed samples
AbxSampleInd <- AxillaGroinStoolSpecies$Sample[AxillaGroinStoolSpecies$PostNatalAntibiotics == "Infant.Abx"]

# Create PCA plot for antibiotic-exposed samples
abx_pca_plot <- fviz_pca_ind(sample_type_pca,
                           geom.ind = "point", 
                           pointsize = 1.5,
                           point.alph = 0.1,
                           title = "Unsupervised Principal Coordinate Analysis",
                           subtitle = "Samples Colored by Sample Type and Infant Age : Antibiotic Exposed",
                           col.ind = grps, 
                           addEllipses = TRUE, 
                           ellipse.alpha = 0.2,
                           ellipse.type = "confidence",
                           ellipse.level = 0.95,
                           legend.title = "Sample", 
                           legend.size = 11,
                           mean.point = TRUE,
                           palette = col,
                           axes.linetype = "blank",
                           select.ind = list(name = AbxSampleInd)
)

# Save PCA plot for antibiotic-exposed samples
pdf(file.path(RESULTS_DIR, "SampleTypeWeekPCA_ABX.pdf")) 
print(abx_pca_plot)
dev.off() 

# Analyze gestational age effects
print("Analyzing gestational age effects...")

# Create gestational age categories for analysis
AxillaGroinStoolGenusGA <- AxillaGroinStoolSpecies
AxillaGroinStoolGenusGA$GA2Wk <- ifelse(
  AxillaGroinStoolGenusGA$GestationTime >= 24 & AxillaGroinStoolGenusGA$GestationTime <= 25.99, "24-26", 
  ifelse(AxillaGroinStoolGenusGA$GestationTime >= 26 & AxillaGroinStoolGenusGA$GestationTime <= 27.99, "26-28",
         ifelse(AxillaGroinStoolGenusGA$GestationTime >= 28 & AxillaGroinStoolGenusGA$GestationTime <= 29.99, "28-30", 
                ifelse(AxillaGroinStoolGenusGA$GestationTime >= 30 & AxillaGroinStoolGenusGA$GestationTime <= 31.99, "30-32",
                       ifelse(AxillaGroinStoolGenusGA$GestationTime >= 32 & AxillaGroinStoolGenusGA$GestationTime <= 33.99, "32-34",
                              ifelse(AxillaGroinStoolGenusGA$GestationTime >= 34 & AxillaGroinStoolGenusGA$GestationTime <= 35.99, "34-36", 
                                     ifelse(AxillaGroinStoolGenusGA$GestationTime >= 36 & AxillaGroinStoolGenusGA$GestationTime <= 37.99, "36-38", 
                                            "NA")))))))

# Filter out NA values
AxillaGroinStoolGenusGA <- subset(AxillaGroinStoolGenusGA, !is.na(AxillaGroinStoolGenusGA$GA2Wk))

# Create combined factor for analysis
AxillaGroinStoolGenusGA$AgeTypeWeek <- paste(
  AxillaGroinStoolGenusGA$SampleType, 
  AxillaGroinStoolGenusGA$SampleCollectionWeek,
  AxillaGroinStoolGenusGA$GA2Wk, 
  sep = "-"
)

# Filter for no antibiotics samples and stool type
AxillaGroinStoolGenusGANoAbx <- subset(
  AxillaGroinStoolGenusGA, 
  AxillaGroinStoolGenusGA$PostNatalAntibiotics == "No.Infant.Abx" &
  AxillaGroinStoolGenusGA$SampleType == "Stool"
)

# Prepare data for PCA
metadata <- AxillaGroinStoolGenusGANoAbx[, 1:20]  # Adjust column range as needed
cts <- as.matrix(AxillaGroinStoolGenusGANoAbx[, -(1:20)])  # Adjust column range as needed
rownames(cts) <- metadata$Sample
cts_l2 <- glog2(cts)
grps <- AxillaGroinStoolGenusGANoAbx$AgeTypeWeek

# Perform PCA for gestational age analysis
col <- rep(tableau.colors, 2)
ga_pca <- PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
ga_pca_plot <- fviz_pca_ind(ga_pca,
                         geom.ind = "point", 
                         pointsize = 1.5,
                         point.alph = 0.1,
                         title = "Unsupervised Principal Coordinate Analysis",
                         subtitle = "Samples Colored by Gestational Age and Collection Week : No Antibiotics",
                         col.ind = grps, 
                         addEllipses = TRUE, 
                         ellipse.alpha = 0.2,
                         ellipse.type = "confidence",
                         ellipse.level = 0.95,
                         legend.title = "Sample", 
                         legend.size = 11,
                         mean.point = TRUE,
                         palette = col,
                         axes.linetype = "blank"
)

# Save PCA plot for gestational age analysis
pdf(file.path(RESULTS_DIR, "GestationalAgePCA.pdf"))
print(ga_pca_plot)
dev.off()

# Perform MRPP tests for gestational age groups
# Filter samples by gestational age range
ga_28_30 <- subset(AxillaGroinStoolGenusGANoAbx, 
                AxillaGroinStoolGenusGANoAbx$GestationTime >= 28 & 
                AxillaGroinStoolGenusGANoAbx$GestationTime <= 29.99 &
                AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")

ga_30_32 <- subset(AxillaGroinStoolGenusGANoAbx, 
                AxillaGroinStoolGenusGANoAbx$GestationTime >= 30 & 
                AxillaGroinStoolGenusGANoAbx$GestationTime <= 31.99 &
                AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")

ga_32_34 <- subset(AxillaGroinStoolGenusGANoAbx, 
                AxillaGroinStoolGenusGANoAbx$GestationTime >= 32 & 
                AxillaGroinStoolGenusGANoAbx$GestationTime <= 33.99 &
                AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")

ga_34_36 <- subset(AxillaGroinStoolGenusGANoAbx, 
                AxillaGroinStoolGenusGANoAbx$GestationTime >= 34 & 
                AxillaGroinStoolGenusGANoAbx$GestationTime <= 35.99 &
                AxillaGroinStoolGenusGANoAbx$SampleType == "Stool")

# Perform MRPP tests
data_cols <- 21:ncol(AxillaGroinStoolGenusGANoAbx)  # Adjust as needed for your data
mrpp_28_30 <- mrpp(ga_28_30[, data_cols], ga_28_30$SampleCollectionWeek, distance = "bray")
mrpp_30_32 <- mrpp(ga_30_32[, data_cols], ga_30_32$SampleCollectionWeek, distance = "bray")
mrpp_32_34 <- mrpp(ga_32_34[, data_cols], ga_32_34$SampleCollectionWeek, distance = "bray")
mrpp_34_36 <- mrpp(ga_34_36[, data_cols], ga_34_36$SampleCollectionWeek, distance = "bray")

print("MRPP results for gestational age groups:")
print("28-30 weeks:")
print(mrpp_28_30)
print("30-32 weeks:")
print(mrpp_30_32)
print("32-34 weeks:")
print(mrpp_32_34)
print("34-36 weeks:")
print(mrpp_34_36)

# Save the workspace
save.image(file = file.path(RESULTS_DIR, "NICU_AnalysisComplete.RData"))

print("Analysis complete!")
print(paste("Results saved to:", RESULTS_DIR))
