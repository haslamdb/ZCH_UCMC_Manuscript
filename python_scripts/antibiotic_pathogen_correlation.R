#!/usr/bin/env Rscript
# Correlate antibiotic exposures with BSI pathogen abundances
# Analysis plan: antibiotic_BSI_analysis_plan.md

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript")

# Create output directory
dir.create("antibiotic_pathogen_analysis", showWarnings = FALSE)

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("ANTIBIOTIC-PATHOGEN CORRELATION ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n\n")

# ============================================================================
# STEP 1: Load antibiotic data and create cumulative week 3 exposures
# ============================================================================

cat("STEP 1: Loading antibiotic data and creating cumulative week 3 exposures\n")
cat(paste(rep("-", 80), collapse = ""))
cat("\n")

antibiotics <- read.csv("antibiotics_combined.csv")

cat(sprintf("Loaded antibiotic data: %d patients\n", nrow(antibiotics)))
cat(sprintf("  - UCMC: %d patients\n", sum(antibiotics$Location == "UCMC")))
cat(sprintf("  - ZCH: %d patients\n", sum(antibiotics$Location == "ZCH")))

# FIX: Normalize Subject IDs (N1-N9 -> N01-N09 for ZCH subjects)
# This ensures proper matching with microbiome data
antibiotics$Subject <- gsub("^N(\\d)$", "N0\\1", antibiotics$Subject)
cat("\n✓ Normalized Subject IDs (N1->N01, etc.)\n")

# Get antibiotic column names
abx_cols <- grep("^(Ampicillin|Penicillin|Nafcillin|Cefotaxime|Moxalactam|Ceftazidime|Ceftriaxone|Cefepime|Gentamicin|Tobramycin|Vancomycin|Fluconazole|Amphotericin|Meropenem|Azithromycin|Piperacillin)",
                 colnames(antibiotics), value = TRUE)

# Separate week 1 and week 2 (which is actually week 3) columns
w1_cols <- grep("_w1$", abx_cols, value = TRUE)
w2_cols <- grep("_w2$", abx_cols, value = TRUE)

cat(sprintf("\nFound %d week 1 antibiotics and %d week 3 antibiotics\n",
            length(w1_cols), length(w2_cols)))

# Create cumulative week 3 by summing week 1 + week 2 (actually week 3)
antibiotics_expanded <- antibiotics

for (i in seq_along(w1_cols)) {
  w1_col <- w1_cols[i]
  # Get corresponding w2 column by replacing _w1 with _w2
  w2_col <- gsub("_w1$", "_w2", w1_col)

  if (w2_col %in% colnames(antibiotics)) {
    # Create cumulative column
    cumulative_col <- gsub("_w1$", "_cumulative_w3", w1_col)
    antibiotics_expanded[[cumulative_col]] <- antibiotics[[w1_col]] + antibiotics[[w2_col]]

    cat(sprintf("  Created %s = %s + %s\n", cumulative_col, w1_col, w2_col))
  }
}

cat("\n✓ Created cumulative week 3 antibiotic exposures\n")

# ============================================================================
# STEP 2: Load microbiome data and extract pathogen abundances
# ============================================================================

cat("\n")
cat("STEP 2: Loading microbiome data and extracting pathogen abundances\n")
cat(paste(rep("-", 80), collapse = ""))
cat("\n")

# Load R workspace
load("/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514")

cat("✓ Loaded R workspace with microbiome data\n")

# 8 BSI pathogens from Figure 4
bsi_pathogens <- c(
  "Staphylococcus.aureus",
  "Klebsiella.pneumoniae",
  "Klebsiella.oxytoca",
  "Enterococcus.faecium",
  "Enterococcus.faecalis",
  "Serratia.marcescens",
  "Escherichia.coli",
  "Streptococcus.pyogenes"
)

# Extract relevant data
# Filter to Week 1 and Week 3 samples
# PatientID is the Subject column (N01, N02, etc.)
microbiome_data <- NICUSpeciesNR %>%
  filter(SampleCollectionWeek %in% c("Week.1", "Week.3")) %>%
  select(Sample, PatientID, Location, SampleType, SampleCollectionWeek, all_of(bsi_pathogens)) %>%
  rename(Subject = PatientID)

# Rename locations to match antibiotic data
microbiome_data <- microbiome_data %>%
  mutate(Location = recode(Location,
                          "Cincinnati" = "UCMC",
                          "Hangzhou" = "ZCH"))

cat(sprintf("\nExtracted pathogen data: %d samples\n", nrow(microbiome_data)))
cat(sprintf("  - Week 1: %d samples\n", sum(microbiome_data$SampleCollectionWeek == "Week.1")))
cat(sprintf("  - Week 3: %d samples\n", sum(microbiome_data$SampleCollectionWeek == "Week.3")))
cat(sprintf("  - Sample types: %s\n", paste(unique(microbiome_data$SampleType), collapse = ", ")))

# For each Subject x Week, average across sample types (Axilla, Groin, Stool)
# This gives one value per subject per week
pathogen_avg <- microbiome_data %>%
  group_by(Subject, Location, SampleCollectionWeek) %>%
  summarise(across(all_of(bsi_pathogens), ~mean(.x, na.rm = TRUE)), .groups = 'drop')

# Create separate datasets for Week 1 and Week 3
pathogen_w1 <- pathogen_avg %>%
  filter(SampleCollectionWeek == "Week.1") %>%
  select(-SampleCollectionWeek)

pathogen_w3 <- pathogen_avg %>%
  filter(SampleCollectionWeek == "Week.3") %>%
  select(-SampleCollectionWeek)

# Rename pathogen columns to indicate week
colnames(pathogen_w1)[3:ncol(pathogen_w1)] <- paste0(colnames(pathogen_w1)[3:ncol(pathogen_w1)], "_w1")
colnames(pathogen_w3)[3:ncol(pathogen_w3)] <- paste0(colnames(pathogen_w3)[3:ncol(pathogen_w3)], "_w3")

cat(sprintf("\n✓ Averaged pathogen abundances across sample types\n"))
cat(sprintf("  - Week 1: %d subjects\n", nrow(pathogen_w1)))
cat(sprintf("  - Week 3: %d subjects\n", nrow(pathogen_w3)))

# ============================================================================
# STEP 3: Merge antibiotic and pathogen data
# ============================================================================

cat("\n")
cat("STEP 3: Merging antibiotic and pathogen datasets\n")
cat(paste(rep("-", 80), collapse = ""))
cat("\n")

# Merge with antibiotic data
merged_data <- antibiotics_expanded %>%
  left_join(pathogen_w1, by = c("Subject", "Location")) %>%
  left_join(pathogen_w3, by = c("Subject", "Location"))

cat(sprintf("\n✓ Merged datasets: %d subjects\n", nrow(merged_data)))

# Check how many subjects have both antibiotic and pathogen data
complete_w1 <- sum(!is.na(merged_data[[paste0(bsi_pathogens[1], "_w1")]]))
complete_w3 <- sum(!is.na(merged_data[[paste0(bsi_pathogens[1], "_w3")]]))

cat(sprintf("  - Subjects with week 1 pathogen data: %d\n", complete_w1))
cat(sprintf("  - Subjects with week 3 pathogen data: %d\n", complete_w3))

# Save merged dataset
write.csv(merged_data, "antibiotic_pathogen_analysis/antibiotic_pathogen_merged.csv", row.names = FALSE)
cat("\n✓ Saved merged dataset: antibiotic_pathogen_analysis/antibiotic_pathogen_merged.csv\n")

# ============================================================================
# STEP 4: Calculate Spearman correlations
# ============================================================================

cat("\n")
cat("STEP 4: Calculating Spearman correlations\n")
cat(paste(rep("-", 80), collapse = ""))
cat("\n")

# Function to calculate correlation with p-value
calc_correlation <- function(x, y) {
  # Remove pairs with missing values
  valid <- !is.na(x) & !is.na(y)

  if (sum(valid) < 5) {  # Need at least 5 pairs
    return(list(rho = NA, p.value = NA, n = sum(valid)))
  }

  test <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
  return(list(rho = test$estimate, p.value = test$p.value, n = sum(valid)))
}

# Week 1 correlations: week 1 antibiotics vs week 1 pathogens
cat("\nCalculating Week 1 correlations...\n")
w1_results <- data.frame()

for (abx in w1_cols) {
  for (pathogen in bsi_pathogens) {
    pathogen_col <- paste0(pathogen, "_w1")

    if (pathogen_col %in% colnames(merged_data)) {
      cor_result <- calc_correlation(merged_data[[abx]], merged_data[[pathogen_col]])

      w1_results <- rbind(w1_results, data.frame(
        Antibiotic = abx,
        Pathogen = pathogen,
        Week = "Week1",
        N = cor_result$n,
        Rho = cor_result$rho,
        P_value = cor_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Week 3 correlations: cumulative week 3 antibiotics vs week 3 pathogens
cat("Calculating Week 3 (cumulative) correlations...\n")
w3_results <- data.frame()

cumulative_w3_cols <- grep("_cumulative_w3$", colnames(antibiotics_expanded), value = TRUE)

for (abx in cumulative_w3_cols) {
  for (pathogen in bsi_pathogens) {
    pathogen_col <- paste0(pathogen, "_w3")

    if (pathogen_col %in% colnames(merged_data)) {
      cor_result <- calc_correlation(merged_data[[abx]], merged_data[[pathogen_col]])

      w3_results <- rbind(w3_results, data.frame(
        Antibiotic = abx,
        Pathogen = pathogen,
        Week = "Week3_Cumulative",
        N = cor_result$n,
        Rho = cor_result$rho,
        P_value = cor_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Combine results
all_correlations <- rbind(w1_results, w3_results)

# Apply FDR correction within each week
all_correlations <- all_correlations %>%
  group_by(Week) %>%
  mutate(P_adjusted = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(Significant = ifelse(P_adjusted < 0.05, "Yes", "No"))

cat(sprintf("\n✓ Calculated %d correlations\n", nrow(all_correlations)))
cat(sprintf("  - Week 1: %d correlations\n", sum(all_correlations$Week == "Week1")))
cat(sprintf("  - Week 3 (cumulative): %d correlations\n", sum(all_correlations$Week == "Week3_Cumulative")))

# Save all correlations
write.csv(all_correlations, "antibiotic_pathogen_analysis/antibiotic_pathogen_correlations.csv", row.names = FALSE)
cat("\n✓ Saved all correlations: antibiotic_pathogen_analysis/antibiotic_pathogen_correlations.csv\n")

# Save significant correlations
significant_cors <- all_correlations %>% filter(P_adjusted < 0.05) %>% arrange(P_adjusted)
write.csv(significant_cors, "antibiotic_pathogen_analysis/antibiotic_pathogen_significant.csv", row.names = FALSE)
cat(sprintf("✓ Saved %d significant correlations (FDR < 0.05): antibiotic_pathogen_analysis/antibiotic_pathogen_significant.csv\n",
            nrow(significant_cors)))

# Print summary
cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("CORRELATION SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat(sprintf("\nTotal correlations tested: %d\n", nrow(all_correlations)))
cat(sprintf("Significant (FDR < 0.05): %d (%.1f%%)\n",
            nrow(significant_cors), 100*nrow(significant_cors)/nrow(all_correlations)))

if (nrow(significant_cors) > 0) {
  cat("\nTop 10 significant correlations:\n")
  top_cors <- head(significant_cors, 10)
  for (i in 1:nrow(top_cors)) {
    cat(sprintf("  %s - %s (%s): rho=%.3f, p=%.2e, FDR=%.2e, n=%d\n",
                gsub("_w1|_cumulative_w3", "", top_cors$Antibiotic[i]),
                gsub("\\.", " ", top_cors$Pathogen[i]),
                top_cors$Week[i],
                top_cors$Rho[i],
                top_cors$P_value[i],
                top_cors$P_adjusted[i],
                top_cors$N[i]))
  }
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("\nGenerated files:\n")
cat("  - antibiotic_pathogen_analysis/antibiotic_pathogen_merged.csv\n")
cat("  - antibiotic_pathogen_analysis/antibiotic_pathogen_correlations.csv\n")
cat("  - antibiotic_pathogen_analysis/antibiotic_pathogen_significant.csv\n")
cat("\nNext steps:\n")
cat("  - Review significant correlations\n")
cat("  - Generate correlation heatmap\n")
cat("  - Create scatter plots for top associations\n")
cat("\n")
