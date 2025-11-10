#!/usr/bin/env Rscript
# Correlate antibiotic exposures with BSI pathogen abundances
# Extended analysis: stratified by body site (Axilla, Groin, Stool)
# Added: Candida parapsilosis

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
cat("ANTIBIOTIC-PATHOGEN CORRELATION ANALYSIS (BY BODY SITE)\n")
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

# 8 BSI pathogens from Figure 4 + Candida parapsilosis
bsi_pathogens <- c(
  "Staphylococcus.aureus",
  "Klebsiella.pneumoniae",
  "Klebsiella.oxytoca",
  "Enterococcus.faecium",
  "Enterococcus.faecalis",
  "Serratia.marcescens",
  "Escherichia.coli",
  "Streptococcus.pyogenes",
  "Candida.parapsilosis"
)

cat(sprintf("\nAnalyzing %d pathogens (8 BSI + Candida parapsilosis)\n", length(bsi_pathogens)))

# Extract relevant data
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

# Create datasets for each body site separately (NOT averaged)
# This gives one value per subject per week per site
pathogen_by_site <- microbiome_data %>%
  select(Subject, Location, SampleType, SampleCollectionWeek, all_of(bsi_pathogens))

# Create separate datasets for Week 1 and Week 3, keeping site information
pathogen_w1_by_site <- pathogen_by_site %>%
  filter(SampleCollectionWeek == "Week.1") %>%
  select(-SampleCollectionWeek)

pathogen_w3_by_site <- pathogen_by_site %>%
  filter(SampleCollectionWeek == "Week.3") %>%
  select(-SampleCollectionWeek)

# Rename pathogen columns to indicate week
colnames(pathogen_w1_by_site)[4:ncol(pathogen_w1_by_site)] <-
  paste0(colnames(pathogen_w1_by_site)[4:ncol(pathogen_w1_by_site)], "_w1")
colnames(pathogen_w3_by_site)[4:ncol(pathogen_w3_by_site)] <-
  paste0(colnames(pathogen_w3_by_site)[4:ncol(pathogen_w3_by_site)], "_w3")

cat(sprintf("\n✓ Organized pathogen data by body site\n"))
cat(sprintf("  - Week 1: %d samples across 3 sites\n", nrow(pathogen_w1_by_site)))
cat(sprintf("  - Week 3: %d samples across 3 sites\n", nrow(pathogen_w3_by_site)))

# ============================================================================
# STEP 3: Merge antibiotic and pathogen data (for each site separately)
# ============================================================================

cat("\n")
cat("STEP 3: Merging antibiotic and pathogen datasets by body site\n")
cat(paste(rep("-", 80), collapse = ""))
cat("\n")

# Merge with antibiotic data - one row per subject per site per week
merged_w1 <- antibiotics_expanded %>%
  inner_join(pathogen_w1_by_site, by = c("Subject", "Location"))

merged_w3 <- antibiotics_expanded %>%
  inner_join(pathogen_w3_by_site, by = c("Subject", "Location"))

# Combine week 1 and week 3
merged_data_by_site <- bind_rows(
  merged_w1 %>% mutate(Week = "Week1"),
  merged_w3 %>% mutate(Week = "Week3")
)

cat(sprintf("\n✓ Merged datasets by body site:\n"))
for (site in c("Axilla", "Groin", "Stool")) {
  n_w1 <- sum(merged_w1$SampleType == site)
  n_w3 <- sum(merged_w3$SampleType == site)
  cat(sprintf("  - %s: %d week 1, %d week 3 samples\n", site, n_w1, n_w3))
}

# Save merged dataset
write.csv(merged_data_by_site,
          "antibiotic_pathogen_analysis/antibiotic_pathogen_merged_by_site.csv",
          row.names = FALSE)
cat("\n✓ Saved: antibiotic_pathogen_analysis/antibiotic_pathogen_merged_by_site.csv\n")

# ============================================================================
# STEP 4: Calculate Spearman correlations (stratified by site)
# ============================================================================

cat("\n")
cat("STEP 4: Calculating Spearman correlations by body site\n")
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

all_results <- data.frame()

# Analyze each body site separately
for (site in c("Axilla", "Groin", "Stool")) {
  cat(sprintf("\n--- %s ---\n", site))

  site_w1 <- merged_w1 %>% filter(SampleType == site)
  site_w3 <- merged_w3 %>% filter(SampleType == site)

  # Week 1 correlations
  cat(sprintf("  Week 1 correlations (n=%d)...\n", nrow(site_w1)))
  for (abx in w1_cols) {
    for (pathogen in bsi_pathogens) {
      pathogen_col <- paste0(pathogen, "_w1")

      if (pathogen_col %in% colnames(site_w1)) {
        cor_result <- calc_correlation(site_w1[[abx]], site_w1[[pathogen_col]])

        all_results <- rbind(all_results, data.frame(
          SampleType = site,
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

  # Week 3 correlations (cumulative antibiotics)
  cat(sprintf("  Week 3 cumulative correlations (n=%d)...\n", nrow(site_w3)))
  cumulative_w3_cols <- grep("_cumulative_w3$", colnames(antibiotics_expanded), value = TRUE)

  for (abx in cumulative_w3_cols) {
    for (pathogen in bsi_pathogens) {
      pathogen_col <- paste0(pathogen, "_w3")

      if (pathogen_col %in% colnames(site_w3)) {
        cor_result <- calc_correlation(site_w3[[abx]], site_w3[[pathogen_col]])

        all_results <- rbind(all_results, data.frame(
          SampleType = site,
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
}

# Apply FDR correction within each site and week combination
all_results <- all_results %>%
  group_by(SampleType, Week) %>%
  mutate(P_adjusted = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(Significant = ifelse(P_adjusted < 0.05, "Yes", "No"))

cat(sprintf("\n✓ Calculated %d correlations across all sites\n", nrow(all_results)))

# Summary by site
for (site in c("Axilla", "Groin", "Stool")) {
  site_results <- all_results %>% filter(SampleType == site)
  n_sig <- sum(site_results$P_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("  - %s: %d correlations, %d significant (FDR<0.05)\n",
              site, nrow(site_results), n_sig))
}

# Save all correlations
write.csv(all_results,
          "antibiotic_pathogen_analysis/antibiotic_pathogen_correlations_by_site.csv",
          row.names = FALSE)
cat("\n✓ Saved: antibiotic_pathogen_analysis/antibiotic_pathogen_correlations_by_site.csv\n")

# Save significant correlations
significant_cors <- all_results %>%
  filter(P_adjusted < 0.05) %>%
  arrange(P_adjusted)

write.csv(significant_cors,
          "antibiotic_pathogen_analysis/antibiotic_pathogen_significant_by_site.csv",
          row.names = FALSE)
cat(sprintf("✓ Saved %d significant correlations: antibiotic_pathogen_analysis/antibiotic_pathogen_significant_by_site.csv\n",
            nrow(significant_cors)))

# ============================================================================
# STEP 5: Summary statistics
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("CORRELATION SUMMARY (BY BODY SITE)\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat(sprintf("\nTotal correlations tested: %d\n", nrow(all_results)))
cat(sprintf("Significant (FDR < 0.05): %d (%.1f%%)\n",
            nrow(significant_cors), 100*nrow(significant_cors)/nrow(all_results)))

if (nrow(significant_cors) > 0) {
  cat("\n--- Top 20 significant correlations ---\n")
  top_cors <- head(significant_cors, 20)
  for (i in 1:nrow(top_cors)) {
    cat(sprintf("%2d. [%s] %s - %s (%s): rho=%.3f, FDR=%.2e, n=%d\n",
                i,
                top_cors$SampleType[i],
                gsub("_w1|_cumulative_w3", "", top_cors$Antibiotic[i]),
                gsub("\\.", " ", top_cors$Pathogen[i]),
                top_cors$Week[i],
                top_cors$Rho[i],
                top_cors$P_adjusted[i],
                top_cors$N[i]))
  }

  # Summary by site and week
  cat("\n--- Significant correlations by site and week ---\n")
  summary_table <- significant_cors %>%
    group_by(SampleType, Week) %>%
    summarise(N_significant = n(), .groups = 'drop') %>%
    arrange(SampleType, Week)
  print(summary_table)
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("\nGenerated files:\n")
cat("  - antibiotic_pathogen_analysis/antibiotic_pathogen_merged_by_site.csv\n")
cat("  - antibiotic_pathogen_analysis/antibiotic_pathogen_correlations_by_site.csv\n")
cat("  - antibiotic_pathogen_analysis/antibiotic_pathogen_significant_by_site.csv\n")
cat("\nNext step: Generate site-specific visualizations\n")
cat("\n")
