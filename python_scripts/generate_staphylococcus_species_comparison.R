#!/usr/bin/env Rscript
# Compare Staphylococcus species abundance between UCMC and ZCH
# By body site (Axilla, Groin, Stool) and time point (Week 1, Week 3)
# Similar to Figure 4 style

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript/python_scripts")

# Create output directory
dir.create("revision_figures", showWarnings = FALSE)

# Load the saved R workspace
cat("Loading workspace...\n")
load("/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514")
cat("Workspace loaded successfully!\n")

cat("\n=== Staphylococcus Species Comparison Analysis ===\n")

# Key Staphylococcus species to analyze
staph_species <- c(
  "Staphylococcus.aureus",      # Already in Figure 4, include for comparison
  "Staphylococcus.epidermidis", # Most common CoNS
  "Staphylococcus.capitis",     # Common in NICU
  "Staphylococcus.hominis",     # Common skin colonizer
  "Staphylococcus.haemolyticus",# Common CoNS
  "Staphylococcus.warneri"      # Less common but relevant
)

cat("\nAnalyzing", length(staph_species), "Staphylococcus species:\n")
for (sp in staph_species) {
  cat("  -", gsub("\\.", " ", sp), "\n")
}

# Prepare data
cat("\n=== Preparing data ===\n")
staph_data <- NICUSpeciesNR %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         SampleCollectionWeek %in% c("Week.1", "Week.3"))

cat(sprintf("Data prepared: %d samples\n", nrow(staph_data)))
cat(sprintf("  Cincinnati: %d samples\n", sum(staph_data$Location == "Cincinnati")))
cat(sprintf("  Hangzhou: %d samples\n", sum(staph_data$Location == "Hangzhou")))

# Function to create species-specific plot with faceting
create_staph_plot <- function(data, species, panel_label) {

  # Check if species exists in data
  if (!species %in% colnames(data)) {
    cat("WARNING:", species, "not found in data. Skipping.\n")
    return(NULL)
  }

  # Select relevant data
  plot_data <- data %>%
    select(Sample, Location, SampleType, SampleCollectionWeek, all_of(species)) %>%
    rename(abundance = all_of(species)) %>%
    filter(!is.na(abundance)) %>%
    mutate(abundance_log = ifelse(abundance == 0, 1, abundance))  # Replace 0 with 1 for log scale

  # Calculate sample sizes and positioning for labels
  sample_counts <- plot_data %>%
    group_by(Location, SampleType, SampleCollectionWeek) %>%
    summarise(
      n = n(),
      median_abundance = median(abundance_log, na.rm = TRUE),
      mean_abundance = mean(abundance_log, na.rm = TRUE),
      n_nonzero = sum(abundance > 0),
      prevalence = n_nonzero / n * 100,
      max_abundance = max(abundance_log, na.rm = TRUE),
      .groups = 'drop'
    )

  # Calculate y position for sample count labels
  sample_counts <- sample_counts %>%
    group_by(SampleType, SampleCollectionWeek) %>%
    mutate(y_pos = max(max_abundance, na.rm = TRUE) * 1.3) %>%
    ungroup()

  cat("\n", gsub("\\.", " ", species), "- Sample counts:\n")
  print(sample_counts)

  # Calculate p-values for location comparisons at each body site and week
  pvalue_results <- data.frame()

  for (week in c("Week.1", "Week.3")) {
    for (site in c("Axilla", "Groin", "Stool")) {
      subset_data <- plot_data %>%
        filter(SampleCollectionWeek == week, SampleType == site)

      if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
        # Use full data (including zeros) for comparison
        ucmc_data <- subset_data %>% filter(Location == "Cincinnati") %>% pull(abundance)
        zch_data <- subset_data %>% filter(Location == "Hangzhou") %>% pull(abundance)

        test_result <- tryCatch({
          wilcox.test(ucmc_data, zch_data)
        }, error = function(e) {
          list(p.value = NA)
        })

        p_value <- test_result$p.value
        max_y <- max(subset_data$abundance_log, na.rm = TRUE)

        p_label <- ifelse(is.na(p_value), "",
                  ifelse(p_value < 0.001, "p<0.001***",
                  ifelse(p_value < 0.01, sprintf("p=%.3f**", p_value),
                  ifelse(p_value < 0.05, sprintf("p=%.3f*", p_value),
                         sprintf("p=%.2f", p_value)))))

        pvalue_results <- rbind(pvalue_results, data.frame(
          SampleType = site,
          SampleCollectionWeek = week,
          x = 1.5,  # Position between Cincinnati and Hangzhou
          y_pos = max_y * 1.5,
          p_value = p_value,
          p_label = p_label
        ))
      }
    }
  }

  cat("\nP-values:\n")
  print(pvalue_results)

  # Create plot with faceting (SampleType ~ SampleCollectionWeek)
  # X-axis: Location (Cincinnati vs Hangzhou)
  p <- ggplot(plot_data,
              aes(x = Location, y = abundance_log)) +
    geom_boxplot(lwd = 1, aes(color = factor(Location), fill = NA),
                 outlier.shape = NA) +  # Don't show outliers, we'll add points
    stat_summary(fun = mean, geom = "point", shape = 5, size = 4.8,
                data = plot_data) +
    scale_colour_manual(values = c("Cincinnati" = "#E69F00", "Hangzhou" = "#56B4E9")) +
    geom_point(size = 2.4, aes(color = factor(Location))) +
    # Add sample count labels
    geom_text(data = sample_counts,
              aes(x = Location, y = y_pos, label = paste0("n=", n)),
              size = 7,
              fontface = "bold",
              inherit.aes = FALSE) +
    # Add p-value labels
    geom_text(data = pvalue_results,
              aes(x = x, y = y_pos, label = p_label),
              size = 7,
              fontface = "bold",
              inherit.aes = FALSE) +
    facet_grid(SampleType ~ SampleCollectionWeek) +
    scale_y_log10() +
    ylab(paste0(gsub("\\.", " ", species), "\n(log10 abundance)")) +
    xlab(NULL) +
    ggtitle(paste0(panel_label, ". ", gsub("\\.", " ", species))) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 20, face = "bold.italic"),
      strip.text = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 22, face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  return(list(plot = p, counts = sample_counts, pvalues = pvalue_results))
}

# Generate plots for each Staphylococcus species
cat("\n=== Generating individual plots ===\n")

staph_plots <- list()
all_counts <- list()
all_pvalues <- list()
panel_labels <- letters[1:length(staph_species)]

for (i in seq_along(staph_species)) {
  species <- staph_species[i]
  panel_label <- panel_labels[i]

  result <- create_staph_plot(staph_data, species, panel_label)

  if (!is.null(result)) {
    staph_plots[[species]] <- result$plot
    all_counts[[species]] <- result$counts
    all_pvalues[[species]] <- result$pvalues

    # Save individual plot
    clean_name <- gsub("\\.", "_", species)
    ggsave(paste0("revision_figures/", clean_name, "_comparison.pdf"),
           plot = result$plot, width = 12, height = 10)

    cat("✓ Saved", species, "comparison plot\n")
  }
}

# Combine all sample counts into one table
cat("\n=== Saving statistics tables ===\n")
combined_counts <- bind_rows(all_counts, .id = "Species")
write.csv(combined_counts,
          "revision_figures/Staphylococcus_species_sample_counts.csv",
          row.names = FALSE)
cat("✓ Saved sample counts table\n")

# Combine all p-values
combined_pvalues <- bind_rows(all_pvalues, .id = "Species")
write.csv(combined_pvalues,
          "revision_figures/Staphylococcus_species_pvalues.csv",
          row.names = FALSE)
cat("✓ Saved p-values table\n")

# ============================================================================
# Create comprehensive multi-panel figure
# ============================================================================

cat("\n=== Creating comprehensive multi-panel figure ===\n")

# Simplify plots for combined figure
create_staph_plot_simplified <- function(data, species, panel_label) {

  if (!species %in% colnames(data)) {
    return(NULL)
  }

  plot_data <- data %>%
    select(Sample, Location, SampleType, SampleCollectionWeek, all_of(species)) %>%
    rename(abundance = all_of(species)) %>%
    filter(!is.na(abundance)) %>%
    mutate(abundance_log = ifelse(abundance == 0, 1, abundance))

  # Calculate sample sizes and positioning
  sample_counts <- plot_data %>%
    group_by(Location, SampleType, SampleCollectionWeek) %>%
    summarise(n = n(),
              max_abundance = max(abundance_log, na.rm = TRUE),
              .groups = 'drop') %>%
    group_by(SampleType, SampleCollectionWeek) %>%
    mutate(y_pos = max(max_abundance, na.rm = TRUE) * 1.3) %>%
    ungroup()

  # Calculate p-values for simplified plot
  pvalue_results <- data.frame()

  for (week in c("Week.1", "Week.3")) {
    for (site in c("Axilla", "Groin", "Stool")) {
      subset_data <- plot_data %>%
        filter(SampleCollectionWeek == week, SampleType == site)

      if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
        ucmc_data <- subset_data %>% filter(Location == "Cincinnati") %>% pull(abundance)
        zch_data <- subset_data %>% filter(Location == "Hangzhou") %>% pull(abundance)

        test_result <- tryCatch({
          wilcox.test(ucmc_data, zch_data)
        }, error = function(e) {
          list(p.value = NA)
        })

        p_value <- test_result$p.value
        max_y <- max(subset_data$abundance_log, na.rm = TRUE)

        p_label <- ifelse(is.na(p_value), "",
                  ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                  ifelse(p_value < 0.05, "*", ""))))

        if (p_label != "") {
          pvalue_results <- rbind(pvalue_results, data.frame(
            SampleType = site,
            SampleCollectionWeek = week,
            x = 1.5,
            y_pos = max_y * 1.5,
            p_label = p_label
          ))
        }
      }
    }
  }

  p <- ggplot(plot_data,
              aes(x = Location, y = abundance_log)) +
    geom_boxplot(lwd = 0.7, aes(color = Location, fill = NA),
                 outlier.shape = NA, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, color = "black") +
    geom_point(size = 1.5, aes(color = Location), alpha = 0.5) +
    scale_colour_manual(values = c("Cincinnati" = "#E69F00", "Hangzhou" = "#56B4E9")) +
    # Add sample count labels
    geom_text(data = sample_counts,
              aes(x = Location, y = y_pos, label = paste0("n=", n)),
              size = 3,
              fontface = "bold",
              inherit.aes = FALSE) +
    # Add p-value labels (only if significant)
    {if (nrow(pvalue_results) > 0)
      geom_text(data = pvalue_results,
                aes(x = x, y = y_pos, label = p_label),
                size = 3.5,
                fontface = "bold",
                inherit.aes = FALSE)
    } +
    facet_grid(SampleType ~ SampleCollectionWeek) +
    scale_y_log10() +
    ylab(NULL) +
    xlab(NULL) +
    ggtitle(paste0(panel_label, ". ", gsub("\\.", " ", species))) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12, face = "bold.italic", hjust = 0),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    )

  return(p)
}

# Create simplified plots
simplified_plots <- list()
for (i in seq_along(staph_species)) {
  species <- staph_species[i]
  panel_label <- panel_labels[i]

  p <- create_staph_plot_simplified(staph_data, species, panel_label)
  if (!is.null(p)) {
    simplified_plots[[i]] <- p
  }
}

# Arrange in grid (3 rows x 2 columns)
n_plots <- length(simplified_plots)
ncols <- 2
nrows <- ceiling(n_plots / ncols)

comprehensive_fig <- arrangeGrob(
  grobs = simplified_plots,
  ncol = ncols,
  nrow = nrows,
  top = textGrob(
    "Supplementary Figure: Staphylococcus Species Comparison by Location, Body Site, and Week\n",
    gp = gpar(fontsize = 18, fontface = "bold")
  ),
  bottom = textGrob(
    "Box plots show median and IQR. Diamonds indicate mean. Log10 scale used for abundance.",
    gp = gpar(fontsize = 11)
  ),
  padding = unit(0.8, "line")
)

# Calculate figure dimensions
fig_width <- ncols * 8
fig_height <- nrows * 6

# Save comprehensive figure
ggsave("revision_figures/Staphylococcus_species_comprehensive.pdf",
       plot = comprehensive_fig,
       width = fig_width,
       height = fig_height,
       limitsize = FALSE)

cat("✓ Saved comprehensive multi-panel figure\n")

# ============================================================================
# Generate summary statistics
# ============================================================================

cat("\n=== Generating summary statistics ===\n")

summary_stats <- data.frame()

for (species in staph_species) {
  if (species %in% colnames(staph_data)) {

    species_data <- staph_data %>%
      select(Sample, Location, SampleType, SampleCollectionWeek, all_of(species)) %>%
      rename(abundance = all_of(species))

    # Overall prevalence by location
    cinci_prev <- sum(species_data$abundance[species_data$Location == "Cincinnati"] > 0) /
                  sum(species_data$Location == "Cincinnati") * 100
    hangz_prev <- sum(species_data$abundance[species_data$Location == "Hangzhou"] > 0) /
                  sum(species_data$Location == "Hangzhou") * 100

    # Mean abundance (non-zero only)
    cinci_mean <- mean(species_data$abundance[species_data$Location == "Cincinnati" &
                                               species_data$abundance > 0], na.rm = TRUE)
    hangz_mean <- mean(species_data$abundance[species_data$Location == "Hangzhou" &
                                               species_data$abundance > 0], na.rm = TRUE)

    # Overall Wilcoxon test
    overall_test <- wilcox.test(abundance ~ Location, data = species_data)

    summary_stats <- rbind(summary_stats, data.frame(
      Species = species,
      Cincinnati_Prevalence = cinci_prev,
      Hangzhou_Prevalence = hangz_prev,
      Cincinnati_Mean_NonZero = ifelse(is.nan(cinci_mean), 0, cinci_mean),
      Hangzhou_Mean_NonZero = ifelse(is.nan(hangz_mean), 0, hangz_mean),
      Overall_P_value = overall_test$p.value,
      Significant = ifelse(overall_test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# Sort by prevalence difference
summary_stats <- summary_stats %>%
  mutate(Prevalence_Diff = abs(Cincinnati_Prevalence - Hangzhou_Prevalence)) %>%
  arrange(desc(Prevalence_Diff))

cat("\nSummary statistics:\n")
print(summary_stats)

write.csv(summary_stats,
          "revision_figures/Staphylococcus_species_summary_stats.csv",
          row.names = FALSE)
cat("✓ Saved summary statistics\n")

# ============================================================================
# Final Summary
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("STAPHYLOCOCCUS SPECIES COMPARISON COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n\n")

cat("Generated files:\n")
cat("  Individual plots (6 PDFs):\n")
for (species in staph_species) {
  clean_name <- gsub("\\.", "_", species)
  cat(sprintf("    - revision_figures/%s_comparison.pdf\n", clean_name))
}
cat("\n  Comprehensive multi-panel figure:\n")
cat("    - revision_figures/Staphylococcus_species_comprehensive.pdf\n")
cat("\n  Statistics tables:\n")
cat("    - revision_figures/Staphylococcus_species_sample_counts.csv\n")
cat("    - revision_figures/Staphylococcus_species_pvalues.csv\n")
cat("    - revision_figures/Staphylococcus_species_summary_stats.csv\n")

cat("\nKey findings:\n")
sig_species <- summary_stats %>% filter(Significant == "Yes")
if (nrow(sig_species) > 0) {
  cat(sprintf("  %d species with significant differences (p<0.05):\n", nrow(sig_species)))
  for (i in 1:nrow(sig_species)) {
    cat(sprintf("    - %s: %.1f%% Cincinnati vs %.1f%% Hangzhou (p=%.4f)\n",
                gsub("\\.", " ", sig_species$Species[i]),
                sig_species$Cincinnati_Prevalence[i],
                sig_species$Hangzhou_Prevalence[i],
                sig_species$Overall_P_value[i]))
  }
} else {
  cat("  No species with overall significant differences\n")
  cat("  (Check individual body site/week p-values in pvalues table)\n")
}

cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
