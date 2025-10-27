#!/usr/bin/env Rscript
# Generate Supplementary Figure: Longitudinal Candida Species Abundance
# Addressing Reviewer 2's question about 20 Candidemia cases at ZCH
# Similar format to Figure 4 (longitudinal changes in BSI-associated species)

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript/python_scripts")

# Create output directory
dir.create("revision_figures", showWarnings = FALSE)

# Load the saved R workspace from the main analysis
cat("Loading workspace...\n")
load("/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514")
cat("Workspace loaded successfully!\n")

cat("\n=== Generating Candida Supplementary Figure ===\n")

# ============================================================================
# First, identify top Candida species by prevalence and abundance
# ============================================================================

# Get all Candida species from the dataset
candida_cols <- grep("^Candida\\.", colnames(NICUSpeciesNR), value = TRUE)

cat("\nFound", length(candida_cols), "Candida species in dataset:\n")
print(candida_cols)

# Filter data for relevant samples
candida_data <- NICUSpeciesNR %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         SampleCollectionWeek %in% c("Week.1", "Week.3"))

# Calculate prevalence and mean abundance for each Candida species
candida_summary <- data.frame(
  species = candida_cols,
  prevalence = sapply(candida_cols, function(sp) {
    sum(candida_data[[sp]] > 0, na.rm = TRUE) / nrow(candida_data) * 100
  }),
  mean_abundance = sapply(candida_cols, function(sp) {
    mean(candida_data[[sp]], na.rm = TRUE)
  }),
  max_abundance = sapply(candida_cols, function(sp) {
    max(candida_data[[sp]], na.rm = TRUE)
  })
) %>%
  arrange(desc(prevalence))

cat("\nCandida species summary (sorted by prevalence):\n")
print(candida_summary)

# Save summary table
write.csv(candida_summary,
          "revision_figures/Candida_species_summary.csv",
          row.names = FALSE)

# Select top 4-6 species (prioritize clinically relevant + prevalent)
# C. albicans, C. parapsilosis, C. tropicalis are most clinically relevant
# Pick top species by prevalence or clinical relevance
top_candida <- candida_summary %>%
  filter(prevalence > 5 | species %in% c("Candida.albicans",
                                          "Candida.parapsilosis",
                                          "Candida.tropicalis",
                                          "Candida.auris")) %>%
  head(6) %>%
  pull(species)

cat("\nSelected Candida species for figure:\n")
print(top_candida)

# If fewer than 4 species meet criteria, just take top 6 by prevalence
if (length(top_candida) < 4) {
  top_candida <- head(candida_summary$species, 6)
  cat("\nUsing top 6 by prevalence instead:\n")
  print(top_candida)
}

# ============================================================================
# Create longitudinal plots for each Candida species (similar to Figure 4)
# ============================================================================

# Function to create organism-specific plot with sample counts
create_candida_plot <- function(data, organism, title, simplified = FALSE) {

  # Select relevant data
  plot_data <- data %>%
    select(Sample, Location, SampleType, SampleCollectionWeek, all_of(organism)) %>%
    rename(abundance = all_of(organism)) %>%
    filter(!is.na(abundance)) %>%
    mutate(abundance_log = ifelse(abundance == 0, 1, abundance))  # Replace 0 with 1 for log scale

  # Calculate sample sizes
  sample_counts <- plot_data %>%
    group_by(Location, SampleType, SampleCollectionWeek) %>%
    summarise(
      n = n(),
      median_abundance = median(abundance_log, na.rm = TRUE),
      y_pos = max(abundance_log, na.rm = TRUE) * 1.5,
      .groups = 'drop'
    )

  cat("\nSample counts for", title, ":\n")
  print(sample_counts)

  # Create plot
  if (simplified) {
    # Simplified version for multi-panel figure
    p <- ggplot(plot_data,
                aes(x = SampleCollectionWeek, y = abundance_log, fill = SampleCollectionWeek)) +
      geom_boxplot(lwd = 0.8, aes(color = factor(SampleCollectionWeek), fill = NA),
                   outlier.size = 1.5) +
      stat_summary(fun = mean, geom = "point", shape = 5, size = 4) +
      geom_point(size = 1.5, aes(color = factor(SampleCollectionWeek)), alpha = 0.4) +
      facet_grid(Location ~ SampleType) +
      scale_y_log10() +
      scale_color_manual(values = c("Week.1" = "#00BFC4", "Week.3" = "#F8766D")) +
      ylab("Abundance (log)") +
      xlab(NULL) +
      ggtitle(title) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
  } else {
    # Full version for individual plots
    p <- ggplot(plot_data,
                aes(x = SampleCollectionWeek, y = abundance_log, fill = SampleCollectionWeek)) +
      geom_boxplot(lwd = 1, aes(color = factor(SampleCollectionWeek), fill = NA),
                   outlier.size = 2) +
      stat_summary(fun = mean, geom = "point", shape = 5, size = 6) +
      geom_point(size = 2, aes(color = factor(SampleCollectionWeek)), alpha = 0.5) +
      # Add sample count labels
      geom_text(data = sample_counts,
                aes(x = SampleCollectionWeek, y = y_pos, label = paste0("n=", n)),
                size = 4,
                fontface = "bold",
                inherit.aes = FALSE) +
      facet_grid(Location ~ SampleType) +
      scale_y_log10() +
      scale_color_manual(values = c("Week.1" = "#00BFC4", "Week.3" = "#F8766D")) +
      ylab("Abundance (log scale)\n") +
      xlab(NULL) +
      ggtitle(paste("Longitudinal changes in", title)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
  }

  return(p)
}

# Generate individual plots for each Candida species
cat("\nGenerating individual plots...\n")
for (organism in top_candida) {
  # Clean name for title
  clean_name <- gsub("\\.", " ", organism)

  # Create plot
  p <- create_candida_plot(candida_data, organism, clean_name, simplified = FALSE)

  # Save individual plot
  filename <- paste0("revision_figures/SuppFig_Candida_", gsub("\\.", "_", organism), ".pdf")
  ggsave(filename, plot = p, width = 14, height = 8)
  cat("  ✓ Saved", filename, "\n")
}

# ============================================================================
# Create comprehensive multi-panel figure with all Candida species
# ============================================================================

cat("\nGenerating comprehensive multi-panel Candida figure...\n")

# Create simplified plots for combined figure
simplified_plots <- lapply(top_candida, function(organism) {
  clean_name <- gsub("\\.", " ", organism)
  create_candida_plot(candida_data, organism, clean_name, simplified = TRUE)
})

# Determine layout based on number of species
n_species <- length(top_candida)
if (n_species <= 4) {
  ncol_layout <- 2
  nrow_layout <- 2
} else if (n_species <= 6) {
  ncol_layout <- 2
  nrow_layout <- 3
} else {
  ncol_layout <- 3
  nrow_layout <- ceiling(n_species / 3)
}

# Combine into multi-panel figure
combined_fig <- arrangeGrob(
  grobs = simplified_plots,
  ncol = ncol_layout,
  nrow = nrow_layout,
  top = textGrob("Supplementary Figure: Longitudinal changes in Candida species abundance\n",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Save comprehensive figure
ggsave("revision_figures/SuppFig_Candida_comprehensive.pdf",
       plot = combined_fig,
       width = 14,
       height = 5 * nrow_layout,
       limitsize = FALSE)

cat("✓ Saved comprehensive Candida figure\n")

# ============================================================================
# Calculate statistics for Candida species by location
# ============================================================================

cat("\n=== Calculating Candida statistics by location ===\n")

candida_stats <- data.frame()

for (organism in top_candida) {
  for (week in c("Week.1", "Week.3")) {
    subset_data <- candida_data %>%
      filter(SampleCollectionWeek == week) %>%
      select(Sample, Location, all_of(organism)) %>%
      rename(abundance = all_of(organism))

    if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
      # Wilcoxon test
      test_result <- tryCatch({
        wilcox.test(abundance ~ Location, data = subset_data)
      }, error = function(e) {
        list(p.value = NA)
      })

      cinci_data <- subset_data %>% filter(Location == "Cincinnati")
      hangz_data <- subset_data %>% filter(Location == "Hangzhou")

      candida_stats <- rbind(candida_stats, data.frame(
        Species = gsub("\\.", " ", organism),
        Week = week,
        n_Cincinnati = nrow(cinci_data),
        median_Cincinnati = median(cinci_data$abundance, na.rm = TRUE),
        prevalence_Cincinnati = sum(cinci_data$abundance > 0) / nrow(cinci_data) * 100,
        n_Hangzhou = nrow(hangz_data),
        median_Hangzhou = median(hangz_data$abundance, na.rm = TRUE),
        prevalence_Hangzhou = sum(hangz_data$abundance > 0) / nrow(hangz_data) * 100,
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\nCandida statistics by location:\n")
print(candida_stats)

# Save statistics
write.csv(candida_stats,
          "revision_figures/Candida_statistics_by_location.csv",
          row.names = FALSE)

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("CANDIDA SUPPLEMENTARY FIGURE GENERATION COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("\nGenerated files in revision_figures/:\n")
cat("  - SuppFig_Candida_comprehensive.pdf (multi-panel,", n_species, "species)\n")
cat("  -", n_species, "individual Candida species PDFs\n")
cat("  - Candida_species_summary.csv\n")
cat("  - Candida_statistics_by_location.csv\n")
cat("\nCandida species included:\n")
for (i in seq_along(top_candida)) {
  cat(sprintf("  %d. %s\n", i, gsub("\\.", " ", top_candida[i])))
}
cat("\nThis addresses Reviewer 2's question about Candida in the microbiome\n")
cat("(20 Candidemia cases at ZCH)\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
