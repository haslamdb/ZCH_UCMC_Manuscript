#!/usr/bin/env Rscript
# Generate comprehensive Figure 4 style plot for top 6 Candida species
# Uses same format as generate_figure4_comprehensive_revised.R
# but replaces bacterial species with most abundant Candida species

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(ggthemes)  # for scale_color_tableau()

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript/python_scripts")

# Create output directory
dir.create("revision_figures", showWarnings = FALSE)

# Load the saved R workspace from the main analysis
cat("Loading workspace...\n")
load("/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514")
cat("Workspace loaded successfully!\n")

cat("\n=== Generating Comprehensive Candida Figure (Figure 4 Style) ===\n")

# ============================================================================
# Identify top 6 Candida species by mean abundance
# ============================================================================

# Get all Candida species from the dataset
candida_cols <- grep("^Candida\\.", colnames(NICUSpeciesNR), value = TRUE)

cat("\nFound", length(candida_cols), "Candida species in dataset:\n")
print(candida_cols)

# Filter data for relevant samples and rename locations
candida_data <- NICUSpeciesNR %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         SampleCollectionWeek %in% c("Week.1", "Week.3")) %>%
  mutate(Location = recode(Location,
                           "Cincinnati" = "UCMC",
                           "Hangzhou" = "ZCH"))

# Calculate mean abundance for each Candida species
candida_summary <- data.frame(
  species = candida_cols,
  mean_abundance = sapply(candida_cols, function(sp) {
    mean(candida_data[[sp]], na.rm = TRUE)
  }),
  prevalence = sapply(candida_cols, function(sp) {
    sum(candida_data[[sp]] > 0, na.rm = TRUE) / nrow(candida_data) * 100
  }),
  max_abundance = sapply(candida_cols, function(sp) {
    max(candida_data[[sp]], na.rm = TRUE)
  })
) %>%
  arrange(desc(mean_abundance))

cat("\nCandida species summary (sorted by mean abundance):\n")
print(candida_summary)

# Select top 6 by mean abundance
top_candida <- head(candida_summary$species, 6)

cat("\nSelected top 6 Candida species by mean abundance:\n")
for (i in seq_along(top_candida)) {
  cat(sprintf("  %d. %s (mean: %.2f, prevalence: %.1f%%)\n",
              i,
              gsub("\\.", " ", top_candida[i]),
              candida_summary$mean_abundance[i],
              candida_summary$prevalence[i]))
}

# Save summary table
write.csv(candida_summary,
          "revision_figures/Candida_abundance_ranking.csv",
          row.names = FALSE)

# ============================================================================
# Generate Figure 4 style plots for each Candida species
# ============================================================================

# Define color palette (matching Figure 4 - UCMC/ZCH)
col <- c("#E69F00", "#56B4E9")  # UCMC (orange), ZCH (blue)

# Key organisms for Candida Figure 4 (top 6 species a-f)
key_organisms <- top_candida
panel_labels <- letters[1:6]

# Prepare data
DummySpeciesNR <- candida_data

# Function to create individual organism plot with sample counts and p-values
# (identical to Figure 4 format)
create_organism_plot_original_style <- function(data, organism, panel_label) {

  # Get data for this organism
  plot_data <- data %>%
    select(Sample, Location, SampleType, SampleCollectionWeek, all_of(organism)) %>%
    rename(abundance = all_of(organism)) %>%
    filter(!is.na(abundance))

  # Calculate sample sizes for each group
  sample_counts <- plot_data %>%
    group_by(Location, SampleCollectionWeek) %>%
    summarise(
      n = n(),
      max_abundance = max(abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(y_pos = max_abundance * 3)  # Position above boxes

  # Calculate p-values between locations for each week
  pvalue_results <- data.frame()

  for (week in c("Week.1", "Week.3")) {
    subset_data <- plot_data %>%
      filter(SampleCollectionWeek == week)

    if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
      test_result <- tryCatch({
        wilcox.test(abundance ~ Location, data = subset_data)
      }, error = function(e) {
        list(p.value = NA)
      })

      # Get max y value for p-value annotation position
      max_y <- max(subset_data$abundance, na.rm = TRUE)

      pvalue_results <- rbind(pvalue_results, data.frame(
        SampleCollectionWeek = week,
        p_value = test_result$p.value,
        y_pos = max_y * 5,  # Position p-value above sample counts
        p_label = ifelse(is.na(test_result$p.value), "",
                  ifelse(test_result$p.value < 0.001, "p<0.001***",
                  ifelse(test_result$p.value < 0.01, sprintf("p=%.3f**", test_result$p.value),
                  ifelse(test_result$p.value < 0.05, sprintf("p=%.3f*", test_result$p.value),
                         sprintf("p=%.2f", test_result$p.value))))),
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("\nP-values for", organism, ":\n")
  print(pvalue_results)

  # Clean up week labels (remove dots)
  plot_data <- plot_data %>%
    mutate(SampleCollectionWeek = gsub("\\.", " ", SampleCollectionWeek))

  sample_counts <- sample_counts %>%
    mutate(SampleCollectionWeek = gsub("\\.", " ", SampleCollectionWeek))

  pvalue_results <- pvalue_results %>%
    mutate(SampleCollectionWeek = gsub("\\.", " ", SampleCollectionWeek))

  # Create plot matching original style
  p <- ggplot(plot_data,
              aes(x = Location, y = as.numeric(abundance))) +
    geom_boxplot(lwd = 1, aes(color = factor(Location)),
                 fill = NA,  # No gray fill in boxes
                 outlier.size = 1.8) +  # Reduced to 60% of 3
    stat_summary(fun = mean, geom = "point", shape = 5, size = 4.8) +  # Reduced to 60% of 8
    scale_colour_manual(values = col) +
    geom_point(size = 2.4, aes(color = factor(Location))) +  # Reduced to 60% of 4
    scale_y_log10() +
    xlab(NULL) +
    ylab(gsub("\\.", " ", organism)) +  # Removed \n to reduce spacing
    facet_grid(rows = vars(SampleCollectionWeek)) +
    scale_color_tableau() +
    theme_bw() +
    theme(
      axis.title.y = element_text(face = "italic", size = 18),  # Increased 150% from 12
      axis.text.x = element_text(size = 14.3),  # Increased 130% from 11
      axis.text.y = element_text(size = 13.2),  # Increased 120% from 11
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),  # Remove gray background from strip
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    # Add sample count labels
    geom_text(data = sample_counts,
              aes(x = Location, y = y_pos, label = paste0("n=", n)),
              size = 4,
              fontface = "bold",
              inherit.aes = FALSE) +
    # Add p-value labels (back to original size)
    geom_text(data = pvalue_results,
              aes(x = 1.5, y = y_pos, label = p_label),
              size = 4,  # Back to original size
              fontface = "bold",
              inherit.aes = FALSE)

  # Add panel label
  p <- p + ggtitle(paste0(panel_label, "."))

  return(p)
}

# Generate plots for all 6 Candida organisms
cat("\nGenerating plots for all 6 Candida species...\n")

plot_list <- list()

for (i in seq_along(key_organisms)) {
  organism <- key_organisms[i]
  panel_label <- panel_labels[i]

  cat(sprintf("\nProcessing %s. %s\n", panel_label, organism))

  p <- create_organism_plot_original_style(DummySpeciesNR, organism, panel_label)
  plot_list[[i]] <- p
}

# Create comprehensive 3-column grid layout (2 rows, 3 per row)
cat("\nCreating comprehensive Candida figure...\n")

comprehensive_fig <- arrangeGrob(
  grobs = plot_list,
  ncol = 3,
  nrow = 2,
  top = textGrob("Supplementary Figure. Longitudinal changes in Candida species abundance by location\n",
                 gp = gpar(fontsize = 16, fontface = "bold")),
  padding = unit(0.8, "line"),  # Increased padding between plots (was 0.5)
  widths = unit(rep(3.78, 3), "in"),  # 90% of 4.2" = 3.78"
  heights = unit(rep(4.5, 2), "in")   # Same height
)

# Save comprehensive figure (3 columns, 2 rows with more spacing)
ggsave("revision_figures/SuppFig_Candida_Figure4_style_comprehensive.pdf",
       plot = comprehensive_fig,
       width = 13.5,  # 3 columns * 4.5" per column
       height = 10,
       limitsize = FALSE)

cat("\n✓ Saved comprehensive Candida figure in Figure 4 style\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("COMPREHENSIVE CANDIDA FIGURE GENERATION COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("\nGenerated files:\n")
cat("  - revision_figures/SuppFig_Candida_Figure4_style_comprehensive.pdf\n")
cat("  - revision_figures/Candida_abundance_ranking.csv\n")
cat("\nFormat:\n")
cat("  - 3-column layout (2 rows, 3 plots per row)\n")
cat("  - 6 species total (a-f) - top 6 by mean abundance\n")
cat("  - Original box plot style with Location on x-axis\n")
cat("  - Faceted by SampleCollectionWeek (Week 1, Week 3 - no dots)\n")
cat("  - Boxplots with no gray fill (transparent)\n")
cat("  - Sample counts (n=) added above each box\n")
cat("  - P-values shown between locations (original font size)\n")
cat("  - Y-axis species names 150% larger and italicized\n")
cat("  - Gray background removed from strip labels\n")
cat("  - Original color scheme (UCMC orange, ZCH blue)\n")
cat("\nTop 6 Candida species by mean abundance:\n")
for (i in seq_along(top_candida)) {
  cat(sprintf("  %s. %s\n", panel_labels[i], gsub("\\.", " ", top_candida[i])))
}
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
