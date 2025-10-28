#!/usr/bin/env Rscript
# Generate comprehensive Figure 4 with original box plot style
# Matches the format from MicrobiomeAnalysisNICU script

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

cat("\n=== Generating Comprehensive Figure 4 (Original Style) ===\n")

# Define color palette (matching your original)
col <- c("#E69F00", "#56B4E9")  # Cincinnati (orange), Hangzhou (blue)

# Key organisms for Figure 4 (all 8 species a-h)
key_organisms <- c("Staphylococcus.aureus",        # a
                   "Klebsiella.pneumoniae",        # b
                   "Klebsiella.oxytoca",           # c
                   "Enterococcus.faecium",         # d
                   "Enterococcus.faecalis",        # e
                   "Serratia.marcescens",          # f
                   "Escherichia.coli",             # g
                   "Streptococcus.pyogenes")       # h

panel_labels <- letters[1:8]

# Prepare data
DummySpeciesNR <- NICUSpeciesNR %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         SampleCollectionWeek %in% c("Week.1", "Week.3"))

# Function to create individual organism plot with sample counts and p-values
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

  # Create plot matching original style
  p <- ggplot(plot_data,
              aes(x = Location, y = as.numeric(abundance), fill = Location)) +
    geom_boxplot(lwd = 1, aes(color = factor(Location), fill = NA),
                 outlier.size = 3) +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
    scale_fill_manual(values = col) +
    scale_colour_manual(values = col) +
    geom_point(size = 4, aes(color = factor(Location))) +
    scale_y_log10() +
    xlab(NULL) +
    ylab(paste0(gsub("\\.", " ", organism), "\n")) +
    facet_grid(rows = vars(SampleCollectionWeek)) +
    scale_color_tableau() +
    theme_bw() +
    theme(
      axis.title.y = element_text(face = "italic", size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      strip.text = element_text(size = 12, face = "bold"),
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
    # Add p-value labels
    geom_text(data = pvalue_results,
              aes(x = 1.5, y = y_pos, label = p_label),
              size = 4,
              fontface = "bold",
              inherit.aes = FALSE)

  # Add panel label
  p <- p + ggtitle(paste0(panel_label, "."))

  return(p)
}

# Generate plots for all 8 organisms
cat("\nGenerating plots for all 8 organisms...\n")

plot_list <- list()

for (i in seq_along(key_organisms)) {
  organism <- key_organisms[i]
  panel_label <- panel_labels[i]

  cat(sprintf("\nProcessing %s. %s\n", panel_label, organism))

  p <- create_organism_plot_original_style(DummySpeciesNR, organism, panel_label)
  plot_list[[i]] <- p
}

# Create comprehensive 4x2 grid layout
cat("\nCreating comprehensive figure...\n")

comprehensive_fig <- arrangeGrob(
  grobs = plot_list,
  ncol = 2,
  nrow = 4,
  top = textGrob("Figure 4. Longitudinal changes in BSI-associated species abundance by location\n",
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save comprehensive figure
ggsave("revision_figures/Figure4_comprehensive_original_style.pdf",
       plot = comprehensive_fig,
       width = 14,
       height = 20,
       limitsize = FALSE)

cat("\n✓ Saved comprehensive Figure 4 with original style\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("COMPREHENSIVE FIGURE 4 GENERATION COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("\nGenerated file:\n")
cat("  - revision_figures/Figure4_comprehensive_original_style.pdf\n")
cat("\nFormat:\n")
cat("  - 4x2 grid (8 species: a-h)\n")
cat("  - Original box plot style with Location on x-axis\n")
cat("  - Faceted by SampleCollectionWeek (Week.1, Week.3)\n")
cat("  - Sample counts (n=) added above each box\n")
cat("  - P-values shown between locations for each week\n")
cat("  - Original color scheme (Cincinnati orange, Hangzhou blue)\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
