#!/usr/bin/env Rscript
# Visualize antibiotic-pathogen correlations by body site
# Creates heatmaps and scatter plots stratified by site

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(ggthemes)
library(tibble)

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript")

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("ANTIBIOTIC-PATHOGEN CORRELATION VISUALIZATIONS (BY BODY SITE)\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n\n")

# ============================================================================
# Load data
# ============================================================================

cat("Loading correlation results and merged data...\n")

correlations <- read.csv("antibiotic_pathogen_analysis/antibiotic_pathogen_correlations_by_site.csv")
merged_data <- read.csv("antibiotic_pathogen_analysis/antibiotic_pathogen_merged_by_site.csv")

cat(sprintf("✓ Loaded %d correlations\n", nrow(correlations)))
cat(sprintf("✓ Loaded merged data with %d samples\n", nrow(merged_data)))

# ============================================================================
# Create correlation heatmaps (one for each site × week combination)
# ============================================================================

cat("\nCreating correlation heatmaps by body site...\n")

# Function to create heatmap for one site and week
create_heatmap <- function(cor_data, site, week_label, title_text) {

  # Clean up names for display
  cor_data <- cor_data %>%
    mutate(
      Antibiotic_clean = gsub("_w1$|_cumulative_w3$", "", Antibiotic),
      Antibiotic_clean = gsub("_", " ", Antibiotic_clean),
      Pathogen_clean = gsub("\\.", " ", Pathogen)
    )

  # Create significance labels
  plot_data <- cor_data %>%
    mutate(
      sig_label = case_when(
        P_adjusted < 0.001 ~ "***",
        P_adjusted < 0.01 ~ "**",
        P_adjusted < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  # Create heatmap
  p <- ggplot(plot_data, aes(x = Pathogen_clean, y = Antibiotic_clean, fill = Rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 5, fontface = "bold", vjust = 0.75) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90",
      name = "Spearman rho"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    ) +
    ggtitle(title_text) +
    coord_fixed()

  return(p)
}

# Create heatmaps for each site × week combination
heatmap_plots <- list()
plot_idx <- 1

for (site in c("Axilla", "Groin", "Stool")) {
  # Week 1
  w1_data <- correlations %>% filter(SampleType == site, Week == "Week1")
  if (nrow(w1_data) > 0) {
    p_w1 <- create_heatmap(w1_data, site, "Week 1",
                           sprintf("%s - Week 1", site))
    heatmap_plots[[plot_idx]] <- p_w1
    plot_idx <- plot_idx + 1
  }

  # Week 3
  w3_data <- correlations %>% filter(SampleType == site, Week == "Week3_Cumulative")
  if (nrow(w3_data) > 0) {
    p_w3 <- create_heatmap(w3_data, site, "Week 3",
                           sprintf("%s - Week 3 Cumulative", site))
    heatmap_plots[[plot_idx]] <- p_w3
    plot_idx <- plot_idx + 1
  }
}

# Save combined heatmap grid (3 sites × 2 weeks = 6 panels)
combined_heatmap <- arrangeGrob(
  grobs = heatmap_plots,
  ncol = 2,
  nrow = 3,
  top = textGrob("Antibiotic-Pathogen Correlations by Body Site\n(* p<0.05, ** p<0.01, *** p<0.001, FDR-adjusted)",
                 gp = gpar(fontsize = 16, fontface = "bold")),
  left = textGrob("Body Site", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)

ggsave("antibiotic_pathogen_analysis/Correlation_Heatmap_By_Site.pdf",
       plot = combined_heatmap,
       width = 20,
       height = 24,
       limitsize = FALSE)

cat("✓ Saved combined heatmap: antibiotic_pathogen_analysis/Correlation_Heatmap_By_Site.pdf\n")

# ============================================================================
# Create scatter plots for significant associations
# ============================================================================

cat("\nCreating scatter plots for significant associations...\n")

significant_cors <- correlations %>%
  filter(P_adjusted < 0.05) %>%
  arrange(P_adjusted)

if (nrow(significant_cors) == 0) {
  cat("No significant correlations to plot.\n")
} else {
  cat(sprintf("Found %d significant correlations to plot\n", nrow(significant_cors)))

  # Color palette
  col <- c("UCMC" = "#E69F00", "ZCH" = "#56B4E9")

  plot_list <- list()

  for (i in 1:nrow(significant_cors)) {
    cor_row <- significant_cors[i, ]

    # Get column names
    abx_col <- cor_row$Antibiotic
    pathogen_col <- paste0(cor_row$Pathogen, "_",
                          ifelse(cor_row$Week == "Week1", "w1", "w3"))
    site <- cor_row$SampleType

    # Clean names for labels
    abx_label <- gsub("_w1$|_cumulative_w3$", "", abx_col)
    abx_label <- gsub("_", " ", abx_label)
    if (cor_row$Week == "Week3_Cumulative") {
      abx_label <- paste0(abx_label, " (cumulative wk3)")
    } else {
      abx_label <- paste0(abx_label, " (wk1)")
    }

    pathogen_label <- gsub("\\.", " ", cor_row$Pathogen)
    week_label <- ifelse(cor_row$Week == "Week1", "wk1", "wk3")

    # Extract data for this comparison (filter by site)
    plot_data <- merged_data %>%
      filter(SampleType == site) %>%
      select(Subject, Location, all_of(abx_col), all_of(pathogen_col)) %>%
      rename(
        antibiotic = all_of(abx_col),
        pathogen = all_of(pathogen_col)
      ) %>%
      filter(!is.na(pathogen))

    # Create scatter plot
    p <- ggplot(plot_data, aes(x = antibiotic, y = pathogen, color = Location)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, aes(fill = Location), alpha = 0.2) +
      scale_color_manual(values = col) +
      scale_fill_manual(values = col) +
      scale_y_log10() +
      xlab(paste("Days of", abx_label)) +
      ylab(paste(pathogen_label, "abundance (", week_label, ")", sep = "")) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
      ) +
      ggtitle(sprintf("%s samples", site)) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("ρ = %.3f\np = %.2e\nFDR = %.2e\nn = %d",
                             cor_row$Rho, cor_row$P_value,
                             cor_row$P_adjusted, cor_row$N),
               hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold")

    plot_list[[i]] <- p
  }

  # Save individual plots
  for (i in 1:length(plot_list)) {
    filename <- sprintf("antibiotic_pathogen_analysis/Scatter_Significant_Site_%d.pdf", i)
    ggsave(filename, plot = plot_list[[i]], width = 8, height = 6)
    cat(sprintf("✓ Saved scatter plot %d: %s\n", i, filename))
  }

  # Save combined scatter plot
  if (length(plot_list) > 0) {
    ncols <- min(2, length(plot_list))
    nrows <- ceiling(length(plot_list) / ncols)

    combined_scatter <- arrangeGrob(
      grobs = plot_list,
      ncol = ncols,
      nrow = nrows,
      top = textGrob("Significant Antibiotic-Pathogen Associations by Body Site (FDR < 0.05)",
                     gp = gpar(fontsize = 16, fontface = "bold"))
    )

    ggsave("antibiotic_pathogen_analysis/Scatter_Significant_By_Site_Combined.pdf",
           plot = combined_scatter,
           width = 8 * ncols,
           height = 6 * nrows,
           limitsize = FALSE)

    cat("✓ Saved combined scatter plot: antibiotic_pathogen_analysis/Scatter_Significant_By_Site_Combined.pdf\n")
  }
}

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("VISUALIZATION COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
cat("\nGenerated files:\n")
cat("  Heatmaps:\n")
cat("    - antibiotic_pathogen_analysis/Correlation_Heatmap_By_Site.pdf\n")
if (nrow(significant_cors) > 0) {
  cat("  Scatter plots:\n")
  for (i in 1:nrow(significant_cors)) {
    cat(sprintf("    - antibiotic_pathogen_analysis/Scatter_Significant_Site_%d.pdf\n", i))
  }
  cat("    - antibiotic_pathogen_analysis/Scatter_Significant_By_Site_Combined.pdf\n")
}
cat("\n")
cat(paste(rep("=", 80), collapse = ""))
cat("\n")
