#!/usr/bin/env Rscript
# Generate comprehensive Figure 4 with faceting by SampleType (rows) and Week (columns)
# This matches the original StaphLocation style more closely

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

cat("\n=== Generating Comprehensive Figure 4 (Faceted by Site and Week) ===\n")

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

# Prepare data - filter for 3 sample types
DummySpeciesNR <- NICUSpeciesNR %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         SampleCollectionWeek %in% c("Week.1", "Week.3"))

# Function to create organism plot with faceting by SampleType (rows) and Week (columns)
create_organism_plot_faceted <- function(data, organism, panel_label) {

  # Get data for this organism
  plot_data <- data %>%
    select(Sample, Location, SampleType, SampleCollectionWeek, all_of(organism)) %>%
    rename(abundance = all_of(organism)) %>%
    filter(!is.na(abundance))

  # Clean up week labels (remove dots)
  plot_data <- plot_data %>%
    mutate(SampleCollectionWeek = gsub("\\.", " ", SampleCollectionWeek))

  # Calculate sample sizes for each group
  sample_counts <- plot_data %>%
    group_by(Location, SampleType, SampleCollectionWeek) %>%
    summarise(
      n = n(),
      max_abundance = max(abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(y_pos = max_abundance * 3)  # Position above boxes

  # Calculate p-values between locations for each SampleType x Week combination
  pvalue_results <- data.frame()

  for (week in c("Week 1", "Week 3")) {
    for (site in c("Axilla", "Groin", "Stool")) {
      subset_data <- plot_data %>%
        filter(SampleCollectionWeek == week, SampleType == site)

      if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
        test_result <- tryCatch({
          wilcox.test(abundance ~ Location, data = subset_data)
        }, error = function(e) {
          list(p.value = NA)
        })

        pvalue_results <- rbind(pvalue_results, data.frame(
          SampleType = site,
          SampleCollectionWeek = week,
          p_value = test_result$p.value,
          p_label = ifelse(is.na(test_result$p.value), "",
                    ifelse(test_result$p.value < 0.001, "p<0.001***",
                    ifelse(test_result$p.value < 0.01, sprintf("p=%.3f**", test_result$p.value),
                    ifelse(test_result$p.value < 0.05, sprintf("p=%.3f*", test_result$p.value),
                           sprintf("p=%.2f", test_result$p.value))))),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  cat("\nP-values for", organism, ":\n")
  print(pvalue_results)

  # Calculate y-axis limits with padding at top for p-values and sample counts
  y_min <- min(plot_data$abundance[plot_data$abundance > 0], na.rm = TRUE)
  y_max <- max(plot_data$abundance, na.rm = TRUE)

  # Add more padding at top (on log scale) to prevent text overlap
  y_limits <- c(y_min * 0.8, y_max * 100)  # 100x gives more room at top on log scale

  # Create plot with faceting by SampleType (rows) and SampleCollectionWeek (columns)
  p <- ggplot(plot_data,
              aes(x = Location, y = as.numeric(abundance))) +
    geom_boxplot(lwd = 0.75, aes(color = factor(Location)),  # Reduced to 75% of 1
                 fill = NA,  # No gray fill in boxes
                 outlier.size = 1.35) +  # 75% of 1.8 (which was 60% of 3)
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3.6) +  # 75% of 4.8
    scale_colour_manual(values = col) +
    geom_point(size = 1.8, aes(color = factor(Location))) +  # 75% of 2.4
    scale_y_log10(limits = y_limits) +  # Add limits with padding
    xlab(NULL) +
    ylab(gsub("\\.", " ", organism)) +  # Remove dots, tight spacing
    facet_grid(SampleType ~ SampleCollectionWeek) +  # Rows = SampleType, Columns = Week
    scale_color_tableau() +
    theme_bw() +
    theme(
      axis.title.y = element_text(face = "italic", size = 18),
      axis.text.x = element_text(size = 10.725),  # 75% of 14.3
      axis.text.y = element_text(size = 13.2),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
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
    # Add p-value labels at top of plot (using Inf for y position)
    geom_text(data = pvalue_results,
              aes(x = 1.5, y = Inf, label = p_label),
              size = 4,
              fontface = "bold",
              vjust = 1.5,  # Position just below the top
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

  p <- create_organism_plot_faceted(DummySpeciesNR, organism, panel_label)
  plot_list[[i]] <- p
}

# Create comprehensive 4-column grid layout (2 rows, 4 per row)
cat("\nCreating comprehensive figure...\n")

comprehensive_fig <- arrangeGrob(
  grobs = plot_list,
  ncol = 4,
  nrow = 2,
  top = textGrob("Figure 4. Longitudinal changes in BSI-associated species abundance by location and body site\n",
                 gp = gpar(fontsize = 16, fontface = "bold")),
  padding = unit(0.8, "line"),
  widths = unit(rep(3.78, 4), "in"),
  heights = unit(rep(5.4, 2), "in")  # Increased by 20% from 4.5 to 5.4
)

# Save comprehensive figure with different name
ggsave("revision_figures/Figure4_comprehensive_faceted_site_and_week.pdf",
       plot = comprehensive_fig,
       width = 18,
       height = 12,  # Increased by 20% from 10 to 12
       limitsize = FALSE)

cat("\n✓ Saved comprehensive Figure 4 with site and week faceting\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("COMPREHENSIVE FIGURE 4 (SITE+WEEK FACETING) GENERATION COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("\nGenerated file:\n")
cat("  - revision_figures/Figure4_comprehensive_faceted_site_and_week.pdf\n")
cat("\nFormat:\n")
cat("  - 4-column layout (2 rows, 4 plots per row)\n")
cat("  - 8 species total (a-h)\n")
cat("  - Original box plot style with Location on x-axis\n")
cat("  - Faceted by SampleType (rows: Axilla, Groin, Stool)\n")
cat("  - Faceted by SampleCollectionWeek (columns: Week 1, Week 3)\n")
cat("  - Boxplots with no gray fill (transparent)\n")
cat("  - Sample counts (n=) added above each box\n")
cat("  - P-values shown between locations\n")
cat("  - Point sizes reduced to 60% of original\n")
cat("  - Original color scheme (Cincinnati orange, Hangzhou blue)\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
