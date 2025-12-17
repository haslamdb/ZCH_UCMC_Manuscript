#!/usr/bin/env Rscript
# Shannon Diversity Figure - Location comparison with p-values
# Faceted by Collection Week (columns) and Sample Type (rows)
# X-axis: UCMC vs ZCH

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(ggthemes)

# Set paths
PROJECT_DIR <- "/home/david/projects/ZCH_UCMC_Manuscript"
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

# Create results directory if needed
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Color palette (UCMC orange, ZCH blue)
col <- c("#E69F00", "#56B4E9")

# Read diversity metrics from master metadata
DATA_FILE <- file.path(PROJECT_DIR, "nicu_resistome_analysis/results/qc/master_metadata.tsv")
cat("Reading data from:", DATA_FILE, "\n")

master_data <- read.delim(DATA_FILE, header = TRUE, stringsAsFactors = FALSE)

# Check if we need to calculate Shannon diversity or if it exists
# The master metadata may not have Shannon - let's check and merge if needed
diversity_file <- file.path(PROJECT_DIR, "nicu_resistome_analysis/results/exploratory/diversity_metrics.tsv")

if (file.exists(diversity_file)) {
  cat("Reading diversity metrics from:", diversity_file, "\n")
  diversity_data <- read.delim(diversity_file, header = TRUE, stringsAsFactors = FALSE)

  # The diversity file has duplicated rows - keep only unique sample_name + Location combinations
  # that match the master metadata
  diversity_data <- diversity_data %>%
    select(sample_name, shannon, richness, simpson) %>%
    distinct()

  # Merge with master metadata
  Diversity <- master_data %>%
    left_join(diversity_data, by = c("sample_name" = "sample_name"))
} else {
  stop("Diversity metrics file not found. Please run diversity calculation first.")
}

# Filter for main sample types and valid locations
Diversity_filtered <- Diversity %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         Location %in% c("UCMC", "ZCH"),
         !is.na(shannon))

# Set factor levels for proper ordering
Diversity_filtered$Location <- factor(Diversity_filtered$Location, levels = c("UCMC", "ZCH"))
Diversity_filtered$SampleType <- factor(Diversity_filtered$SampleType, levels = c("Axilla", "Groin", "Stool"))
Diversity_filtered$SampleCollectionWeek <- factor(Diversity_filtered$SampleCollectionWeek,
                                                   levels = c("Week.1", "Week.3"))

# Print sample sizes
cat("\nSample sizes by group:\n")
print(table(Diversity_filtered$Location, Diversity_filtered$SampleType, Diversity_filtered$SampleCollectionWeek))

# Create individual plots for each facet combination
plot_list <- list()
plot_idx <- 1

for (week in c("Week.1", "Week.3")) {
  for (stype in c("Axilla", "Groin", "Stool")) {

    # Subset data for this facet
    plot_data <- Diversity_filtered %>%
      filter(SampleCollectionWeek == week, SampleType == stype)

    if (nrow(plot_data) == 0) {
      next
    }

    # Calculate sample sizes
    sample_counts <- plot_data %>%
      group_by(Location) %>%
      summarize(n = n(), .groups = "drop") %>%
      mutate(y_pos = min(plot_data$shannon, na.rm = TRUE) - 0.3)

    # Calculate p-value
    ucmc_data <- plot_data$shannon[plot_data$Location == "UCMC"]
    zch_data <- plot_data$shannon[plot_data$Location == "ZCH"]

    if (length(ucmc_data) > 0 && length(zch_data) > 0) {
      test_result <- wilcox.test(ucmc_data, zch_data)
      p_value <- test_result$p.value

      p_label <- ifelse(p_value < 0.001, "p<0.001***",
                ifelse(p_value < 0.01, sprintf("p=%.3f**", p_value),
                ifelse(p_value < 0.05, sprintf("p=%.3f*", p_value),
                       sprintf("p=%.2f", p_value))))
    } else {
      p_label <- ""
    }

    # Create p-value data frame
    pvalue_df <- data.frame(
      x = 1.5,
      y_pos = max(plot_data$shannon, na.rm = TRUE) + 0.3,
      p_label = p_label
    )

    # Create panel title
    panel_title <- paste0(stype, " - ", gsub("\\.", " ", week))

    # Create the plot
    p <- ggplot(plot_data, aes(x = Location, y = shannon)) +
      # Boxplot styling
      geom_boxplot(lwd = 1, aes(color = factor(Location)),
                   fill = NA,
                   outlier.shape = NA) +

      # Mean as diamond
      stat_summary(fun = mean, geom = "point", shape = 5, size = 4.8) +

      # Individual data points
      geom_point(size = 2.4, aes(color = factor(Location)),
                 position = position_jitter(width = 0.15, seed = 42)) +

      # Color scheme
      scale_colour_manual(values = col) +

      # Labels
      xlab(NULL) +
      ylab("Shannon Diversity Index") +
      ggtitle(panel_title) +

      theme_bw() +

      # Theme customization
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      ) +

      # Sample size labels
      geom_text(data = sample_counts,
                aes(x = Location, y = y_pos, label = paste0("n=", n)),
                size = 5,
                fontface = "bold",
                inherit.aes = FALSE) +

      # P-value annotation
      geom_text(data = pvalue_df,
                aes(x = x, y = y_pos, label = p_label),
                size = 5,
                fontface = "bold",
                inherit.aes = FALSE)

    plot_list[[plot_idx]] <- p
    plot_idx <- plot_idx + 1
  }
}

# Arrange plots in grid (2 columns for Week.1 and Week.3, 3 rows for sample types)
ncols <- 2
nrows <- 3

fig <- arrangeGrob(
  grobs = plot_list,
  ncol = ncols,
  nrow = nrows,
  top = textGrob("Shannon Diversity by Location\n",
                 gp = gpar(fontsize = 24, fontface = "bold")),
  padding = unit(0.8, "line")
)

# Calculate dimensions and save
fig_width <- ncols * 4.5
fig_height <- nrows * 4

output_file <- file.path(RESULTS_DIR, "ShannonDiversity_Location.pdf")
ggsave(output_file,
       plot = fig,
       width = fig_width,
       height = fig_height,
       limitsize = FALSE)
cat("\nSaved plot to:", output_file, "\n")

# Also save as PNG
output_png <- file.path(RESULTS_DIR, "ShannonDiversity_Location.png")
ggsave(output_png,
       plot = fig,
       width = fig_width,
       height = fig_height,
       dpi = 300,
       limitsize = FALSE)
cat("Saved plot to:", output_png, "\n")

# Print statistical summary
cat("\n=== Statistical Summary (Wilcoxon rank-sum tests) ===\n")
for (week in c("Week.1", "Week.3")) {
  for (stype in c("Axilla", "Groin", "Stool")) {
    subset_data <- Diversity_filtered %>%
      filter(SampleCollectionWeek == week, SampleType == stype)

    if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
      ucmc_data <- subset_data$shannon[subset_data$Location == "UCMC"]
      zch_data <- subset_data$shannon[subset_data$Location == "ZCH"]

      if (length(ucmc_data) > 0 && length(zch_data) > 0) {
        test_result <- wilcox.test(ucmc_data, zch_data)
        cat(sprintf("%s - %s: p = %.4f (n_UCMC=%d, n_ZCH=%d)\n",
                    week, stype, test_result$p.value,
                    length(ucmc_data), length(zch_data)))
      }
    }
  }
}

cat("\nDone!\n")
