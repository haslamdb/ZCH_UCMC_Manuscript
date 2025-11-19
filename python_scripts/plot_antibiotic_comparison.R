#!/usr/bin/env Rscript
# Compare antibiotic usage between UCMC and ZCH
# Uses same styling as Figure 4

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(ggthemes)

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript")

# Create output directory
dir.create("revision_figures", showWarnings = FALSE)

cat("\n=== Loading Antibiotic Data ===\n")

# Load data
uc_data <- read.csv("UC_antibiotics_parsed.csv")
zch_data <- read.csv("ZCH_antibiotics_parsed.csv")

# Add location column
uc_data$Location <- "UCMC"
zch_data$Location <- "ZCH"

# Standardize MRN column names
colnames(uc_data)[colnames(uc_data) == "mrn_v2"] <- "MRN"
colnames(zch_data)[colnames(zch_data) == "mrn_v2_zju"] <- "MRN"

# Get antibiotic columns (exclude record_id, MRN, Location)
antibiotic_cols <- grep("^(Ampicillin|Penicillin|Nafcillin|Cefotaxime|Moxalactam|Ceftazidime|Ceftriaxone|Cefepime|Gentamicin|Tobramycin|Vancomycin|Fluconazole|Amphotericin|Meropenem|Azithromycin|Piperacillin)",
                        colnames(uc_data), value = TRUE)

# Separate week 1 and week 2 (which is actually week 3 data) columns
w1_cols <- grep("_w1$", antibiotic_cols, value = TRUE)
w2_cols <- grep("_w2$", antibiotic_cols, value = TRUE)

cat(sprintf("Week 1 columns: %d\n", length(w1_cols)))
cat(sprintf("Week 2 (week 3 data) columns: %d\n", length(w2_cols)))

# Create cumulative week 3 columns for both datasets
cat("\n=== Creating cumulative week 3 columns ===\n")
for (i in seq_along(w1_cols)) {
  w1_col <- w1_cols[i]
  # Find matching w2 column (replace _w1 with _w2)
  w2_col <- gsub("_w1$", "_w2", w1_col)

  if (w2_col %in% colnames(uc_data)) {
    # Create cumulative column name
    cumulative_col <- gsub("_w1$", "_cumulative_w3", w1_col)

    # Add cumulative columns to UC data
    uc_data[[cumulative_col]] <- uc_data[[w1_col]] + uc_data[[w2_col]]

    # Add cumulative columns to ZCH data
    zch_data[[cumulative_col]] <- zch_data[[w1_col]] + zch_data[[w2_col]]
  }
}

# Get updated antibiotic columns including cumulative
antibiotic_cols <- grep("^(Ampicillin|Penicillin|Nafcillin|Cefotaxime|Moxalactam|Ceftazidime|Ceftriaxone|Cefepime|Gentamicin|Tobramycin|Vancomycin|Fluconazole|Amphotericin|Meropenem|Azithromycin|Piperacillin)",
                        colnames(uc_data), value = TRUE)

# Now separate into week 1 and cumulative week 3 columns
w1_cols <- grep("_w1$", antibiotic_cols, value = TRUE)
cumulative_w3_cols <- grep("_cumulative_w3$", antibiotic_cols, value = TRUE)

cat(sprintf("Created %d cumulative week 3 columns\n", length(cumulative_w3_cols)))

# Combine datasets
combined_data <- bind_rows(
  uc_data[, c("MRN", "Location", antibiotic_cols)],
  zch_data[, c("MRN", "Location", antibiotic_cols)]
)

cat(sprintf("Combined data: %d UCMC + %d ZCH = %d total patients\n",
            nrow(uc_data), nrow(zch_data), nrow(combined_data)))

# Define color palette (same as Figure 4)
col <- c("#E69F00", "#56B4E9")  # UCMC (orange), ZCH (blue)

# Function to create individual antibiotic plot (Figure 4 style)
create_antibiotic_plot <- function(data, antibiotic, panel_label) {

  # Reshape data for this antibiotic
  # Keep ALL patients (including zeros) so both locations always appear
  plot_data <- data %>%
    select(MRN, Location, all_of(antibiotic)) %>%
    rename(days = all_of(antibiotic)) %>%
    filter(!is.na(days))

  if (sum(plot_data$days) == 0) {
    return(NULL)  # Skip if no usage at all
  }

  # For plotting, we'll still only show points for days > 0
  plot_data_points <- plot_data %>% filter(days > 0)

  # Calculate sample sizes (only count patients who received it)
  sample_counts <- plot_data_points %>%
    group_by(Location) %>%
    summarise(
      n = n(),
      max_days = max(days, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(y_pos = max_days * 1.3)  # Position above boxes

  # Calculate p-value (Wilcoxon test)
  ucmc_days <- plot_data %>% filter(Location == "UCMC") %>% pull(days)
  zch_days <- plot_data %>% filter(Location == "ZCH") %>% pull(days)

  # For the full comparison (including zeros from original data)
  ucmc_all <- data %>% filter(Location == "UCMC") %>% pull(!!antibiotic)
  zch_all <- data %>% filter(Location == "ZCH") %>% pull(!!antibiotic)

  if (length(ucmc_all) > 0 && length(zch_all) > 0) {
    test_result <- tryCatch({
      wilcox.test(ucmc_all, zch_all)
    }, error = function(e) {
      list(p.value = NA)
    })

    p_value <- test_result$p.value
    max_y <- max(plot_data$days, na.rm = TRUE)

    p_label <- ifelse(is.na(p_value), "",
                ifelse(p_value < 0.001, "p<0.001***",
                ifelse(p_value < 0.01, sprintf("p=%.3f**", p_value),
                ifelse(p_value < 0.05, sprintf("p=%.3f*", p_value),
                       sprintf("p=%.2f", p_value)))))

    pvalue_df <- data.frame(
      x = 1.5,
      y_pos = max_y * 1.5,
      p_label = p_label
    )
  } else {
    pvalue_df <- data.frame(x = numeric(0), y_pos = numeric(0), p_label = character(0))
  }

  # Get actual counts
  ucmc_count <- sample_counts %>% filter(Location == "UCMC") %>% pull(n)
  zch_count <- sample_counts %>% filter(Location == "ZCH") %>% pull(n)
  ucmc_count <- ifelse(length(ucmc_count) == 0, 0, ucmc_count)
  zch_count <- ifelse(length(zch_count) == 0, 0, zch_count)

  cat(sprintf("  %s: UCMC n=%d, ZCH n=%d, p=%s\n",
              antibiotic,
              ucmc_count,
              zch_count,
              if(nrow(pvalue_df) > 0) pvalue_df$p_label else "N/A"))

  # Clean antibiotic name for display
  # Convert _w1 to _wk1 and _cumulative_w3 to (cumulative wk3)
  clean_name <- gsub("_cumulative_w3", " (cumulative wk3)", antibiotic)
  clean_name <- gsub("_w1", " (wk1)", clean_name)
  clean_name <- gsub("_", " ", clean_name)

  # Special handling for Piperacillin Tazobactam - put on two lines
  if (grepl("Piperacillin", clean_name)) {
    clean_name <- gsub("Piperacillin Tazobactam", "Piperacillin\nTazobactam", clean_name)
  }

  # Create plot matching Figure 4 style
  # Use plot_data (includes zeros) for boxplot, but plot_data_points for individual points
  p <- ggplot(plot_data, aes(x = Location, y = days)) +
    geom_boxplot(lwd = 1, aes(color = factor(Location)),
                 fill = NA,
                 outlier.shape = NA) +  # Don't show outliers, we'll add points separately
    stat_summary(fun = mean, geom = "point", shape = 5, size = 4.8,
                data = plot_data_points) +  # Only show mean for non-zero
    scale_colour_manual(values = col) +
    geom_point(data = plot_data_points, size = 2.4,
              aes(color = factor(Location))) +  # Only show points for non-zero
    xlab(NULL) +
    ylab(clean_name) +
    scale_color_tableau() +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 20, face = "bold"),  # Antibiotic name - increased to 20
      axis.text.x = element_text(size = 18),  # Location labels (UCMC/ZCH)
      axis.text.y = element_text(size = 18),  # Y-axis tick labels
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    # Add sample count labels (increased size from 4 to 7)
    geom_text(data = sample_counts,
              aes(x = Location, y = y_pos, label = paste0("n=", n)),
              size = 7,
              fontface = "bold",
              inherit.aes = FALSE)

  # Add p-value if available (increased size from 4 to 7)
  if (nrow(pvalue_df) > 0 && pvalue_df$p_label != "") {
    p <- p + geom_text(data = pvalue_df,
                      aes(x = x, y = y_pos, label = p_label),
                      size = 7,
                      fontface = "bold",
                      inherit.aes = FALSE)
  }

  # Add panel label
  p <- p + ggtitle(paste0(panel_label, "."))

  return(p)
}

# Identify antibiotics with usage - separate by week
cat("\n=== Identifying antibiotics with usage ===\n")

# Week 1 antibiotics
w1_with_usage <- c()
for (antibiotic in w1_cols) {
  total_usage <- sum(combined_data[[antibiotic]] > 0, na.rm = TRUE)
  if (total_usage > 0) {
    w1_with_usage <- c(w1_with_usage, antibiotic)
  }
}

# Cumulative Week 3 antibiotics
w3_with_usage <- c()
for (antibiotic in cumulative_w3_cols) {
  total_usage <- sum(combined_data[[antibiotic]] > 0, na.rm = TRUE)
  if (total_usage > 0) {
    w3_with_usage <- c(w3_with_usage, antibiotic)
  }
}

cat(sprintf("Week 1: Found %d/%d antibiotics with usage\n",
            length(w1_with_usage), length(w1_cols)))
cat(sprintf("Cumulative Week 3: Found %d/%d antibiotics with usage\n",
            length(w3_with_usage), length(cumulative_w3_cols)))

# Generate plots for Week 1
cat("\n=== Generating Week 1 plots ===\n")
w1_plot_list <- list()
w1_panel_labels <- letters[1:length(w1_with_usage)]

for (i in seq_along(w1_with_usage)) {
  antibiotic <- w1_with_usage[i]
  panel_label <- w1_panel_labels[i]

  p <- create_antibiotic_plot(combined_data, antibiotic, panel_label)
  if (!is.null(p)) {
    w1_plot_list[[length(w1_plot_list) + 1]] <- p
  }
}

# Generate plots for Cumulative Week 3
cat("\n=== Generating Cumulative Week 3 plots ===\n")
w3_plot_list <- list()
w3_panel_labels <- letters[1:length(w3_with_usage)]

for (i in seq_along(w3_with_usage)) {
  antibiotic <- w3_with_usage[i]
  panel_label <- w3_panel_labels[i]

  p <- create_antibiotic_plot(combined_data, antibiotic, panel_label)
  if (!is.null(p)) {
    w3_plot_list[[length(w3_plot_list) + 1]] <- p
  }
}

# Create Week 1 figure
cat("\n=== Creating Week 1 figure ===\n")

n_plots_w1 <- length(w1_plot_list)
ncols <- 4
nrows_w1 <- ceiling(n_plots_w1 / ncols)

w1_fig <- arrangeGrob(
  grobs = w1_plot_list,
  ncol = ncols,
  nrow = nrows_w1,
  top = textGrob("Supplementary Figure A. Antibiotic Usage Comparison (Week 1): UCMC vs ZCH\n",
                 gp = gpar(fontsize = 24, fontface = "bold")),
  padding = unit(0.8, "line"),
  widths = unit(rep(3.78, ncols), "in"),
  heights = unit(rep(3.5, nrows_w1), "in")
)

# Calculate figure dimensions for Week 1
fig_width_w1 <- ncols * 4.5
fig_height_w1 <- nrows_w1 * 4

# Save Week 1 figure
ggsave("revision_figures/Antibiotic_Comparison_UCMC_vs_ZCH_Week1.pdf",
       plot = w1_fig,
       width = fig_width_w1,
       height = fig_height_w1,
       limitsize = FALSE)

cat("✓ Saved Week 1 antibiotic comparison figure\n")

# Create Cumulative Week 3 figure
cat("\n=== Creating Cumulative Week 3 figure ===\n")

n_plots_w3 <- length(w3_plot_list)
nrows_w3 <- ceiling(n_plots_w3 / ncols)

w3_fig <- arrangeGrob(
  grobs = w3_plot_list,
  ncol = ncols,
  nrow = nrows_w3,
  top = textGrob("Supplementary Figure B. Antibiotic Usage Comparison (Cumulative Week 3): UCMC vs ZCH\n",
                 gp = gpar(fontsize = 24, fontface = "bold")),
  padding = unit(0.8, "line"),
  widths = unit(rep(3.78, ncols), "in"),
  heights = unit(rep(3.5, nrows_w3), "in")
)

# Calculate figure dimensions for Week 3
fig_width_w3 <- ncols * 4.5
fig_height_w3 <- nrows_w3 * 4

# Save Week 3 figure
ggsave("revision_figures/Antibiotic_Comparison_UCMC_vs_ZCH_CumulativeWeek3.pdf",
       plot = w3_fig,
       width = fig_width_w3,
       height = fig_height_w3,
       limitsize = FALSE)

cat("✓ Saved Cumulative Week 3 antibiotic comparison figure\n")

# ============================================================================
# Generate statistics table
# ============================================================================

cat("\n=== Calculating statistics ===\n")

stats_results <- data.frame()

# Combine week 1 and cumulative week 3 antibiotics for statistics
all_antibiotics_for_stats <- c(w1_with_usage, w3_with_usage)

for (antibiotic in all_antibiotics_for_stats) {
  ucmc_data <- combined_data %>% filter(Location == "UCMC") %>% pull(!!antibiotic)
  zch_data <- combined_data %>% filter(Location == "ZCH") %>% pull(!!antibiotic)

  # Prevalence
  ucmc_prev <- sum(ucmc_data > 0) / length(ucmc_data) * 100
  zch_prev <- sum(zch_data > 0) / length(zch_data) * 100

  # Mean days (among those who received)
  ucmc_mean <- mean(ucmc_data[ucmc_data > 0])
  zch_mean <- mean(zch_data[zch_data > 0])

  # Count
  ucmc_n <- sum(ucmc_data > 0)
  zch_n <- sum(zch_data > 0)

  # Wilcoxon test
  if (length(ucmc_data) > 0 && length(zch_data) > 0) {
    test_result <- tryCatch({
      wilcox.test(ucmc_data, zch_data)
    }, error = function(e) {
      list(p.value = NA)
    })
    p_value <- test_result$p.value
  } else {
    p_value <- NA
  }

  stats_results <- rbind(stats_results, data.frame(
    Antibiotic = antibiotic,
    UCMC_N = ucmc_n,
    UCMC_Prevalence_Percent = ucmc_prev,
    UCMC_Mean_Days = ifelse(is.nan(ucmc_mean), 0, ucmc_mean),
    ZCH_N = zch_n,
    ZCH_Prevalence_Percent = zch_prev,
    ZCH_Mean_Days = ifelse(is.nan(zch_mean), 0, zch_mean),
    P_value = p_value,
    Significant = ifelse(!is.na(p_value) && p_value < 0.05, "Yes", "No")
  ))
}

# Sort by p-value
stats_results <- stats_results %>% arrange(P_value)

# Save statistics
write.csv(stats_results,
          "revision_figures/Antibiotic_Statistics_UCMC_vs_ZCH.csv",
          row.names = FALSE)

cat("✓ Saved statistics table\n")

# Print summary
cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("ANTIBIOTIC COMPARISON COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat(sprintf("\nTotal patients: %d (%d UCMC, %d ZCH)\n",
            nrow(combined_data),
            sum(combined_data$Location == "UCMC"),
            sum(combined_data$Location == "ZCH")))
cat(sprintf("Week 1 antibiotics with usage: %d/%d\n",
            length(w1_with_usage), length(w1_cols)))
cat(sprintf("Cumulative Week 3 antibiotics with usage: %d/%d\n",
            length(w3_with_usage), length(cumulative_w3_cols)))
cat(sprintf("Significant differences (p<0.05): %d\n",
            sum(stats_results$P_value < 0.05, na.rm = TRUE)))

cat("\nGenerated files:\n")
cat("  - revision_figures/Antibiotic_Comparison_UCMC_vs_ZCH_Week1.pdf\n")
cat("  - revision_figures/Antibiotic_Comparison_UCMC_vs_ZCH_CumulativeWeek3.pdf\n")
cat("  - revision_figures/Antibiotic_Statistics_UCMC_vs_ZCH.csv\n")

cat("\nTop significant differences:\n")
sig_results <- stats_results %>%
  filter(P_value < 0.05) %>%
  head(10)
for (i in 1:min(10, nrow(sig_results))) {
  cat(sprintf("  %s: UCMC %.1f%% vs ZCH %.1f%% (p=%.4f)\n",
              sig_results$Antibiotic[i],
              sig_results$UCMC_Prevalence_Percent[i],
              sig_results$ZCH_Prevalence_Percent[i],
              sig_results$P_value[i]))
}

cat(paste(rep("=", 70), collapse = ""))
cat("\n")
