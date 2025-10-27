#!/usr/bin/env Rscript
# Generate revised versions of Figure 2a and Figure 4 with sample counts added
# This script reads the processed data and generates publication-ready figures

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

# Set working directory
setwd("/home/david/projects/ZCH_UCMC_Manuscript/python_scripts")

# Create output directory
dir.create("revision_figures", showWarnings = FALSE)

# Load the saved R workspace from the main analysis
# Using the most recent workspace with starting dataframes (no .RData extension)
cat("Loading workspace...\n")
load("/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514")
cat("Workspace loaded successfully!\n")

cat("Generating revised figures with sample counts...\n")

# ============================================================================
# FIGURE 2a: Shannon Diversity - Location x SampleType x Week
# ============================================================================

cat("\n=== Figure 2a: Shannon Diversity ===\n")

# Filter data for Figure 2a
fig2a_data <- Diversity %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         Location %in% c("Cincinnati", "Hangzhou"),
         SampleCollectionWeek %in% c("Week.1", "Week.3")) %>%
  filter(!is.na(Shannon))

# Calculate comprehensive statistics for each group
sample_counts_2a <- fig2a_data %>%
  group_by(Location, SampleType, SampleCollectionWeek) %>%
  summarise(
    n = n(),
    median_shannon = median(Shannon, na.rm = TRUE),
    q25 = quantile(Shannon, 0.25, na.rm = TRUE),
    q75 = quantile(Shannon, 0.75, na.rm = TRUE),
    iqr = IQR(Shannon, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    median_iqr_label = sprintf("%.2f (%.2f-%.2f)", median_shannon, q25, q75)
  )

cat("\nSample sizes and statistics for Figure 2a:\n")
print(sample_counts_2a)

# Calculate p-values for ALL comparisons (between locations at each site/week)
cat("\n=== Calculating p-values for ALL comparisons ===\n")
pvalue_results <- data.frame()

for (week in c("Week.1", "Week.3")) {
  for (site in c("Axilla", "Groin", "Stool")) {
    subset_data <- fig2a_data %>%
      filter(SampleCollectionWeek == week, SampleType == site)

    if (nrow(subset_data) > 0 && length(unique(subset_data$Location)) == 2) {
      test_result <- wilcox.test(Shannon ~ Location, data = subset_data)

      cinci_data <- subset_data %>% filter(Location == "Cincinnati")
      hangz_data <- subset_data %>% filter(Location == "Hangzhou")

      pvalue_results <- rbind(pvalue_results, data.frame(
        Week = week,
        SampleType = site,
        n_Cincinnati = nrow(cinci_data),
        median_Cincinnati = median(cinci_data$Shannon, na.rm = TRUE),
        n_Hangzhou = nrow(hangz_data),
        median_Hangzhou = median(hangz_data$Shannon, na.rm = TRUE),
        p_value = test_result$p.value,
        p_label = ifelse(test_result$p.value < 0.001, "p<0.001",
                  ifelse(test_result$p.value < 0.01, sprintf("p=%.3f**", test_result$p.value),
                  ifelse(test_result$p.value < 0.05, sprintf("p=%.3f*", test_result$p.value),
                         sprintf("p=%.3f", test_result$p.value))))
      ))
    }
  }
}

cat("\nAll p-values for Figure 2a:\n")
print(pvalue_results)

# Save comprehensive statistics table
write.csv(pvalue_results,
          "revision_figures/Figure2a_all_pvalues.csv",
          row.names = FALSE)

# Create the plot with sample counts
# Calculate positions for sample count labels
count_positions <- sample_counts_2a %>%
  mutate(y_pos = median_shannon * 1.15)  # Position labels above boxes

fig2a <- ggplot(fig2a_data,
                aes(x = Location, y = Shannon, fill = Location)) +
  geom_boxplot(lwd = 1, aes(color = factor(Location), fill = NA),
               outlier.size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
  geom_point(size = 4, aes(color = factor(Location))) +
  # Add sample count labels using pre-calculated data
  geom_text(data = count_positions,
            aes(x = Location, y = y_pos, label = paste0("n=", n)),
            size = 5,
            fontface = "bold",
            inherit.aes = FALSE) +
  facet_grid(SampleCollectionWeek ~ SampleType) +
  ylab("Shannon Diversity Index\n") +
  xlab(NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave("revision_figures/Figure2a_Shannon_with_sample_counts.pdf",
       plot = fig2a, width = 12, height = 10, limitsize = FALSE)

cat("\n✓ Saved Figure 2a with sample counts\n")

# Create enhanced version with ALL p-values annotated
# Prepare p-value positions for annotation
pvalue_positions <- pvalue_results %>%
  mutate(
    SampleCollectionWeek = Week,
    y_pos = pmax(median_Cincinnati, median_Hangzhou) * 1.25,
    x_pos = 1.5  # Between Cincinnati and Hangzhou
  )

fig2a_enhanced <- ggplot(fig2a_data,
                aes(x = Location, y = Shannon, fill = Location)) +
  geom_boxplot(lwd = 1, aes(color = factor(Location), fill = NA),
               outlier.size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 8) +
  geom_point(size = 4, aes(color = factor(Location))) +
  # Add sample count labels
  geom_text(data = count_positions,
            aes(x = Location, y = y_pos, label = paste0("n=", n)),
            size = 5,
            fontface = "bold",
            inherit.aes = FALSE) +
  # Add p-value annotations
  geom_text(data = pvalue_positions,
            aes(x = x_pos, y = y_pos, label = p_label),
            size = 4.5,
            fontface = "bold",
            color = "black",
            inherit.aes = FALSE) +
  facet_grid(SampleCollectionWeek ~ SampleType) +
  ylab("Shannon Diversity Index\n") +
  xlab(NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave("revision_figures/Figure2a_Shannon_with_all_stats.pdf",
       plot = fig2a_enhanced, width = 12, height = 10, limitsize = FALSE)

cat("✓ Saved Figure 2a with ALL p-values annotated\n")

# Alternative version with counts in strip labels
fig2a_data_labeled <- fig2a_data %>%
  group_by(SampleCollectionWeek, SampleType) %>%
  mutate(facet_label = paste0(SampleType, " (Week ",
                               gsub("Week.", "", SampleCollectionWeek),
                               ")")) %>%
  ungroup()

# Save sample count table
write.csv(sample_counts_2a,
          "revision_figures/Figure2a_sample_counts.csv",
          row.names = FALSE)

# ============================================================================
# FIGURE 4: Longitudinal Changes in BSI-Associated Species
# ============================================================================

cat("\n=== Figure 4: Longitudinal Changes ===\n")

# Key organisms for Figure 4 (from the figure legend)
# The legend mentions 8 species but only lists 5 (a-e)
# Adding 3 more based on BSI epidemiology analysis (Table 2 significant organisms)
key_organisms <- c("Staphylococcus.aureus",        # a
                   "Klebsiella.pneumoniae",        # b
                   "Klebsiella.oxytoca",           # c
                   "Enterococcus.faecium",         # d
                   "Enterococcus.faecalis",        # e
                   "Serratia.marcescens",          # f (6.2% UCMC vs 0.8% ZCH)
                   "Escherichia.coli",             # g (common pathogen)
                   "Streptococcus.pyogenes")       # h (0% UCMC vs 4.2% ZCH)

cat("\nGenerating box plots for", length(key_organisms), "organisms (Figure 4 parts a-h)\n")

# Prepare data
fig4_data <- NICUSpeciesNR %>%
  filter(SampleType %in% c("Axilla", "Groin", "Stool"),
         SampleCollectionWeek %in% c("Week.1", "Week.3"))

# Function to create organism-specific plot with sample counts
create_organism_plot_with_counts <- function(data, organism, title) {

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
    facet_grid(SampleType ~ Location) +
    scale_y_log10() +
    ylab(paste0(gsub("\\.", " ", title), "\n(log10 abundance)")) +
    xlab(NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14, face = "bold.italic"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  return(list(plot = p, counts = sample_counts))
}

# Generate plots for each organism
organism_plots <- list()
all_counts <- list()

for (org in key_organisms) {
  result <- create_organism_plot_with_counts(fig4_data, org, org)
  organism_plots[[org]] <- result$plot
  all_counts[[org]] <- result$counts

  # Save individual plot
  clean_name <- gsub("\\.", "_", org)
  ggsave(paste0("revision_figures/Figure4_", clean_name, "_with_counts.pdf"),
         plot = result$plot, width = 10, height = 8)

  cat("✓ Saved", org, "plot with sample counts\n")
}

# Combine all sample counts into one table
combined_counts <- bind_rows(all_counts, .id = "Organism")
write.csv(combined_counts,
          "revision_figures/Figure4_all_sample_counts.csv",
          row.names = FALSE)

# ============================================================================
# Create comprehensive combined Figure 4 with all 8 species
# ============================================================================

cat("\n=== Creating comprehensive Figure 4 with all 8 species ===\n")

library(gridExtra)
library(grid)

# Create a combined plot with all 8 organisms in a grid layout
# We'll use 4 rows x 2 columns for better readability

# Simplify the individual plots for the combined figure
create_organism_plot_simplified <- function(data, organism, panel_label) {

  plot_data <- data %>%
    select(Sample, Location, SampleType, SampleCollectionWeek, all_of(organism)) %>%
    rename(abundance = all_of(organism)) %>%
    filter(!is.na(abundance)) %>%
    mutate(abundance_log = ifelse(abundance == 0, 1, abundance))

  # Calculate sample sizes for annotation
  sample_counts <- plot_data %>%
    group_by(Location, SampleType, SampleCollectionWeek) %>%
    summarise(n = n(), .groups = 'drop') %>%
    group_by(Location, SampleCollectionWeek) %>%
    summarise(total_n = sum(n), .groups = 'drop') %>%
    mutate(label = paste0("n=", total_n))

  # Create simplified plot
  p <- ggplot(plot_data,
              aes(x = SampleCollectionWeek, y = abundance_log, fill = Location)) +
    geom_boxplot(lwd = 0.5, outlier.size = 1, alpha = 0.7) +
    facet_grid(SampleType ~ Location, scales = "free_y") +
    scale_y_log10() +
    scale_fill_manual(values = c("Cincinnati" = "#E69F00", "Hangzhou" = "#56B4E9")) +
    ylab(NULL) +
    xlab(NULL) +
    ggtitle(paste0(panel_label, ". ", gsub("\\.", " ", organism))) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 11, face = "bold.italic", hjust = 0),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    )

  return(p)
}

# Create all 8 plots
panel_labels <- c("a", "b", "c", "d", "e", "f", "g", "h")
simplified_plots <- list()

for (i in seq_along(key_organisms)) {
  org <- key_organisms[i]
  label <- panel_labels[i]
  simplified_plots[[i]] <- create_organism_plot_simplified(fig4_data, org, label)
  cat("Created panel", label, "for", org, "\n")
}

# Arrange all plots in a grid (4 rows x 2 columns)
combined_fig4 <- arrangeGrob(
  grobs = simplified_plots,
  ncol = 2,
  nrow = 4,
  top = textGrob("Figure 4. Longitudinal changes in BSI-associated species abundance\n",
                 gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob("Log10 relative abundance shown for each body site (Axilla, Groin, Stool) by location and sample collection week",
                    gp = gpar(fontsize = 9, fontface = "italic"))
)

# Save the comprehensive figure
ggsave("revision_figures/Figure4_comprehensive_all_8_species.pdf",
       plot = combined_fig4,
       width = 16,
       height = 20,
       limitsize = FALSE)

cat("\n✓ Saved comprehensive Figure 4 with all 8 species!\n")

# ============================================================================
# Statistical Tests with Sample Sizes Reported
# ============================================================================

cat("\n=== Statistical Tests ===\n")

# For Figure 2a - test differences between locations at each timepoint/site
cat("\nFigure 2a - Shannon diversity comparisons:\n")

for (week in c("Week.1", "Week.3")) {
  for (site in c("Axilla", "Groin", "Stool")) {
    subset_data <- fig2a_data %>%
      filter(SampleCollectionWeek == week, SampleType == site)

    if (nrow(subset_data) > 0) {
      test_result <- wilcox.test(Shannon ~ Location, data = subset_data)

      n_cinci <- sum(subset_data$Location == "Cincinnati")
      n_hangz <- sum(subset_data$Location == "Hangzhou")

      cat(sprintf("%s %s: Cincinnati (n=%d) vs Hangzhou (n=%d), p=%.4f %s\n",
                  week, site, n_cinci, n_hangz, test_result$p.value,
                  ifelse(test_result$p.value < 0.05, "*", "")))
    }
  }
}

cat("\n✓ All revised figures with sample counts have been generated!\n")
cat("✓ Sample count tables saved to revision_figures/\n")
cat("\nOutput files:\n")
cat("- revision_figures/Figure2a_Shannon_with_sample_counts.pdf\n")
cat("- revision_figures/Figure2a_sample_counts.csv\n")
cat("\n*** COMPREHENSIVE FIGURE 4 (recommended for manuscript): ***\n")
cat("- revision_figures/Figure4_comprehensive_all_8_species.pdf (4x2 grid, all species)\n")
cat("\nIndividual Figure 4 panels (for supplementary or detailed review):\n")
for (org in key_organisms) {
  clean_name <- gsub("\\.", "_", org)
  cat(paste0("- revision_figures/Figure4_", clean_name, "_with_counts.pdf\n"))
}
cat("- revision_figures/Figure4_all_sample_counts.csv\n")
