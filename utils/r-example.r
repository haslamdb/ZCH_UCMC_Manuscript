# Example R script using the unified configuration
# examples/r_example.R

# Source the configuration loader
source("../utils/config.R")

# Load required libraries
library(tidyverse)
library(ggplot2)

# Main function
main <- function() {
  # Load the configuration
  config <- load_config()
  
  # Access configuration values
  cat("Loading data from:", config$input_files$microbiome_counts, "\n")
  cat("Saving results to:", config$paths$results_dir, "\n")
  cat("Key organisms to analyze:", paste(config$key_organisms, collapse=", "), "\n")
  
  # Set random seed for reproducibility
  set.seed(config$statistics$random_seed)
  
  # Load data
  tryCatch({
    microbiome_df <- read.csv(config$input_files$microbiome_counts, row.names=1)
    cat("Loaded microbiome data with dimensions:", dim(microbiome_df)[1], "x", dim(microbiome_df)[2], "\n")
    
    metadata_df <- read.csv(config$input_files$sample_key, row.names=1)
    cat("Loaded metadata with dimensions:", dim(metadata_df)[1], "x", dim(metadata_df)[2], "\n")
    
    # Show a summary of key organisms
    key_organisms <- config$key_organisms
    organisms_found <- key_organisms[key_organisms %in% colnames(microbiome_df)]
    
    cat("\nFound", length(organisms_found), "out of", length(key_organisms), "key organisms in the dataset:\n")
    for (organism in organisms_found) {
      non_zero <- sum(microbiome_df[, organism] > 0)
      percent <- non_zero / nrow(microbiome_df) * 100
      cat("  -", organism, ": Present in", non_zero, "samples (", round(percent, 1), "%)\n")
    }
    
    # Create a simple visualization using the config's visualization parameters
    fig_size <- config$visualization$figure_sizes$default
    dpi <- config$visualization$dpi
    color_palette <- config$visualization$color_palettes$categorical
    
    # Example: Create a bar plot of organism prevalence
    prevalence <- sapply(organisms_found, function(org) mean(microbiome_df[, org] > 0) * 100)
    prevalence_df <- data.frame(
      Organism = organisms_found,
      Prevalence = prevalence
    )
    
    # Sort by prevalence
    prevalence_df <- prevalence_df[order(-prevalence_df$Prevalence), ]
    prevalence_df$Organism <- factor(prevalence_df$Organism, levels = prevalence_df$Organism)
    
    # Create plot
    p <- ggplot(prevalence_df, aes(x = Organism, y = Prevalence, fill = Organism)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = color_palette[1:length(organisms_found)]) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      labs(
        title = "Prevalence of Key Organisms",
        y = "Prevalence (%)",
        x = NULL
      )
    
    # Save the figure to the results directory
    figures_dir <- config$paths$figures_dir
    dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
    
    output_file <- file.path(figures_dir, "organism_prevalence_r.png")
    ggsave(output_file, plot = p, width = fig_size[1], height = fig_size[2], dpi = dpi)
    
    cat("\nCreated visualization:", output_file, "\n")
    
  }, error = function(e) {
    cat("Error processing data:", conditionMessage(e), "\n")
  })
}

# Run the main function
main()
