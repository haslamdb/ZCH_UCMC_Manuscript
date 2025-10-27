# Helper function to add sample counts to ggplot2 boxplots
# This can be sourced in the main analysis script

# Function to add sample counts above box plots
add_n_to_boxplot <- function(ggplot_obj, data, x_var, facet_vars = NULL) {
  # Calculate sample sizes
  if (is.null(facet_vars)) {
    n_df <- data %>%
      group_by(across(all_of(x_var))) %>%
      summarise(n = n(), .groups = 'drop')
  } else {
    n_df <- data %>%
      group_by(across(all_of(c(x_var, facet_vars)))) %>%
      summarise(n = n(), .groups = 'drop')
  }

  # Add sample counts to plot
  ggplot_obj +
    geom_text(data = n_df,
              aes(label = paste0("n=", n)),
              y = Inf,
              vjust = 1.5,
              size = 4,
              fontface = "bold")
}

# Alternative function that adds counts in the x-axis labels
add_n_to_labels <- function(ggplot_obj, data, x_var) {
  # Calculate sample sizes
  n_df <- data %>%
    group_by(across(all_of(x_var))) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(label = paste0(get(x_var), "\n(n=", n, ")"))

  # Get the x variable levels with counts
  labels <- setNames(n_df$label, n_df[[x_var]])

  # Update x-axis labels
  ggplot_obj +
    scale_x_discrete(labels = labels)
}

# Function to create a summary table of sample sizes for faceted plots
create_sample_size_table <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = tail(group_vars, 1), values_from = n, values_fill = 0)
}

# Function to add sample counts to facet labels (most useful for Figure 2a and Figure 4)
# This modifies the data to include counts in the facet variable labels
add_n_to_facet_labels <- function(data, facet_var) {
  data %>%
    group_by(across(all_of(facet_var))) %>%
    mutate(!!facet_var := paste0(get(facet_var), "\n(n=", n(), ")")) %>%
    ungroup()
}

# Comprehensive function for adding sample counts to any combination
add_sample_counts_comprehensive <- function(data, x_var, y_var, facet_row = NULL, facet_col = NULL) {
  library(dplyr)
  library(tidyr)

  # Build grouping variables
  group_vars <- c(x_var)
  if (!is.null(facet_row)) group_vars <- c(group_vars, facet_row)
  if (!is.null(facet_col)) group_vars <- c(group_vars, facet_col)

  # Calculate sample sizes
  n_df <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n = n(),
      y_max = max(get(y_var), na.rm = TRUE),
      .groups = 'drop'
    )

  return(n_df)
}
