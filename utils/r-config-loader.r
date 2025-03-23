# utils/config.R - Configuration loader for R scripts

#' Load and parse configuration from a YAML file
#'
#' @param config_path Path to the YAML configuration file
#' @return A list containing the configuration
#' @import yaml
load_config <- function(config_path = "config.yaml") {
  # Check if yaml package is installed
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required. Please install it with 'install.packages(\"yaml\")'")
  }
  
  # Resolve script path
  script_path <- getwd()
  
  # If config_path is not absolute, make it relative to the script directory
  if (!file.exists(config_path)) {
    # Try relative to project root
    project_dir <- script_path
    config_path <- file.path(project_dir, config_path)
    
    # If still not found, look in parent directory (assuming utils subdir)
    if (!file.exists(config_path)) {
      project_dir <- dirname(script_path)
      config_path <- file.path(project_dir, config_path)
    }
  }
  
  # Check if the config file exists
  if (!file.exists(config_path)) {
    stop(paste0("Configuration file not found: ", config_path))
  }
  
  # Load the YAML config
  config <- yaml::read_yaml(config_path)
  
  # Resolve path variables
  config <- expand_path_variables(config)
  
  # Create directories if they don't exist
  if (!is.null(config$paths)) {
    for (path_name in names(config$paths)) {
      if (grepl("_dir$", path_name) && config$paths[[path_name]] != ".") {
        dir.create(config$paths[[path_name]], showWarnings = FALSE, recursive = TRUE)
      }
    }
  }
  
  return(config)
}

#' Resolve path variables in the configuration
#'
#' @param config The configuration list
#' @return The configuration list with resolved path variables
expand_path_variables <- function(config) {
  # First, collect all path variables
  path_vars <- list()
  
  # Start with paths section
  if (!is.null(config$paths)) {
    for (key in names(config$paths)) {
      path_vars[[key]] <- resolve_path(config$paths[[key]], path_vars)
      config$paths[[key]] <- path_vars[[key]]
    }
  }
  
  # Process all sections recursively
  config <- resolve_all_paths(config, path_vars)
  
  return(config)
}

#' Recursively resolve all path variables in the configuration
#'
#' @param config The configuration object (list, vector, or scalar)
#' @param path_vars The path variables lookup table
#' @return The configuration with resolved path variables
resolve_all_paths <- function(config, path_vars) {
  if (is.list(config)) {
    for (key in names(config)) {
      if (is.list(config[[key]]) || is.vector(config[[key]])) {
        config[[key]] <- resolve_all_paths(config[[key]], path_vars)
      } else if (is.character(config[[key]])) {
        config[[key]] <- resolve_path(config[[key]], path_vars)
      }
    }
  } else if (is.vector(config) && !is.null(names(config))) {
    for (i in seq_along(config)) {
      if (is.character(config[i])) {
        config[i] <- resolve_path(config[i], path_vars)
      }
    }
  } else if (is.vector(config) && is.character(config)) {
    for (i in seq_along(config)) {
      config[i] <- resolve_path(config[i], path_vars)
    }
  }
  
  return(config)
}

#' Resolve a single path string with variables
#'
#' @param path The path string
#' @param path_vars The path variables lookup table
#' @return The resolved path string
resolve_path <- function(path, path_vars) {
  if (!is.character(path)) {
    return(path)
  }
  
  # Replace ${variable} with its value
  pattern <- "\\${([^}]*)}"
  matches <- gregexpr(pattern, path)
  
  if (matches[[1]][1] != -1) {
    for (match_idx in seq_along(matches[[1]])) {
      match_pos <- matches[[1]][match_idx]
      match_len <- attr(matches[[1]], "match.length")[match_idx]
      
      var_name <- substr(path, match_pos + 2, match_pos + match_len - 2)
      
      if (var_name %in% names(path_vars)) {
        path <- sub(paste0("\\${", var_name, "}"), path_vars[[var_name]], path, fixed = TRUE)
      }
    }
  }
  
  return(path)
}

#' Get a nested value from the configuration
#'
#' @param config The configuration list
#' @param ... Key path elements
#' @return The value at the specified path
get_config_value <- function(config, ...) {
  keys <- list(...)
  
  # If the first argument is a character vector, use it as the keys
  if (length(keys) == 1 && is.character(keys[[1]]) && length(keys[[1]]) > 1) {
    keys <- as.list(keys[[1]])
  }
  
  current <- config
  for (key in keys) {
    if (is.null(current) || !is.list(current) || !(key %in% names(current))) {
      return(NULL)
    }
    current <- current[[key]]
  }
  
  return(current)
}

# Example usage
if (identical(environment(), globalenv())) {
  config <- load_config()
  cat("Loaded configuration:\n")
  cat("Project directory:", config$paths$project_dir, "\n")
  cat("Results directory:", config$paths$results_dir, "\n")
  cat("Key organisms:", paste(config$key_organisms, collapse=", "), "\n")
}
