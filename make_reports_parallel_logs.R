#!/usr/bin/env Rscript
# ✅ Required libraries for make_reports_parallel.R and report_util.R
required_packages <- c(
  "tidyverse", "rmarkdown", "furrr", "progressr", "future", "tictoc", "readr", "yaml",
  "stats", "multcomp", "dplyr", "ggplot2", "scales", "ggsignif", "plotly", "nnet"
)

# Install any missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing_packages) > 0) {
  message("📦 Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

# Load packages silently
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

 # read arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || "--help" %in% args) {
  cat("
Usage:
  Rscript make_reports.R <config.yaml>

Arguments:
  config.yaml   A YAML file containing all required and optional input paths and parameters.

Example config.yaml:

comparisons_file: Metadata/comparisons.csv
samples_file: Metadata/sorted_LLNEXT_samples_binary.csv
exist_file: Data/exist.csv
timepoints_file: Metadata/LLNext_ind_timepoints.csv Optional
extra_cols: [Sex, Age] Optional
output_dir: reports Default
\n")
  quit(status = 0)
}

# Read config file
config_path <- args[1]
if (!file.exists(config_path)) stop("❌ Config file not found: ", config_path)

config <- yaml::read_yaml(config_path)

# Extract values from config
cmp_file        <- config$comparisons_file
samples_file    <- config$samples_file
exist_file      <- config$exist_file
timepoints_file <- config$timepoints_file %||% NULL  # allow missing
extra_cols      <- config$extra_cols %||% character()  # allow missing
output_dir      <- config$output_dir %||% "reports"


# Validate required files
required <- list(cmp_file, samples_file, exist_file)
names(required) <- c("comparisons_file", "samples_file", "exist_file") 

missing <- names(required)[!file.exists(unlist(required))]

if (length(missing)) {
  stop(paste("❌ Missing required files:", paste(missing, collapse = ", ")))
}

# Validate timepoints file if provided
if (!is.null(timepoints_file) && !file.exists(timepoints_file)) {
  stop("❌ Timepoints file does not exist: ", timepoints_file)
}


# read inputs
comparisons <- read_csv(cmp_file)
samples     <- read_csv(samples_file) %>%
  rename_with(~ "SampleName", .cols = 1)
exist       <- read.csv(exist_file, header = TRUE, row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
extra_syms <- syms(extra_cols)

# Get the directory of this script
args_full <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(normalizePath(sub("--file=", "", args_full[grep("--file=", args_full)])))
template <- file.path(script_path, "template/template_phipseq.Rmd")
library_meta <- file.path(script_path, "library_meta/all_libraries_with_important_info.rds")
lib_metadata_df <- readRDS(library_meta)



# Set a palette of up to 12 distinct colors
available_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "bisque2", "#a6761d", "#666666",
  "#a6cee3", "#1f78b4", "#b2df8a", "salmon2")

# Color map to remember assignments
group_color_map <- list()


# Use all available cores (or adjust to your needs)
plan(multisession, workers = 4)  # #plan(multisession, workers = parallel::detectCores())  # use multicore on Unix, or multisession for cross-platform

# Create output directory if needed
is_absolute_path <- function(path) grepl("^(/|[A-Za-z]:)", path)
output_dir <- if (is_absolute_path(output_dir)) output_dir else file.path(getwd(), output_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)
out_tables <- file.path(output_dir, "Tables")
if (!dir.exists(out_tables)) dir.create(out_tables)

# Create a log file and initialize output tracking
log_file <- "rendering_log.txt"
summary_file <- "render_summary.csv"
writeLines("Rendering Log\n==============", con = log_file)

# Run rendering in parallel and track summary
results <- future_pmap_dfr(
  list(name = comparisons[[1]], g1 = comparisons[[2]], g2 = comparisons[[3]]),
  function(name, g1, g2) {
    message(sprintf("Rendering %s...", name))
    tic(name)
    
    
    tryCatch({
      # Color handling (unchanged)
      groups <- c(g1, g2)
      for (g in groups) {
        if (!g %in% names(group_color_map)) {
          used_colors <- unlist(group_color_map)
          available <- setdiff(available_colors, used_colors)
          if (length(available) == 0) {
            warning("Ran out of distinct colors, reusing colors!")
            available <- available_colors
          }
          group_color_map[[g]] <- available[1]
        }
      }
      custom_colors <- setNames(unlist(group_color_map[groups]), groups)
        
      # Metadata filtering
      meta_i <- samples %>%
        filter(.data[[g1]] == 1 | .data[[g2]] == 1) %>%
        transmute(SampleName, group_test = if_else(.data[[g1]] == 1, g1, g2), !!!extra_syms) %>%
        mutate(group_test = factor(group_test, levels = c(g1, g2)))
        
      exist_i <- exist %>%
        select(any_of(meta_i$SampleName)) %>%
        filter(rowSums(.) > 0 & rowSums(.) < ncol(.))
        
      # Timepoints (if applicable)
      if (!is.null(timepoints_file) && file.exists(timepoints_file)) {
        timepoints_df <- read.csv(timepoints_file, header = TRUE, check.names = FALSE)
        relevant_cols <- c("ind_id", g1, g2)
        available_cols <- intersect(relevant_cols, colnames(timepoints_df))
          
        timepoints_df <- timepoints_df %>%
          rename(ind_id = 1) %>%
          select(any_of(available_cols)) %>%
          filter(!is.na(.data[[g1]]) | !is.na(.data[[g2]])) %>%
          rowwise() %>%
          mutate(
            g1_exists = !is.na(.data[[g1]]) && .data[[g1]] %in% colnames(exist_i),
            g2_exists = !is.na(.data[[g2]]) && .data[[g2]] %in% colnames(exist_i)
          ) %>%
          ungroup() %>%
          filter(g1_exists | g2_exists)
      }
      
      render(
        input = template,
        output_file =  paste0(name, ".html"),
        output_dir  =  output_dir,
        params = list(
          metadata      = meta_i,
          exist         = exist_i,
          comparison    = name,
          library_meta  = lib_metadata_df,
          custom_colors = custom_colors,
          timepoints    = timepoints_df,
          out_tables    = out_tables
        ),
        envir = new.env()
      )

      time_taken <- toc(log = FALSE)
        
      # Write to log
      log_msg <- sprintf("[%s] ✅ Rendered %s in %.2f sec → %s", Sys.time(), name, time_taken$toc - time_taken$tic, output_dir)
      write(log_msg, file = log_file, append = TRUE)
        
      # Return a row for the summary
      tibble(
        comparison = name,
        group1 = g1,
        group2 = g2,
        output_file = output_dir,
        time_sec = round(time_taken$toc - time_taken$tic, 2),
        timestamp = Sys.time()
      )
        
    }, error = function(e) {
      err_msg <- sprintf("[%s] ❌ Error rendering %s: %s", Sys.time(), name, e$message)
      write(err_msg, file = log_file, append = TRUE)
        
      tibble(
        comparison = name,
        group1 = g1,
        group2 = g2,
        output_file = NA,
        time_sec = NA,
        timestamp = Sys.time()
      )
    })
  },
  .options = furrr::furrr_options(seed = TRUE)
)
  
# Write summary
write_csv(results, summary_file)
