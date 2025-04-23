#!/usr/bin/env Rscript
library(yaml)

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
library_meta: Metadata/all_libraries_with_important_info.rds
template_file: scripts/template.Rmd
timepoints_file: Metadata/LLNext_ind_timepoints.csv Optional
extra_cols: [Sex, Age] Optional
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
library_meta    <- config$library_meta
template        <- config$template_file
timepoints_file <- config$timepoints_file %||% NULL  # allow missing
extra_cols      <- config$extra_cols %||% character()  # allow missing
output_dir      <- config$output_dir %||% "reports"
is_absolute_path <- function(path) grepl("^(/|[A-Za-z]:)", path)  # Unix (/) or Windows (C:/)
output_dir <- if (is_absolute_path(output_dir)) {
  output_dir
} else {
  file.path(getwd(), output_dir)  # relative to where script is run from
}
if (!dir.exists(output_dir)) dir.create(output_dir)


# Validate required files
required <- list(cmp_file, samples_file, exist_file, library_meta, template)
names(required) <- c("comparisons_file", "samples_file", "exist_file", "library_meta", "template_file")

missing <- names(required)[!file.exists(unlist(required))]

if (length(missing)) {
  stop(paste("❌ Missing required files:", paste(missing, collapse = ", ")))
}

# Validate timepoints file if provided
if (!is.null(timepoints_file) && !file.exists(timepoints_file)) {
  stop("❌ Timepoints file does not exist: ", timepoints_file)
}

# Load libraries to parse, filter and render data
library(tidyverse)
library(rmarkdown)

# read inputs
comparisons <- read_csv(cmp_file)
samples     <- read_csv(samples_file) %>%
  rename_with(~ "SampleName", .cols = 1)
exist       <- read.csv(exist_file, header = TRUE, row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
lib_metadata_df <- readRDS(library_meta)
extra_syms <- syms(extra_cols)

# Get the directory of this script
args_full <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(normalizePath(sub("--file=", "", args_full[grep("--file=", args_full)])))
template <- file.path(script_path, "template_phipseq.Rmd")

# Set a palette of up to 12 distinct colors
available_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "bisque1", "#a6761d", "#666666",
  "#a6cee3", "#1f78b4", "#b2df8a", "salmon2")

# Color map to remember assignments
group_color_map <- list()


# Create reports
for(i in c(1)){#seq_len(nrow(comparisons))) { 
  comp <- comparisons[i, ]
  name <- pull(comp[1])
  g1   <- pull(comp[2])
  g2   <- pull(comp[3])


  
  groups <- c(g1, g2)
  for (g in groups) {
    if (!g %in% names(group_color_map)) {
      used_colors <- unlist(group_color_map)
      available <- setdiff(available_colors, used_colors)
      if (length(available) == 0) {
        warning("Ran out of distinct colors, reusing colors!!! Add more :)")
        available <- available_colors
      }
      group_color_map[[g]] <- available[1]
    }
  }
  custom_colors <- setNames(unlist(group_color_map[groups]), groups)
  
  # pull out just those samples that belong to g1 or g2
  meta_i <- samples %>%
    filter( .data[[g1]] == 1 | .data[[g2]] == 1 ) %>%
    transmute(
      SampleName,
      group_test = if_else(.data[[g1]] == 1, g1, g2),
      !!!extra_syms
    ) %>%
    mutate(
      # force the factor levels so g1 is first, g2 second
      group_test = factor(group_test, levels = c(g1, g2))
    )

  # filter your existence matrix to just those SampleName
  exist_i <- exist %>%
    select(any_of(meta_i$SampleName)) %>% 
    dplyr::filter(rowSums(.) > 0 & rowSums(.) < ncol(.))

  if (!is.null(timepoints_file) && file.exists(timepoints_file)) {
    relevant_cols <- c("ind_id", g1, g2)
    timepoints_df <- read.csv(timepoints_file, header = TRUE, check.names = FALSE)
    available_cols <- intersect(relevant_cols, colnames(timepoints_df))
    timepoints_df <- timepoints_df %>% 
      rename(ind_id = 1) %>% 
      dplyr::select(any_of(available_cols)) %>%
      filter(!is.na(.data[[g1]]) | !is.na(.data[[g2]])) %>%
      rowwise() %>%
      mutate(
        g1_exists = !is.na(.data[[g1]]) && .data[[g1]] %in% colnames(exist_i),
        g2_exists = !is.na(.data[[g2]]) && .data[[g2]] %in% colnames(exist_i)
      ) %>%
      ungroup() %>%
      filter(g1_exists | g2_exists)
  }
  
  # render directly, passing the objects
  render(
    input          = template,
    output_file =  paste0(name, ".html"),
    output_dir  =  output_dir,
    params = list(
      metadata     = meta_i,
      exist        = exist_i,
      comparison   = name,
      library_meta = lib_metadata_df,
      custom_colors = custom_colors,
      timepoints    = timepoints_df  # only if longitudinal!
    ),
    envir = new.env()
  )
}