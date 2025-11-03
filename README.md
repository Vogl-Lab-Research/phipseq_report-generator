# PhIP-seq Report Generator

Phage ImmunoPrecipitation sequencing (PhIP-seq) enables high-throughput profiling of antibody binding
against large libraries of peptides. This repository provides an R-based toolkit that converts raw
PhIP-seq experiment outputs into shareable HTML summaries. The reports highlight key quality metrics,
group comparisons and sample diversity to streamline downstream interpretation.

## Features

- Automated HTML reports with summary plots, diversity metrics and multidimensional scaling
- Batch rendering with optional parallel execution for large datasets
- Configurable input files and custom metadata columns

## Installation

1. Clone this repository.
2. Ensure R is installed along with the required packages (`yaml`, `ggplot2`, `dplyr`, and `testthat`).
3. Adjust your configuration file as described below.

## Usage

From the command line, run:

```bash
Rscript path_to_report-generator/make_reports_parallel_logs.R config.yaml
```

The script reads information from the provided `config.yaml` and writes HTML reports to the configured
output directory.

## Configuration (`config.yaml`)

Define the input files and parameters in a YAML configuration file. Example:

```yaml
comparisons_file: Metadata/comparisons.csv
samples_file: Metadata/sorted_LLNEXT_samples_binary.csv
exist_file: Data/exist.csv
timepoints_file: Metadata/LLNext_ind_timepoints.csv  # Optional for longitudinal studies
extra_cols: ["Sex", "Age"]  # Optional additional metadata columns
prevalence_threshold: 0 # Optional, apply a prevalence filter where either one group must pass the defined threshold
# Output
output_dir: reports  # Directory for saving reports (relative or absolute path)
```

## Input file formats

The `mock_files` directory contains simulated data that serves as an example of the expected input
formats for each file referenced in the configuration. Use these examples as guidance when formatting
your own data files.

## Output

Each report is saved as a standalone HTML file in the specified `output_dir`. These reports can be
opened locally in a web browser or shared with collaborators.

