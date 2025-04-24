# PhIP-seq Report Generator


This tool generates automated HTML reports for PhIP-seq data, including key summary plots, diversity metrics, group comparisons, multidimensional scaling. Designed for high-throughput batch rendering with support for parallel processing.


# How to Run

From the command line, run:

Rscript path_to_report-generator/make_reports_parallel_logs.R config.yaml

# Configuration (config.yaml)

Define your input files and parameters in a config.yaml file. Example:

comparisons_file: Metadata/comparisons.csv 
samples_file: Metadata/sorted_LLNEXT_samples_binary.csv
exist_file: Data/exist.csv
timepoints_file: Metadata/LLNext_ind_timepoints.csv  # Optional for longitudinal studies
extra_cols: ["Sex", "Age"]  # Optional additional metadata columns

# Output
output_dir: reports  # Default directory for saving reports [relative or full path]

# Input Files format

Look at the mock_files folder for examples of simulated data as a guidance of format for the required files
