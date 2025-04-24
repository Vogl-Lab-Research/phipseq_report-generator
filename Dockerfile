# Use an R base image with tidyverse
FROM rocker/tidyverse:latest

# Install additional R packages
RUN R -e "install.packages(c('rmarkdown', 'furrr', 'future', 'progressr', 'tictoc', 'yaml', 'tidyverse', 'rmarkdown', 'furrr', 'progressr', 'future', 'tictoc', 'readr', 'yaml',
  'stats', 'multcomp', 'dplyr', 'ggplot2', 'scales', 'ggsignif', 'plotly', 'nnet'))"

# Set working directory inside the container
WORKDIR /app

# Copy relevant files into the container
COPY make_reports_parallel_logs.R /app/
COPY make_reports.R /app/
COPY template/template_phipseq.Rmd /app/template
COPY library_meta/ /app/library_meta/
COPY utils/ /app/utils/

# Set default command
CMD ["Rscript", "make_reports_parallel_logs.R", "--help"]
