# Use an R base image with tidyverse
FROM rocker/tidyverse:latest

# Install additional R packages
RUN R -e "install.packages(c('rmarkdown', 'furrr', 'future', 'progressr', 'tictoc', 'yaml'))"

# Set working directory inside the container
WORKDIR /app

# Copy relevant files into the container
COPY make_reports_parallel_logs.R /app/
COPY make_reports.R /app/
COPY template_phipseq.Rmd /app/
COPY utils/ /app/utils/

# Set default command
CMD ["Rscript", "make_reports_parallel_logs.R", "--help"]
