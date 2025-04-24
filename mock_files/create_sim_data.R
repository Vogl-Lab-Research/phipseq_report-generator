library(tidyverse)

# Set seed and parameters
set.seed(123)
n_samples <- 50
n_peptides <- 500
group_options <- c("mother_P12", "mother_P28", "mother_B")

# --- 1. Simulate longitudinal sample metadata ---

# Define number of individuals (shared across samples)
n_individuals <- 20
ind_ids <- paste0("ind_", 1:n_individuals)

# Randomly assign 1 to 3 groups per individual, then assign unique SampleName per group
sample_list <- list()
sample_counter <- 1

for (ind in ind_ids) {
  n_groups <- sample(1:3, 1)
  assigned_groups <- sample(group_options, n_groups)
  
  for (grp in assigned_groups) {
    sample_id <- paste0("sample_", sample_counter)
    
    sample_list[[sample_counter]] <- tibble(
      SampleName = sample_id,
      ind_id = ind,
      group = grp,
      Sex = sample(c(0, 1), 1),
      Age = sample(20:70, 1)
    )
    
    sample_counter <- sample_counter + 1
    if (sample_counter > n_samples) break
  }
  if (sample_counter > n_samples) break
}

# Combine sample data
samples_df <- bind_rows(sample_list)

# --- 2. Create binary group membership (samples_file) ---
samples_binary <- samples_df %>%
  pivot_wider(names_from = group, values_from = group, values_fn = length, values_fill = 0) %>%
  mutate(across(all_of(group_options), ~ ifelse(. > 0, 1, 0)))

samples_file <- samples_binary %>%
  select(SampleName, all_of(group_options), Sex, Age)

dir.create("Metadata", showWarnings = FALSE)
write_csv(samples_file, "unit_test/Metadata/samples_meta.csv")

# --- 3. Create comparisons file ---
comparisons <- tribble(
  ~comparison,                   ~group1,       ~group2,
  "mother_P12_vs_mother_P28",    "mother_P12",  "mother_P28",
  "mother_P28_vs_mother_B",      "mother_P28",  "mother_B",
  "mother_P12_vs_mother_B",      "mother_P12",  "mother_B"
)

write_csv(comparisons, "unit_test/Metadata/comparisons.csv")

# --- 4. Create exist.csv with 500 peptides x 50 samples ---
sample_names <- samples_file$SampleName
exist_matrix <- matrix(sample(0:1, n_peptides * length(sample_names), replace = TRUE),
                       nrow = n_peptides,
                       dimnames = list(paste0("pep_", 1:n_peptides), sample_names))

dir.create("Data", showWarnings = FALSE)
write.csv(exist_matrix, "unit_test/Data/exist.csv", quote = FALSE)


# --- 5. Timepoints file ---
timepoints_file <- samples_df %>%
  pivot_wider(
    id_cols = ind_id,
    names_from = group,
    values_from = SampleName,
    values_fn = list
  ) %>%
  mutate(across(
    all_of(group_options),
    ~ map_chr(., ~ if (length(.x) == 0) NA_character_ else .x[1])
  ))

write_csv(timepoints_file, "unit_test/Metadata/samples2ind_timepoints.csv")

message("✅ All simulated test files generated for 50 samples and 500 peptides.")
