library(tidyverse)

setwd("/home/creyna/Vogl-lab_Projects_git/MCI-Dementia/Tables")

# <<< set your folder here >>>
in_dir <- "/home/creyna/Vogl-lab_Projects_git/MCI-Dementia/Tables"

# load rds files to get description, lineages and protein and aa information and sequences
prot_seq <- readRDS("/home/creyna/Vogl-lab_Projects_git/Annotations/protein-peptide_seq_agilent-twist.rds") %>% 
  select(Peptide, full_aa_seq)
peptide_meta <- read.csv("/home/creyna/Vogl-lab_Projects_git/Annotations/combined_libraries_with_lineages_important_info.csv", check.names = F) %>%
  rename(Peptide = 1) %>%  # rename first column to "Peptide"
  select(Peptide, Description, species, len_seq, pos,
         common, genus, family, order, class, phylum, kingdom, domain,aa_seq)  # select only relevant columns


# -- merge tables ---- 
# find all relevant tables
files <- list.files(
  in_dir,
  pattern = "*\\.csv$",
  full.names = TRUE
)

tables <- purrr::map(files, function(f) {
  # extract label between "group_test_" and ".csv"
  label <- stringr::str_match(basename(f), "(.+)\\.csv$")[, 2]
  
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # identify the "_count" columns
  count_cols <- grep("_count$", names(df), value = TRUE)
  # group columns are the names without "_count"
  group_cols <- gsub("_count$", "", count_cols)
  
  # select: Peptide, ratio, pval, and those group columns
  out <- df %>%
    dplyr::select(Peptide, ratio, pvals_not_adj, pvals_bh, dplyr::all_of(group_cols)) %>%
    dplyr::rename(
      !!paste0(label, "_ratio")    := ratio,
      !!paste0(label, "_raw_pval") := pvals_not_adj,
      !!paste0(label, "_adj_pval") := pvals_bh
    ) %>%
    dplyr::distinct(Peptide, .keep_all = TRUE) %>%
    # ensure numeric
    dplyr::mutate(dplyr::across(-Peptide, as.numeric))
  
  out
})

# full-join across all files by Peptide
merged_df <- purrr::reduce(tables, dplyr::full_join, by = "Peptide") %>%
  dplyr::arrange(Peptide)

# ---- remove duplicated group columns (keep first) ----
merged_df <- merged_df %>%
  dplyr::select(!dplyr::any_of(names(.)[duplicated(names(.))]))

# ---- reorder columns: Peptide, all ratios, then all group cols, then all pvals ----
ratio_cols <- grep("_ratio$", names(merged_df), value = TRUE)
raw_pval_cols  <- grep("_raw_pval$", names(merged_df), value = TRUE)
pval_cols  <- grep("_adj_pval$", names(merged_df), value = TRUE)

# group columns = anything that came from the "_count" derivation
# i.e. exclude Peptide, ratios, and pvals
group_cols <- setdiff(
  names(merged_df),
  c("Peptide", ratio_cols, raw_pval_cols, pval_cols)
)

merged_df <- merged_df %>%
  dplyr::select(Peptide,  all_of(group_cols), all_of(ratio_cols), all_of(raw_pval_cols), all_of(pval_cols))

# ---- sort rows by all pval columns (smallest first) ----
merged_df <- merged_df %>%
  dplyr::arrange(across(all_of(pval_cols)))


# inspect / save
#print(merged_df, n = 20, width = Inf)


# --- get final merged table with all important info ------
merged_df <- merged_df %>%
  left_join(peptide_meta, by = "Peptide") %>% 
  relocate(Description, species, len_seq, pos, .after = Peptide) %>% 
  left_join(prot_seq, by = "Peptide")
write.csv(merged_df, file.path(in_dir, "EBVpeptides_merged_ratios_pvals.csv"),quote = F, row.names = F)


# -------- Define group to look at and name if the table -------------------
# wheat_proteins_df  <- merged_df %>%  filter(class %in% "Magnoliopsida") %>%  arrange(VG_B_vs_CS_B_adj_pval)
# write.csv(wheat_proteins_df, file.path(in_dir, "wheat_proteins_blablabla.csv"),quote = F, row.names = F)

#merged_df %>%  filter(species %in% "Staphylococcus aureus"))
#merged_df %>%  filter(species %in% "Lymphocryptovirus humangamma4"))
#merged_df %>%  filter(is_flagellum %in% "True"))