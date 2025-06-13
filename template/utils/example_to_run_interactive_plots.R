#####################################
### Author: Carlos S. Reyna-Blanco###
###                               ###
###   Display interactive plots   ###
#####################################

source("../phipseq_report-generator/template/utils/make_interactive_plots.R")

# Prepare data -----

# Define the size and group 1 and group 2 name
N1 <- 55
N2 <- 55
group1 <- "mother_B"
group2 <- "BM_M3"

# groups and patterns to find that group
flags_to_patterns <- list(
  `Milk allergens` = c("twist_25139", "twist_43321", "twist_25139", "twist_5555", "twist_54532"),
  Enterovirus = c("Enterovirus"),
  Bacteriodes = c("Bacteroides") # different patterns here
  # can add more groups here
)

# Load table with significant peptides group
comparison_df <- read.csv(paste("reports_Jun13_2025/Tables/table_peptidesSignificance_group_test_", group1,"_vs_", group2,".csv", sep=""))

# Look for the patterns in the table and add new columns based on the intersting groups
for(flag in names(flags_to_patterns)){
  patterns <- flags_to_patterns[[flag]]
  # (A) Add a new TRUE/FALSE column "flag" by matching those patterns:
  comparison_df <- comparison_df %>%
    add_flag_by_patterns(
      new_flag    = flag,
      patterns    = patterns,
      target_cols = c("Organism_complete_name", "Description", "Peptide")
    )
}

flags <- c(names(flags_to_patterns)) # no extra annotation

#############################RUN ONLY IF INCLUDING EXTRA SUBGROUP ANNOTATION#####################################
# SUBGROUPS_TO_NAME <- c(
#   'all' = 'Complete library',
#   'is_PNP' = 'Metagen antigens',  'is_patho' = 'Pathogenic strains',
#   'is_probio' = 'Probiotic strains', 'is_MPA' = 'Microbiota strains', 'is_IgA' = 'Antibody-coated strains',
#   'is_bac_flagella' = 'Flagellins', 'is_infect' = 'Infectious pathogens',
#   'is_IEDB_or_cntrl' = 'IEDB/controls')
# actual_subgroup_columns <- setdiff(names(SUBGROUPS_TO_NAME), "all")
# 
# #load extra file with the annotation info and format it
# subgroup_lib_df <- readRDS("../phipseq_report-generator/library_meta/all_libraries_with_important_info.rds") %>%
#   tibble::rownames_to_column(var = "Peptide") %>%
#   mutate(
#     across(
#       all_of(actual_subgroup_columns),
#       ~ case_when(
#         is.na(.)                       ~ FALSE,
#         . %in% c(1, "1", TRUE, "True") ~ TRUE,
#         TRUE                           ~ FALSE
#       )
#     ),
#     all = TRUE
#   ) %>%
#   select(Peptide, all_of(names(SUBGROUPS_TO_NAME))) %>% 
#   rename(setNames(names(SUBGROUPS_TO_NAME), SUBGROUPS_TO_NAME))
# 
# # add extra col to the comparison_df
# comparison_df <- comparison_df %>% 
#   left_join(
#     subgroup_lib_df %>% 
#         select(Peptide, `Infectious pathogens`),  
#     by = "Peptide"
#     )
# 
# flags <- c("Infectious pathogens")
###################################################################################################################

# Call functions----

# Interactive mode
make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             highlight_cols   = flags, 
                             highlight_colors = c(
                               `Milk allergens` = "#1b9e77", Enterovirus = "#d95f02", 
                               Bacteriodes = "#7570b3",
                               #`Infectious pathogens` = "gold"
                               ),
                             default_color    = "gray70", 
                             interactive = T)


# No interactive mode
make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             highlight_cols   = flags, 
                             highlight_colors = c(`Milk allergens` = "#1b9e77", Enterovirus = "#d95f02", Bacteriodes = "#7570b3"),
                             default_color    = "gray70", 
                             interactive = F)

# Show significance
make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             default_color    = "gray70", 
                             interactive = T)
