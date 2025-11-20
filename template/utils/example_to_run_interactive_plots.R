#####################################
### Author: Carlos S. Reyna-Blanco###
###                               ###
###   Display interactive plots   ###
#####################################

source("../phipseq_report-generator/template/utils/make_interactive_plots.R")

# Define the size and group 1 and group 2 name
N1 <- 63
N2 <- 156
group1 <- "control_healty"
group2 <- "SNMG_no_col"

# Prepare data -----
# Load table with significant peptides group
comparison_df <- read.csv(paste("../SNMG//reports//Tables/table_peptidesSignificance_group_test_", group1,"_vs_", group2,".csv", sep=""))


#########################################
#  Only run if you want to define the rules for custom grouping if you want to include taxa from different lineages in one column
########################################
my_custom_rules <- data.frame(
  lineage_col = c("genus", "genus"),
  taxa_name = c("Bacteroides", "Enterovirus"),
  stringsAsFactors = FALSE)
# Create the new custom column that includes taxa from different lineages
comparison_df_custom <- create_flexible_taxa_column(
  df = comparison_df,
  default_lineage_col = "species", # Default is species level
  custom_rules = my_custom_rules,
  new_col_name = "newTaxaCol")
#######################################

#  Run statistical analysis using a lineage column (can be species, order or the new generated taxa col)
pvals_df <- get_top_significant_taxa_df(
  comparison_df = comparison_df_custom,
  lineage_col = "newTaxaCol" # <--- here provide lineage col (e.g. species, order, genus or a new custom taxa column)
)


#############
# groups and patterns to find that group
flags_to_patterns <- list(
  Bos = c("Bos taurus"),
  #`Milk allergens` = c("twist_25139", "twist_43321", "twist_25139", "twist_5555", "twist_54532"),
  CMV = c("Cytomegalovirus humanbeta5"),
  Enterovirus = c("Enterovirus"),
  Bacteroides = c("Bacteroides") # different patterns here
  # can add more groups here
)
flags <- c(names(flags_to_patterns)) # no extra annotation
taxa_to_plot <- unlist(flags_to_patterns)
pvals_taxa <- pvals_df[pvals_df[,1] %in% taxa_to_plot,]$p.adj
col_taxa <- c(Bos = "purple", CMV = "dodgerblue3", Enterovirus="#d95f02", Bacteroides="forestgreen")

# flags_list <- make_flag_lists(pvals_df,
#                               #n = 5,
#                               taxa_labels =  taxa_to_plot, # New: Vector of taxa names to select (e.g., c("species1", "species2"))
#                               colors = col_taxa)


# Look for the patterns in the table and add new columns based on the intersting groups
for(flag in flags){
  patterns <- flags_to_patterns[[flag]]
  #patterns <- significant_taxa_list[[rank_level]][[comp_name]]$flags_to_patterns[[flag]]
  # (A) Add a new TRUE/FALSE column "flag" by matching those patterns:
  comparison_df <- comparison_df %>%
    add_exact_flag( #can use add_flag_by_patterns to be more flexible in the search or if you want to use Description col and look up a word
      new_flag    = flag,
      patterns    = patterns,
      target_cols = c("species", "genus", "order", "class","Description") #can add Peptide if necessary
    )
}

#############################RUN ONLY IF INCLUDING Functional SUBGROUP ANNOTATION#####################################
# subgroup_lib_df <-  readRDS("../phipseq_report-generator/library_meta/combined_libraries_with_lineages_important_info_nonAAseq.rds") %>%
#   tibble::rownames_to_column(var = "Peptide") %>%
#   mutate(
#     across(
#       all_of(SUBGROUPS_TO_INCLUDE),
#       ~ if_else(is.na(.x), FALSE, .x)
#     )
#   ) %>%
#   select(Peptide, all_of(SUBGROUPS_TO_INCLUDE))

#############################RUN ONLY IF INCLUDING EXTRA SUBGROUP ANNOTATION#####################################
# SUBGROUPS_TO_NAME <- c(
#   #'all' = 'Complete library',
#   'is_PNP' = 'Metagen antigens',  'is_patho' = 'Pathogenic strains',
#   'is_probio' = 'Probiotic strains', 'is_MPA' = 'Microbiota strains', 'is_IgA' = 'Antibody-coated strains',
#   'is_bac_flagella' = 'Flagellins', 'is_infect' = 'Infectious pathogens',
#   'is_IEDB_or_cntrl' = 'IEDB/controls')
# #actual_subgroup_columns <- setdiff(names(SUBGROUPS_TO_NAME), "all")
# 
# #load extra file with the annotation info and format it
# subgroup_lib_df <- readRDS("../phipseq_report-generator/library_meta/all_libraries_with_important_info.rds") %>%
#   tibble::rownames_to_column(var = "Peptide") %>%
#   mutate(
#     across(
#       all_of(SUBGROUPS_TO_NAME),
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
source("../phipseq_report-generator/template/utils/make_interactive_plots.R")
make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             highlight_cols   = flags, 
                             highlight_colors = col_taxa,
                             pvals_adj = pvals_taxa,
                             # highlight_colors = c(
                             #   CMV = "dodgerblue3",
                             #   Enterovirus = "#d95f02", 
                             #   Bacteriodes = "#7570b3"
                             #   #`Infectious pathogens` = "gold"
                             #   ),
                             default_color    = "gray70", 
                             interactive = T)


make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             highlight_cols   = (flags), 
                             highlight_colors = (col_taxa),
                             reverse_legend = F,
                             pvals_adj = (pvals_taxa),
                             default_color    = "gray70", 
                             interactive = F)

make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             highlight_cols   = (flags), 
                             highlight_colors = (col_taxa),
                             reverse_legend = F,
                             pvals_adj = (pvals_taxa),
                             default_color    = "gray70", 
                             interactive = T)

# Show significance
make_interactive_scatterplot(comparison_df = comparison_df,
                             group1 = group1, group2 = group2, N = c(N1,N2),
                             default_color    = "gray70", 
                             interactive = T)
