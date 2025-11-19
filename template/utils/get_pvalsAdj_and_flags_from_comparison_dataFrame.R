library(dplyr)

comparison <- read.csv("../IBD-Chile/reports_agilent/Tables/table_peptidesSignificance_group_test_control_vs_disease.csv")


#' Creates a custom lineage column based on flexible rules.
#'
#' @param comparison_df The data frame containing the lineage columns (e.g., Class, Genus, Species).
#' @param default_lineage_col The lineage column to use as the default for all
#'   taxa not matched by the rules (e.g., "Species").
#' @param custom_rules A data frame or list of lists defining the custom groups.
#'   Each row/list must have:
#'   - 'lineage_col': The column name to check (e.g., "Class", "Genus").
#'   - 'taxa_name': The specific name to match (e.g., "Bacteroidales", "Streptococcus").
#' @param new_col_name The name for the resulting custom column (default: "Custom_Taxa").
#' @return The data frame with the new custom lineage column added.
create_flexible_taxa_column <- function(
    df,
    default_lineage_col,
    custom_rules,
    new_col_name = "custom_taxa"
) {
  # --- 1. VALIDATION ---
  
  # Check if the default column exists
  if (!(default_lineage_col %in% names(df))) {
    stop("Default lineage column '", default_lineage_col, "' not found in data frame.")
  }
  
  # Ensure custom_rules is a data frame for easier iteration
  if (!is.data.frame(custom_rules)) {
    custom_rules <- data.frame(custom_rules)
  }
  
  # Check if all columns specified in custom_rules exist
  required_cols <- unique(custom_rules$lineage_col)
  missing_cols <- required_cols[!(required_cols %in% names(df))]
  if (length(missing_cols) > 0) {
    stop("Lineage columns specified in rules are missing: ", paste(missing_cols, collapse = ", "))
  }
  
  # --- 2. BUILD MUTATE LOGIC (using case_when) ---
  
  # Initialize the case_when statement with the TRUE (default) condition
  # We use !!sym(default_lineage_col) to dynamically select the column
  case_conditions <- list(quo(TRUE ~ !!rlang::sym(default_lineage_col)))
  
  # Build the custom rules (Rule 1, Rule 2, etc. - in order)
  for (i in seq_len(nrow(custom_rules))) {
    rule <- custom_rules[i, ]
    
    lineage_sym <- sym(rule$lineage_col)
    taxa_name_val <- rule$taxa_name
    
    # Create the condition: e.g., Class == "Bacteroidales" ~ "Bacteroidales"
    condition <- quo((!!lineage_sym == !!taxa_name_val) ~ !!taxa_name_val)
    case_conditions <- c(list(condition), case_conditions) # Prepend to apply rules first
  }
  
  # --- 3. APPLY MUTATION ---
  
  # Use !!! to splice the list of quoted expressions into the case_when function
  df <- df %>%
    mutate(
      !!sym(new_col_name) := case_when(!!!case_conditions)
    )
  
  return(df)
}




get_top_significant_taxa_df <- function(comparison_df, lineage_col){#, n = 100) {
  # turn the string into a symbol once
  lineage_sym <- sym(lineage_col)
  
  top_taxa_df <- comparison_df %>%
    # drop rows where that lineage is missing
    filter(!is.na(!!lineage_sym)) %>%
    # make sure we have the same 'ratio' name
    mutate(ratio = log2(ratio)) %>%
    # get one row per lineage
    distinct(!!lineage_sym) %>%
    # for each lineage, pull out x = that group, y = all others
    mutate(
      x = purrr::map(!!lineage_sym, ~ comparison_df$ratio[comparison_df[[lineage_col]] == .x]),
      y = purrr::map(!!lineage_sym, ~ comparison_df$ratio[comparison_df[[lineage_col]] != .x]),
      p = purrr::map2_dbl(x, y, ~ wilcox.test(.x, .y)$p.value)
    ) %>%
    select(!!lineage_sym, p) %>%
    # FDR‐adjust
    adjust_pvalue(method = "BH") %>%
    arrange(p.adj) %>%
    # take the top n
    #slice_head(n = n) %>%
    # ---- CLEAN THE NAMES HERE ----
  # Modify the lineage column in place
  mutate(
    !!lineage_sym := str_remove_all(!!lineage_sym, "\\[|\\]"), # drop any brackets
    !!lineage_sym := str_squish(!!lineage_sym)               # collapse multiple spaces, trim
  ) %>%
    # Select the cleaned lineage and the adjusted p-value
    select(!!lineage_sym, p, p.adj)
  
  return(top_taxa_df)
}


make_flag_lists <- function(
    df,
    n = 10,
    taxa_labels = NULL, # New: Vector of taxa names to select (e.g., c("species1", "species2"))
    colors = NULL,
    brewer_palette = "Set3",
    brewer_n = 12
) {
  # Assuming the taxa column is the first one
  taxa_col_name <- names(df)[1]
  
  # Filter the data frame based on the provided labels
  if (!is.null(taxa_labels)){
    df_filtered <- df %>%
      filter(!!sym(taxa_col_name) %in% taxa_labels)
  } else {
    df_filtered <- head(df, n)
  }
  
  # Check if any of the labels were found
  top_names <- df_filtered[[taxa_col_name]]
  
  if (length(top_names) == 0) {
    stop("No specified taxa found in the first column of the input data frame.")
  }
  
  # derive short keys from first three words
  short_keys <- sapply(top_names, function(x) {
    words <- strsplit(x, "\\s+")[[1]]
    paste(head(words, 3), collapse = " ")
  })
  # ensure unique
  short_keys <- make.unique(short_keys)
  
  # build flags_to_patterns list: key = short_key, value = full name
  flags_to_patterns <- setNames(as.list(top_names), short_keys)
  n_taxa <- length(top_names)
  
  # build highlight_colors vector
  if (!is.null(colors)) {
    if (length(colors) != n_taxa) {
      stop("Length of 'colors' (", length(colors), ") must match number of flags (", n_taxa, ")")
    }
    # Use the provided colors, named by the new short keys
    highlight_colors <- setNames(colors, short_keys)
  } else {
    pal <- colorRampPalette(brewer.pal(brewer_n, brewer_palette))(n_taxa)
    highlight_colors <- setNames(pal, short_keys)
  }
  
  # return
  list(
    flags_to_patterns = flags_to_patterns,
    highlight_colors  = highlight_colors,
    p = df_filtered$p,
    p.adj = df_filtered$p.adj
  )
}



################################################################
#  Define the rules for custom grouping
my_custom_rules <- data.frame(
  lineage_col = c("order", "genus"),
  taxa_name = c("Bacteroidales", "Listeria"),
  stringsAsFactors = FALSE
)

# Create the new custom column
comparison_df_custom <- create_flexible_taxa_column(
  df = comparison,
  default_lineage_col = "species", # Default is species level
  custom_rules = my_custom_rules,
  new_col_name = "newTaxaCol"
)

#  Run your statistical analysis using the new column
results_df <- get_top_significant_taxa_df(
  comparison_df = comparison_df_custom,
  lineage_col = "newTaxaCol" # <--- Use the new custom column
)

make_flag_lists(
    results_df,
    n = 5,
    #c("Campylobacter coli", "Norwalk virus", "Bacillus cereus"), # New: Vector of taxa names to select (e.g., c("species1", "species2"))
    colors = NULL,
    brewer_palette = "Set3",
    brewer_n = 12
)
