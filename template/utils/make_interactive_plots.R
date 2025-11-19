#####################################
### Author: Carlos S. Reyna-Blanco###
###                               ###
###   Display interactive plots   ###
#####################################

library("ggplot2")
library("tidyr")
library("dplyr")
library("stringr")
library("plotly")

# Functions ----

# ---- Format p-values nicely ----
format_pval <- function(p, alpha = 0.05) {
  
  # helper to drop trailing zeros (e.g. "1.00" → "1", "0.500" → "0.5")
  drop_zeros <- function(x) sub("\\.?0+$", "", x)
  
  if (is.na(p)) return("NA")
  
  # -------------------------------------------------------------------------
  # 1. Non-significant (p > alpha)
  # -------------------------------------------------------------------------
  if (p > alpha) {
    raw <- formatC(p, digits = 2, format = "f")
    raw <- drop_zeros(raw)
    return(paste0("ns [", raw, "]"))
  }
  
  # -------------------------------------------------------------------------
  # 2. Normal fixed-decimal formatting (0.001 ≤ p ≤ alpha)
  # -------------------------------------------------------------------------
  if (p >= 0.001) {
    raw <- formatC(p, digits = 3, format = "f")
    raw <- drop_zeros(raw)
    return(raw)
  }
  
  # -------------------------------------------------------------------------
  # 3. Scientific notation (< 0.001)
  # -------------------------------------------------------------------------
  raw <- formatC(p, digits = 2, format = "e")
  
  # remove unnecessary zeros: "1.00e-05" → "1e-05"
  raw <- sub("([0-9]+)\\.0+e", "\\1e", raw)
  
  # remove trailing zeros inside "1.10e-04" → "1.1e-04"
  raw <- sub("([0-9]+\\.[0-9]*[1-9])0+e", "\\1e", raw)
  
  return(raw)
}



add_flag_by_patterns <- function(df,
                                 new_flag,
                                 patterns,
                                 target_cols = c("species", "Description")) {
  # build our two regexes:
  regex_sub   <- paste0("(?i)", paste(patterns, collapse = "|"))
  regex_exact <- paste0("(?i)^(", paste(patterns, collapse = "|"), ")$")
  
  # split out any Peptide column
  pep_present  <- "Peptide" %in% target_cols
  nonpep_cols  <- setdiff(target_cols, "Peptide")
  
  df %>%
    mutate(
      # temporary flags
      sub_flag = if (length(nonpep_cols) > 0) {
        if_any(all_of(nonpep_cols), ~ str_detect(.x, regex(regex_sub)))
      } else {
        # no non-Peptide columns requested
        FALSE
      },
      pep_flag = if (pep_present) {
        # exact match only on the Peptide column
        str_detect(.data$Peptide, regex(regex_exact))
      } else {
        FALSE
      },
      # final flag is TRUE if either test passes
      !!new_flag := sub_flag | pep_flag
    ) %>%
    # clean up
    select(-sub_flag, -pep_flag)
}

add_exact_flag <- function(df,
                           new_flag,
                           patterns,
                           target_cols = c("species", "Description")) {
  # build an anchored regex: ^(A|B|C)$
  regex_exact <- paste0("^(", paste(patterns, collapse="|"), ")$")
  
  df %>%
    mutate(
      # across each target col, test for an exact match
      !!new_flag := if_any(all_of(target_cols),
                           ~ stringr::str_detect(.x, regex_exact))
    )
}


make_interactive_scatterplot <- function(comparison_df,
                                         group1, group2, N,
                                         highlight_cols   = NULL,
                                         highlight_colors = NULL,
                                         pvals_adj = NULL,
                                         default_color    = "gray70",
                                         #multiple_color   = "black",
                                         significant_colors = c(
                                           "not significant"                 = "dodgerblue",
                                           "significant prior correction"    = "forestgreen",
                                           "significant post FDR correction" = "firebrick"),
                                         interactive = TRUE) {
  # sanity-check:
  if (!is.null(highlight_cols)) {
    missing_cols <- setdiff(highlight_cols, names(comparison_df))
    if (length(missing_cols)) {
      stop("These highlight_cols are not in your data frame: ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  # build the collapsed factor ------------------------------------------------
  if (!is.null(highlight_cols) && length(highlight_cols) > 0) {
    comparison_df <- comparison_df %>%
      rowwise() %>%
      mutate(
        # collect the names of all TRUE flags in this row:
        .trues = list(highlight_cols[ c_across(all_of(highlight_cols)) ]),
        # now assign highlight:
        highlight = if (length(.trues) == 0) {
          "none"
        } else if (length(.trues) >= 1) {
          .trues[[1]]
        } #else {
        #"multiple"}
      ) %>%
      ungroup() %>%
      select(-.trues)
    
    comparison_df <- comparison_df %>% 
      mutate(highlight = replace_na(highlight, "none"))
    # ensure factor has all levels:
    levels_needed <- c("none",  highlight_cols) #multiple
    comparison_df$highlight <- factor(
      comparison_df$highlight,
      levels = levels_needed
    )
    # build tooltip (only for highlighted points)
    comparison_df <- comparison_df %>%
      mutate(
        log2ratio = log2( ratio ),
        tooltip_txt = if_else(
          highlight == "none",
          NA_character_,
          paste0(
            "Peptide: ",  Peptide,               "<br>",
            "Desc: ",     Description,           "<br>",
            "Species: ",  species,"<br>",
            group1, ": ", !!sym(group1), " / ",
            group2, ": ", !!sym(group2),       "<br>",
            "Highlight: ", highlight
          )
        )
      ) %>%
      filter(
        is.finite(.data[[group1]]),
        is.finite(.data[[group2]])
      ) %>%
      arrange(highlight)
    
    
    if (is.null(pvals_adj)){
      pvals <- sapply(highlight_cols, function(flag) {
        x <- comparison_df %>% filter( !!sym(flag) ) %>% pull(log2ratio)
        #y <- comparison_df$log2ratio
        # better subset y to only non-flag too:
        y <- comparison_df %>% filter(! (!!sym(flag)) ) %>% pull(log2ratio)
        w <- wilcox.test(x, y)
        w$p.value
      })
      pvals_adj <- p.adjust(pvals, method = "BH")
    }
    fmt_p <- sapply(pvals_adj, format_pval)
    legend_labels <- paste0(highlight_cols, " (P=", fmt_p, ")")
    
    # colors: user‐supplied or a simple default palette
    if (is.null(highlight_colors)) {
      # pick a palette for the flags
      palette_vals <- setNames(
        RColorBrewer::brewer.pal(
          n = max(length(highlight_cols), 8),
          name = "Set2"
        )[1:length(highlight_cols)],
        highlight_cols
      )
    } else {
      palette_vals <- highlight_colors
    }
    manual_vals <- c(
      none     = default_color,
      #multiple = multiple_color,
      palette_vals
    )
    
    color_aes   <- aes(color = highlight, text = tooltip_txt)
    color_scale <- scale_color_manual(
      name   = NULL,
      values = manual_vals,
      #limits = levels_needed,    # ← ensures “none” is the first group drawn
      breaks = highlight_cols,
      labels = legend_labels
      #labels = c("Milk allergens", "Enterovirus", "Bacteriodes")
    )
    legend_theme <- theme(
      #legend.position   = "top",            # place above the plot
      #legend.justification = "center",      # center it
      legend.position     = c(0, 1),   # 50% across, 95% up
      legend.justification = c(0, 1),  
      #legend.direction  = "horizontal",     # lay keys out side-by
      legend.background    = element_rect(fill = alpha("white", 0.8), color = "gray80"),
      legend.key.size      = unit(10, "pt"),
      legend.text          = element_text(size = 9),
      legend.title         = element_text(size = 9, face = "bold")
    )
    show_legend <- TRUE
    names(legend_labels) <- highlight_cols
    
  } else {
    # no highlights requested → fall back to your old categories logic
    comparison_df <- comparison_df %>%
      mutate(
        tooltip_txt = ifelse(
          categories == "not significant",
          #categories == "ns", #"not significant",
          
          NA_character_,
          paste0(
            "Peptide: ",  Peptide,               "<br>",
            "Desc: ",     Description,           "<br>",
            "Species: ",  species,"<br>",
            group1, ": ", !!sym(group1), " / ",
            group2, ": ", !!sym(group2)
          )
        )
      ) %>%
      filter(
        is.finite(.data[[group1]]),
        is.finite(.data[[group2]])
      )
    
    color_aes   <- aes(color = categories, text = tooltip_txt)
    color_scale <- scale_color_manual(
      values = significant_colors, 
      labels = c("ns", "significant", "significant FDR"),
      name = NULL)
    legend_theme <- theme(legend.position = "none")
    show_legend  <- FALSE
  }
  
  # build the ggplot + ggplotly -----------------------------------------------
  p <- ggplot(comparison_df,
              aes(x = !!sym(group1), y = !!sym(group2))) +
    geom_point(color_aes, alpha = 0.65) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray4") +
    color_scale +
    labs(
      x = paste0("% ", group1, " in whom a peptide is\nsignificantly bound (n = ", N[1], ")"),
      y = paste0("% ", group2, " in whom a peptide is\nsignificantly bound (n = ", N[2], ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA),
      plot.margin      = margin(t = 10, r = 15, b = 15, l = 10, unit = "pt"),
      axis.text.y.left = element_text(size = 10, face = "italic"),
      axis.title.x     = element_text(face = "bold"),
      axis.title.y     = element_text(face = "bold")
    ) +
    legend_theme
  
  if (interactive){
    interactive_plot <- ggplotly(p, tooltip = "text",
                                 width   = 550, height  = 550)
    
    
    if (!is.null(highlight_cols) && length(highlight_cols) > 0) {
      # only then do the trace‐name patching:
      for (i in seq_along(interactive_plot$x$data)) {
        tr    <- interactive_plot$x$data[[i]]
        nm    <- tr$name
        # hide the greys
        if (nm %in% c("none","multiple")) {
          tr$showlegend <- FALSE
        }
        # relabel the real flags
        else if (nm %in% highlight_cols) {
          tr$name <- legend_labels[[nm]]
        }
        interactive_plot$x$data[[i]] <- tr
      }
    }
    
    return(
      interactive_plot %>% layout(
        showlegend = show_legend,
        legend = list(
          #orientation = "h",       # horizontal keys
          x       = 0,       # center
          xanchor = "left",
          y       = 1,      # 95% up the plot area
          yanchor = "top",
          font        = list(size = 9)
        ),
        margin     = list(l = 80, r = 80, b = 80, t = 80, pad = 0),
        hoverlabel = list(font = list(size = 10)),
        # xaxis      = list(scaleratio = 1, scaleanchor = "y"),
        # yaxis      = list(scaleratio = 1, scaleanchor = "x")
        # -----------------------------------------------------------------
        # --- CHANGES TO ADD THE PLOT FRAME (AXIS LINES) ---
        # -----------------------------------------------------------------
        xaxis      = list(
          scaleratio = 1, 
          scaleanchor = "y",
          showline = TRUE,         #  Add the axis line
          linewidth = 1,           #  Set line thickness (optional, default is fine)
          linecolor = 'black',     #  Set line color
          mirror = TRUE            #  Mirror the line on the top
        ),
        yaxis      = list(
          scaleratio = 1, 
          scaleanchor = "x",
          showline = TRUE,         #  Add the axis line
          linewidth = 1,           #  Set line thickness
          linecolor = 'black',     #  Set line color
          mirror = TRUE            #  Mirror the line on the right
        )
      )
    )
  } else {
    return(p)
  }
}





#########################################
############ Volcano plot ###############
#########################################
make_interactive_volcano <- function(comparison_df, group1, group2,
                                     fc_cut = 1,
                                     p_cut  = 0.05,
                                     highlight_cols   = NULL,
                                     highlight_colors = NULL,
                                     pvals_adj = NULL,
                                     default_color    = "gray70",
                                     significant_colors = c(
                                       "not significant"                 = "dodgerblue",
                                       "significant prior correction"    = "forestgreen",
                                       "significant post FDR correction" = "firebrick"),
                                     interactive = TRUE){
  
  # sanity-check:
  if (!is.null(highlight_cols)) {
    missing_cols <- setdiff(highlight_cols, names(comparison_df))
    if (length(missing_cols)) {
      stop("These highlight_cols are not in your data frame: ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  # build the collapsed factor ------------------------------------------------
  if (!is.null(highlight_cols) && length(highlight_cols) > 0) {
    plotly_width <- 650
    margins <- list(l   = 80, r   = 120, t   = 80, b   = 40, pad = 0)
    comparison_df <- comparison_df %>%
      rowwise() %>%
      mutate(
        # collect the names of all TRUE flags in this row:
        .trues = list(highlight_cols[ c_across(all_of(highlight_cols)) ]),
        # now assign highlight:
        highlight = if (length(.trues) == 0) {
          "none"
        } else if (length(.trues) >= 1) {
          .trues[[1]]
        } #else {
        #"multiple"}
      ) %>%
      ungroup() %>%
      select(-.trues) 
    
    # ensure factor has all levels:
    comparison_df <- comparison_df %>% 
      mutate(highlight = replace_na(highlight, "none"))
    
    levels_needed <- c("none", highlight_cols) #multiple
    comparison_df$highlight <- factor(
      comparison_df$highlight,
      levels = levels_needed
    )
    # build tooltip (only for highlighted points)
    comparison_df <- comparison_df %>%
      mutate(
        log2ratio = log2( ratio ),
        tooltip_txt = if_else(
          highlight == "none",
          NA_character_,
          paste0(
            "Peptide: ",  Peptide,               "<br>",
            "Desc: ",     Description,           "<br>",
            "Species: ",  species,"<br>",
            group1, ": ", !!sym(group1), " / ",
            group2, ": ", !!sym(group2),       "<br>",
            "Highlight: ", highlight
          )
        )
      ) %>%
      filter(is.finite(log2ratio), is.finite(-log10(pvals_not_adj))) %>%  # drop any Inf or NaN 
      arrange(highlight)
    
    if (is.null(pvals_adj)){
      pvals <- sapply(highlight_cols, function(flag) {
        x <- comparison_df %>% filter( !!sym(flag) ) %>% pull(log2ratio)
        #y <- comparison_df$log2ratio
        # better subset y to only non-flag too:
        y <- comparison_df %>% filter(! (!!sym(flag)) ) %>% pull(log2ratio)
        w <- wilcox.test(x, y)
        w$p.value
      })
      pvals_adj <- p.adjust(pvals, method = "BH")
    }
    fmt_p <- sapply(pvals_adj, format_pval)
    legend_labels <- paste0(highlight_cols, " (P=", fmt_p, ")")
    
    
    # colors: user‐supplied or a simple default palette
    if (is.null(highlight_colors)) {
      # pick a palette for the flags
      palette_vals <- setNames(
        RColorBrewer::brewer.pal(
          n = max(length(highlight_cols), 8),
          name = "Set2"
        )[1:length(highlight_cols)],
        highlight_cols
      )
    } else {
      palette_vals <- highlight_colors
    }
    manual_vals <- c(
      none     = default_color,
      #multiple = multiple_color,
      palette_vals
    )
    
    color_aes   <- aes(color = highlight, text = tooltip_txt)
    color_scale <- scale_color_manual(
      name   = NULL,
      values = manual_vals,
      #limits = levels_needed,    # ← ensures “none” is the first group drawn
      breaks = highlight_cols,
      labels = legend_labels
      #labels = c("Milk allergens", "Enterovirus", "Bacteriodes")
    )
    legend_theme <- theme(
      #legend.position   = "top",            # place above the plot
      #legend.justification = "center",      # center it
      legend.position     = c(0, 1),   # 50% across, 95% up
      legend.justification = c(0, 1),  
      #legend.direction  = "horizontal",     # lay keys out side-by
      legend.background    = element_rect(fill = alpha("white", 0.8), color = "gray80"),
      legend.key.size      = unit(10, "pt"),
      legend.text          = element_text(size = 8),
      legend.title         = element_text(size = 9, face = "bold")
    )
    show_legend <- TRUE
    names(legend_labels) <- highlight_cols
    
  } else {
    # no highlights requested → fall back to your old categories logic
    plotly_width <- 500
    margins <- list(l = 0, r = 0, t = 0, b = 0, pad = 2)
    comparison_df <- comparison_df %>%
      mutate(
        tooltip_txt = ifelse(
          #categories == "ns", #"not significant",
          categories == "not significant",
          NA_character_,
          paste0(
            "Peptide: ",  Peptide,               "<br>",
            "Desc: ",     Description,           "<br>",
            "Species: ",  species,"<br>",
            group1, ": ", !!sym(group1), " / ",
            group2, ": ", !!sym(group2)
          )
        )
      ) %>%
      filter(
        is.finite(.data[[group1]]),
        is.finite(.data[[group2]])
      )
    
    color_aes   <- aes(color = categories, text = tooltip_txt)
    color_scale <- scale_color_manual(
      values = significant_colors, 
      labels = c("ns", "significant", "significant FDR"),
      name = NULL)
    legend_theme <- theme(legend.position = "none")
    show_legend  <- FALSE
  }
  
  
  p <- ggplot(comparison_df, aes(x = -log2(ratio), y=-log10(pvals_not_adj))) +
    
    # #geom_point(data = subset(comparison_df, passed_not_adj == "Yes" & log2(ratio) > 0), 
    # #                   aes(color = "significant prior correction group 1"), alpha = 0.6) +
    # #geom_point(data = subset(comparison_df, passed_not_adj == "Yes" & log2(ratio) < 0), 
    # #                   aes(color = "significant prior correction group 2"), alpha = 0.6) +
    
    geom_point(color_aes, alpha = 0.65) +
    geom_hline( yintercept = -log10(p_cut), linetype   = "dashed", color = "gray50") +
    geom_vline( xintercept = c(fc_cut,-fc_cut), linetype   = "dashed", color = "gray50") +
    geom_vline( xintercept = 0, linetype   = "dashed", color = "gray50",alpha = 0.5) +
    # scale_color_manual(
    #   values = volcano_colors, 
    #   name   = NULL) +
    color_scale +
    labs(
      x = paste0("log₂-ratio of antibody responses\nin ", group1, " and ", group2),
      y = paste0("-log₁₀(p-value)")
    ) +
    theme_bw(base_size = 12) +  # Use a minimal theme for elegance and set a base font size
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(colour = "black", fill = NA),  # Keep border if desired
      plot.margin = margin(t = 10, r = 15, b = 15, l = 10, unit = "pt"),  # Adjust margins in point
      axis.text.y.left = element_text(size = 10, face = "italic"),  # Style y-axis labels
      axis.title.x = element_text(face = "bold"),  # Add margin and bold to x-axis title
      axis.title.y = element_text(face = "bold")  # Bold y-axis title
    ) +
    legend_theme
  
  if (interactive){
    interactive_plot <- ggplotly(p, tooltip = c("text", "x", "y"),
                                 width   = plotly_width, height  = 500)
    
    
    if (!is.null(highlight_cols) && length(highlight_cols) > 0) {
      # only then do the trace‐name patching:
      for (i in seq_along(interactive_plot$x$data)) {
        tr    <- interactive_plot$x$data[[i]]
        nm    <- tr$name
        # hide the greys
        if (nm %in% c("none","multiple")) {
          tr$showlegend <- FALSE
        }
        # relabel the real flags
        else if (nm %in% highlight_cols) {
          tr$name <- legend_labels[[nm]]
        }
        interactive_plot$x$data[[i]] <- tr
      }
    }
    
    return(
      interactive_plot %>% layout(
        showlegend = show_legend,
        legend = list(
          #orientation = "h",       # horizontal keys
          x       = 1.02,       # center
          xanchor = "left",
          y       = 1,      # 95% up the plot area
          yanchor = "top",
          # bgcolor     = "rgba(255,255,255,0.8)",
          # bordercolor = "gray80",
          # borderwidth = 0,
          font        = list(size = 9)
        ),
        #margin      = list(l = 0, r = 0, t = 0, b = 0, pad = 2),
        margin = margins, #list(l   = 80, r   = 120, t   = 80, b   = 40, pad = 0),
        hoverlabel  = list(font = list(size = 10)),
        # xaxis       = list(automargin = TRUE),
        # yaxis       = list(automargin = TRUE,
        #                    title      = list(standoff = 10)),
        
        xaxis      = list(
          automargin = TRUE,
          showline = TRUE,         #  Add the axis line
          linewidth = 1,           #  Set line thickness (optional, default is fine)
          linecolor = 'black',     #  Set line color
          mirror = TRUE            #  Mirror the line on the top
        ),
        yaxis      = list(
          automargin = TRUE,
          title      = list(standoff = 10),
          showline = TRUE,         #  Add the axis line
          linewidth = 1,           #  Set line thickness
          linecolor = 'black',     #  Set line color
          mirror = TRUE            #  Mirror the line on the right
        ),
        
        
        autosize    = FALSE
      )
    )
  } else {
    return(p)
  }
}



###############get taxa_vs_rest with pvals adjusted


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



# Updated plot_ratios_by_subgroup with simplified subgroup logic
# • which_subgroups = "all" (subgroups + any new organism), "default (only subgroups", or "added (only new organisms"
plot_ratios_by_subgroup <- function(
    comparison_df,
    group1,
    group2,
    subgroup_lib_df,
    subgroup_colors      = NULL,
    pvals_adj = NULL,
    add_subgroups        = NULL,
    prevalence_threshold = 0,
    which_subgroups      = c("all", "default", "added"),
    x_label = "Subgroups of the antigen library"
) {
  which_subgroups <- match.arg(which_subgroups)
  
  
  # error if no additional subgroups provided but 'added' or 'all' requested
  if (is.null(add_subgroups) && which_subgroups %in% c("added","all")) {
    stop("`add_subgroups` must be provided when which_subgroups is 'added' or 'all'.")
  }
  
  # if add_subgroups given, ensure subgroup_lib_df has those columns
  if (!is.null(add_subgroups)) {
    missing_cols <- setdiff(add_subgroups, names(subgroup_lib_df))
    if (length(missing_cols)) {
      stop("The following 'add_subgroups' are not columns in 'subgroup_lib_df': ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  # additional subgroups
  added <- if (!is.null(add_subgroups)) add_subgroups else character()
  
  # Build full lists only if needed
  all_incl  <- c(SUBGROUPS_TO_INCLUDE, added)
  all_order <- c(SUBGROUPS_ORDER, added)
  all_names <- SUBGROUPS_TO_NAME
  if (length(added)) {
    for (x in added) all_names[x] <- x
  }
  
  # Determine which flags and ordering based on user choice
  keep_flags <- switch(
    which_subgroups,
    default = SUBGROUPS_TO_INCLUDE,
    added   = added,
    all     = all_incl
  )
  keep_order <- switch(
    which_subgroups,
    default = SUBGROUPS_ORDER,
    added   = added,
    all     = all_order
  )
  
  
  comparison_df <- comparison_df %>% 
    mutate(log2ratio = log2(ratio)) %>%
    # bring in your logical flags (one column per keep_flag)
    left_join(subgroup_lib_df, by = "Peptide") %>% 
    filter(
      .data[[group1]] >= prevalence_threshold |
        .data[[group2]] >= prevalence_threshold
    )
  
  # compute a named vector of raw p‐values, one flag at a time
  if (is.null(pvals_adj)){
    pvals_adj <- sapply(keep_flags, function(flag) {
      in_sub <- comparison_df %>% filter(.data[[flag]]) %>% pull(log2ratio)
      # “rest” is everything _not_ in that subgroup:
      out_sub <- comparison_df %>% filter(! .data[[flag]]) %>% pull(log2ratio)
      
      if (sum(!is.na(in_sub)) < 1 || sum(!is.na(out_sub)) < 1) {
        message("Skipping flag '", flag, "' due to insufficient data.")
        return(NA_real_)
      }
      
      wilcox.test(in_sub, out_sub)$p.value
    })
    pvals_adj <- p.adjust(pvals_adj, method = "BH")
  }
  
  
  subgroup_vs_rest <- tibble::tibble(
    subgroup_flag = keep_flags,       # temporary holder of the internal names
    p.adj         = pvals_adj
  ) %>%
    mutate(
      group2 = all_names[subgroup_flag],           # map to your “pretty” names
      group1 = "Complete library*"
    ) %>%
    select(group1, group2, p.adj) %>%
    add_significance("p.adj") %>%
    # enforce the plotting order:
    mutate(group2 = factor(group2, levels = keep_order))
  
  
  combos <- combn(keep_flags, 2, simplify = FALSE)
  pair_p <- lapply(combos, function(pair) {
    f1 <- pair[1]; f2 <- pair[2]
    x <- comparison_df %>% filter(.data[[f1]]) %>% pull(log2ratio)
    y <- comparison_df %>% filter(.data[[f2]]) %>% pull(log2ratio)
    data.frame(
      group1 = f1,
      group2 = f2,
      p.raw  = wilcox.test(x, y)$p.value
    )
  })
  
  pairwise_subgroups <- bind_rows(pair_p) %>%
    mutate(p.adj = p.adjust(p.raw, method = "BH"),
           # map internal flag names → display names
           group1 = all_names[group1],
           group2 = all_names[group2]) %>%
    add_significance("p.adj") %>%
    select(group1, group2, p.adj, p.adj.signif) %>% 
    mutate(
      # enforce the desired plotting order
      group1 = factor(group1, levels = keep_order),
      group2 = factor(group2, levels = keep_order)
    ) %>% 
    bind_rows(subgroup_vs_rest)
  
  # new format
  pairwise_subgroups <- pairwise_subgroups %>%
    mutate(
      p.adj.label = sapply(p.adj, format_pval)  # formatted version
    )
  
  # Plot 1: boxplot
  base_cols <- RColorBrewer::brewer.pal(9, "Paired")
  names(base_cols) <- SUBGROUPS_ORDER[1:9]
  
  # build palette depending on which_subgroups
  if (which_subgroups == "default") {
    palette_vals <- base_cols[keep_order]
  } else if (which_subgroups == "added") {
    palette_vals <- subgroup_colors[names(subgroup_colors) %in% keep_order]
  } else {
    # all: combine default + added
    extra <- subgroup_colors[names(subgroup_colors) %in% keep_order]
    palette_vals <- c(base_cols, extra)[keep_order]
  }
  
  long_ratios <- comparison_df %>%
    pivot_longer(
      cols      = all_of(keep_flags),
      names_to  = "subgroup",
      values_to = "in_subgroup"
    ) %>%
    filter(in_subgroup) %>%   # now one row per (Peptide, subgroup) that *is* in that subgroup
    mutate(
      subgroup = factor(
        all_names[subgroup], levels = keep_order
      )
    )
  
  p1 <- ggplot(long_ratios, aes(subgroup, -log2ratio)) +
    geom_violin(fill = "gray80") +    
    geom_jitter(aes(color = subgroup), width = 0.2, alpha = 0.5, size = 1) +
    stat_summary(fun = median, geom = "crossbar", color = "black", size = 0.5, fatten = 1) + #geom = "errorbar", aes(ymin = ..y.., ymax = ..y..) ) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 1.5) +
    
    #geom_boxplot(aes(fill = subgroup), outlier.shape = 21, width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray3") +
    scale_colour_manual(values = palette_vals, limits = keep_order) +
    scale_fill_manual(values = palette_vals, limits = keep_order) +
    guides(fill = "none", colour = "none") +
    labs(
      x = x_label,
      y = paste("log₂-ratio of antibody responses\nin", group1, "and", group2, sep =" ")
    ) + 
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          #legend.position = "none",
          plot.margin =  margin(2, 2, 2, 2))
  
  # Plot 2: heatmap
  p2 <- ggplot(pairwise_subgroups, aes(group1, group2, fill = p.adj.signif)) +
    geom_tile(color = "grey90", size = 0.2, show.legend = T) +
    #geom_text(aes(label = sprintf("%.1g", p.adj)), size = 1.5, colour = "black") +
    geom_text(aes(label = p.adj.label), size = 1.5, colour = "black") +
    
    scale_fill_manual(
      values = c(
        "****" = "#67001F",
        "***"  = "#B2182B",
        "**"   = "#D6604D",
        "*"    = "#F4A582",
        "ns"   = "dodgerblue1"
      ),
      limits = c("****","***","**","*","ns"),
      drop = FALSE,
      na.value = "white",
      name     = "Significance"
    ) +
    
    scale_x_discrete(limits = c("Complete library*", keep_order), position = "bottom") +
    scale_y_discrete(limits = rev(c("Complete library*", keep_order)), position = "left") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.margin =  margin(2, 2, 2, 2),
      legend.position  = "right",
      legend.margin    = margin(0, 0, 0, -25),
      legend.box.margin = margin(0, 0, 0, -25),  # pull legend 10px left
      legend.key.size  = unit(10, "pt")
    )
  p2
  # Return each plot separately
  return(list(
    boxplot = p1,
    heatmap = p2
  ))
}
