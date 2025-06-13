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

add_flag_by_patterns <- function(df,
                                 new_flag,
                                 patterns,
                                 target_cols = c("Organism_complete_name", "Description")) {
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


make_interactive_scatterplot <- function(comparison_df,
                                         group1, group2, N,
                                         highlight_cols   = NULL,
                                         highlight_colors = NULL,
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
    
    # ensure factor has all levels:
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
            "Organism: ", Organism_complete_name,"<br>",
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
    
    pvals <- sapply(highlight_cols, function(flag) {
      x <- comparison_df %>% filter( !!sym(flag) ) %>% pull(log2ratio)
      y <- comparison_df$log2ratio
      # if you prefer, you could subset y to only non-flag too:
      # y <- comparison_df %>% filter(! (!!sym(flag)) ) %>% pull(log2ratio)
      w <- wilcox.test(x, y)
      w$p.value
    })
    pvals_adj <- p.adjust(pvals, method = "BH")
    fmt_p <- formatC(pvals_adj, format="e", digits=1) # e.g. "1.2e-03" → "1.2×10⁻³"
    #fmt_p <- sub("e([-+]?)([0-9]+)$", "×10\\^\\2", fmt_p)
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
      #limits = levels_needed,    # ensures “none” is the first group drawn
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
    # no highlights requested → fall back to old categories logic
    comparison_df <- comparison_df %>%
      mutate(
        tooltip_txt = ifelse(
          categories == "not significant",
          NA_character_,
          paste0(
            "Peptide: ",  Peptide,               "<br>",
            "Desc: ",     Description,           "<br>",
            "Organism: ", Organism_complete_name,"<br>",
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
    color_scale +
    labs(
      x = paste0("% ", group1, " in whom\na peptide is significantly bound\n(n = ", N[1], ")"),
      y = paste0("% ", group2, " in whom\na peptide is significantly bound\n(n = ", N[2], ")")
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
        xaxis      = list(scaleratio = 1, scaleanchor = "y"),
        yaxis      = list(scaleratio = 1, scaleanchor = "x")
      )
    )
  } else {
    return(p)
  }
}



# Volcano plot ----
make_interactive_volcano <- function(comparison_df, group1, group2,
                                     fc_cut = 1,
                                     p_cut  = 0.05,
                                     highlight_cols   = NULL,
                                     highlight_colors = NULL,
                                     default_color    = "gray70",
                                     significant_colors = c(
                                       "not significant"                 = "dodgerblue",
                                       "significant prior correction"    = "forestgreen",
                                       "significant post FDR correction" = "firebrick"),
                                     interactive = TRUE
                                     ){

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
    
    # ensure factor has all levels:
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
            "Organism: ", Organism_complete_name,"<br>",
            group1, ": ", !!sym(group1), " / ",
            group2, ": ", !!sym(group2),       "<br>",
            "Highlight: ", highlight
          )
        )
      ) %>%
      filter(is.finite(log2ratio), is.finite(-log10(pvals_not_adj))) %>%  # drop any Inf or NaN 
      arrange(highlight)
    
    pvals <- sapply(highlight_cols, function(flag) {
      x <- comparison_df %>% filter( !!sym(flag) ) %>% pull(log2ratio)
      y <- comparison_df$log2ratio
      # if you prefer, you could subset y to only non-flag too:
      # y <- comparison_df %>% filter(! (!!sym(flag)) ) %>% pull(log2ratio)
      w <- wilcox.test(x, y)
      w$p.value
    })
    pvals_adj <- p.adjust(pvals, method = "BH")
    fmt_p <- formatC(pvals_adj, format="e", digits=1) # e.g. "1.2e-03" → "1.2×10⁻³"
    #fmt_p <- sub("e([-+]?)([0-9]+)$", "×10\\^\\2", fmt_p)
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
    comparison_df <- comparison_df %>%
      mutate(
        tooltip_txt = ifelse(
          categories == "not significant",
          NA_character_,
          paste0(
            "Peptide: ",  Peptide,               "<br>",
            "Desc: ",     Description,           "<br>",
            "Organism: ", Organism_complete_name,"<br>",
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
  
  
  p <- ggplot(comparison_df, aes(x = log2(ratio), y=-log10(pvals_not_adj))) +
    
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
                                 width   = 500, height  = 500)
    
    
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
          # bgcolor     = "rgba(255,255,255,0.8)",
          # bordercolor = "gray80",
          # borderwidth = 0,
          font        = list(size = 9)
        ),
        margin      = list(l = 0, r = 0, t = 0, b = 0, pad = 2),
        hoverlabel  = list(font = list(size = 10)),
        xaxis       = list(automargin = TRUE),
        yaxis       = list(automargin = TRUE,
                           title      = list(standoff = 10)),
        autosize    = FALSE
        
      )
    )
  } else {
    return(p)
  }
}

