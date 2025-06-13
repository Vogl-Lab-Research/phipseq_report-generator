library("stats")  # For chi-squared test
library("multcomp")  # For multiple testing correction
library("dplyr")
library("tidyr")
library("ggplot2")
library("scales")
library("stringr")
library("ggsignif")
library("plotly")
library("nnet")
library("patchwork")
library("ggpubr")
library("rstatix")
library("ggvenn")
library("ggmsa")

# Global vars

SUBGROUPS_TO_INCLUDE <- c('all', 
                          'is_PNP', 'is_patho', 'is_probio', 'is_MPA', 'is_IgA',
                          'is_bac_flagella',   'is_infect',
                          'is_IEDB_or_cntrl')
SUBGROUPS_ORDER <- c('Complete library', 
                     'Metagen antigens', 'Pathogenic strains', 'Probiotic strains', 'Microbiota strains',
                     'Antibody-coated strains',  'Flagellins', 'Infectious pathogens',
                     'IEDB/controls')
SUBGROUPS_TO_NAME <- c(
  'all' = 'Complete library',
  'is_PNP' = 'Metagen antigens',  'is_patho' = 'Pathogenic strains', 
  'is_probio' = 'Probiotic strains', 'is_MPA' = 'Microbiota strains', 'is_IgA' = 'Antibody-coated strains', 
  'is_bac_flagella' = 'Flagellins', 'is_infect' = 'Infectious pathogens', 
  'is_IEDB_or_cntrl' = 'IEDB/controls')

#######################################################
################ Helpers ##############################
#######################################################
# add_flag_by_patterns <- function(df,
#                                  new_flag,
#                                  patterns,
#                                  target_cols = c("Organism_complete_name", "Description")) {
#   # 1) Build a single case‐insensitive regex from all patterns:
#   regex_str <- paste0("(?i)", paste(patterns, collapse = "|"))
#   
#   # 2) Mutate a new logical column.  We use if_any(all_of(target_cols), ~ str_detect(...)):
#   df %>%
#     mutate(
#       !!new_flag := if_any(
#         all_of(target_cols),
#         ~ str_detect(.x, regex(regex_str))
#       )
#     )
# }


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

######################################################
############## Count Distribution ####################
######################################################
get_count_percentage_df <- function(features_target, group_col, group_cols, prevalence_threshold = 0,
                                    lib_metadata_df = NULL) {
  df <- features_target %>%
    # only samples that actually have a group
    filter(!is.na(.data[[group_col]])) %>%
    # gather 0/1 peptide calls
    pivot_longer(
      cols     = -c(SampleName, any_of(group_cols)),
      names_to = "Peptide",
      values_to= "Presence"
    ) %>%
    # compute both raw count and percent in one summarise
    group_by(Peptide, Group = .data[[group_col]]) %>%
    summarise(
      count   = sum(Presence, na.rm = TRUE),
      Percent = 100 * count / n(),
      .groups = "drop"
    ) %>%
    # spread into one column per group for Percent and one per group for Count
    pivot_wider(
      names_from  = Group,
      values_from = c(count, Percent),
      names_glue  = "{Group}_{.value}"
    )
    # bring back peptide metadata
  if (!is.null(lib_metadata_df)) {      
    df <- df %>% 
        left_join(
          lib_metadata_df %>%
            tibble::rownames_to_column("Peptide") %>%
            select(Peptide, Description, `full name`, Organism_complete_name),
          by = "Peptide"
        )
    }
  df <- df %>%
    # drop peptides that never hit prevalence threshold in any group
    filter(if_any(ends_with("_Percent"),~ . >= prevalence_threshold)) %>% 
    # strip off the "_Percent" suffix
    rename_with(
      ~ sub("_Percent$", "", .x),
      ends_with("_Percent")
    ) 
  return(df)
}


plot_enrichment_counts <- function(features_target,
                                   group_col, group_cols,
                                   prevalence_threshold = 0,
                                   custom_colors,
                                   binwidth = 1) {
  # get the count/percent table
  percentage_df <- get_count_percentage_df(
    features_target      = features_target,
    group_col            = group_col,
    group_cols           = group_cols, 
    prevalence_threshold = prevalence_threshold
  )
  
  
  # pivot to long of the *_count columns
  count_df <- percentage_df %>%
    select(ends_with("_count")) %>%
    pivot_longer(
      cols      = everything(),
      names_to  = "Cohort",
      values_to = "n_present"
    ) %>%
    mutate(Cohort = sub("_count$", "", Cohort)) %>%
    filter(n_present > 0)   # drop zero‐present peptides 
  
  # compute thresholds per cohort
  thresholds <- count_df %>%
    group_by(Cohort) %>%
    summarise(
      n_samples   = max(n_present),
      thresh      = ceiling(n_samples * 0.05),
      n_peptides5 = sum(n_present >= thresh),
      .groups     = "drop"
    )
  
  # build the plot
  p <- ggplot(count_df, aes(x = n_present, fill = Cohort)) +
    geom_histogram(
      binwidth = binwidth,
      position = "identity",
      alpha    = 0.9,
    ) +
    scale_y_log10(
      breaks = 10^(0:6),                                      # 10^0,10^1,…,10^6
      labels = trans_format("log10", math_format(10^.x)),     # render as 10^x
      expand = expansion(mult = c(0, .15))                    # a little space above
    ) +
    annotation_logticks(sides = "l", scaled = TRUE) +
    scale_fill_manual(values = custom_colors) +
    labs(
      x = "Number of samples peptide is present in",
      y = "Number of significantly bound peptides (log₁₀ scale)",
    ) +
    
    # horizontal arrowed line
    geom_segment(
      data        = thresholds,
      aes(
        x    = thresh,
        xend = n_samples,
        y    = n_peptides5,
        yend = n_peptides5
      ),
      inherit.aes = FALSE,
      linetype    = "dashed",
      color       = "black",
      size        = 0.4,
      arrow       = arrow(length = unit(0.1, "cm"), ends = "both")
    ) +
    # centred label above line
    geom_text(
      data        = thresholds,
      aes(
        x     = (thresh + n_samples) / 2,
        y     = n_peptides5,
        label = paste0(n_peptides5, " peptides in ≥5%")
      ),
      inherit.aes = FALSE,
      vjust       = -0.5,
      size        = 3
    ) +
    
    facet_wrap(~ Cohort, ncol = 2, scales = "free_x") +
    theme_bw(base_size = 10) +
    theme(
      legend.position   = "none",
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      panel.grid.major  = element_line(color = "grey90"),
      panel.grid.minor  = element_blank()
    )
  
  return(p)
}


######################################################
####### plot  enrichment and diversity################
######################################################
plot_groups_boxplots <- function(data, group_col, values_col, custom_colors, pairwise_comparisons, label_axis = NA) {
  # Convert grouping column name to symbol
  group_sym <- sym(group_col)
  values_sym <- sym(values_col)
  
  # Drop rows where either group or value is NA
  data <- data %>%
    filter(!is.na(!!group_sym), !is.na(!!values_sym))
  
  # Summarize counts by that group
  df_counts <- data %>%
    group_by(!!group_sym) %>%
    summarize(sample_count = n(), .groups = "drop")
  
  # Create x-axis labels using the dynamic variable.
  # Since group_col is a string, we access it directly on df_counts.
  x_labels <- setNames(
    paste0(df_counts[[group_col]], "\n(n = ", df_counts$sample_count, ")"),
    df_counts[[group_col]]
  )
  
  # Build the plot:
  p <- ggplot(data, aes(x = !!group_sym, y = !!values_sym, fill = !!group_sym)) +
    geom_boxplot(show.legend = FALSE) +  # Hide legend if desired
    geom_jitter(color = "black", size = 1, width = 0.2, alpha = 0.3, show.legend = FALSE) +  
    scale_fill_manual(values = custom_colors) +  # Assign custom colors
    scale_x_discrete(labels = x_labels) +         # Use the custom labels
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5),
      panel.grid = element_blank()
    ) +
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               comparisons = pairwise_comparisons, 
                               label = "p.signif",  # Display significance level (e.g., * or **)
                               hide.ns = FALSE,     # Option to hide non-significant comparisons
                               size = 4)
  
  # If label_axis is not NA, add custom axis labels. Assuming label_axis is a vector of length 2:
  if (!is.na(label_axis[1]) && !is.na(label_axis[2])) {
    p <- p + labs(
      x = label_axis[1],
      y = label_axis[2]
    )
  }
  
  return(p)
}

#################################
## plot sex/age distribution#####
#################################
test_sex_age_distribution <- function(data, 
                                      group_col, 
                                      age_col = "Age", 
                                      sex_reg = "Sex", 
                                      sex_ctg = "Sex_ctg") {
  # Chi-square test using the grouping variable (group_col) and the categorical sex column (sex_ctg)
  cat("Chi-square test result:\n")
  chisq_result <- chisq.test(table(data[[group_col]], data[[sex_ctg]]))
  print(chisq_result)
  
  # Multinomial regression: Predict the grouping variable based on Age and the numeric Sex variable.
  # The formula is constructed as: `group_col` ~ Age * Sex
  cat("\nMultinomial Regression Summary:\n")
  fmla_multinom <- as.formula(paste0("`", group_col, "` ~ ", age_col, " * ", sex_reg))
  multinom_model <- multinom(fmla_multinom, data = data)
  print(summary(multinom_model))
  
  # Two-way ANOVA: Predict Age based on group and the categorical sex variable.
  # Here the formula is: Age ~ `group_col` * Sex_ctg
  cat("\nANOVA Summary:\n")
  fmla_aov <- as.formula(paste0(age_col, " ~ `", group_col, "` * ", sex_ctg))
  aov_model <- aov(fmla_aov, data = data)
  tidy_aov <- broom::tidy(aov_model)
  tidy_aov$p.value <- round(tidy_aov$p.value, 3)
  
  # Display the ANOVA results in a neat table
  print(knitr::kable(tidy_aov, digits = 3, 
                     caption = "Two-way ANOVA Results for Age by Group and Sex"))
  
  # Optionally, return a list with the results:
  return(list(chisq = chisq_result,
              multinom_summary = summary(multinom_model),
              aov_results = tidy_aov))
}

plot_sex_age_distribution <- function(data, 
                                      group_col, 
                                      age_col = "Age_group", 
                                      sex_col = "Sex_ctg", 
                                      custom_colors) {
  
  # Convert the grouping column name to a symbol for tidy evaluation.
  group_sym <- sym(group_col)
  
  # Step 1: Summarize counts per combination of group, sex and age group.
  # Here we assume that you want to mirror counts so that Male counts become negative.
  data_summary <- data %>%
    # 0) drop any samples with missing Sex or missing Age or missing group
    filter(!is.na(!!group_sym), !is.na(.data[[sex_col]]), !is.na(.data[[age_col]])) %>%
    group_by(!!group_sym, .data[[sex_col]], .data[[age_col]]) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(count = ifelse(.data[[sex_col]] == "Male", -count, count))
  
  
  # Step 2: Compute overall counts by group (for generating facet labels).
  df_counts <- data %>%
    filter(!is.na(!!group_sym), !is.na(.data[[sex_col]]), !is.na(.data[[age_col]])) %>%
    group_by(!!group_sym) %>%
    summarize(sample_count = n(), .groups = "drop")
  
  # Create x-axis labels (here, using the group variable value and overall sample count).
  x_labels <- setNames(
    paste0(df_counts[[group_col]], "\n(n = ", df_counts$sample_count, ")"),
    df_counts[[group_col]]
  )
  
  
  p <- ggplot(data_summary, aes(x = ifelse(!!group_sym == levels(factor(.data[[group_sym]]))[1], count, -count), y = .data[[age_col]], fill = .data[[sex_col]])) +
    geom_bar(stat = "identity", position = "identity", width = 0.85) +
    scale_x_continuous(
      labels = function(x) ifelse(x %in%  seq(-300, 300, by = 3), abs(x), ""),
      breaks = seq(-300, 300, by = 3)
    ) +
    scale_fill_manual(values = custom_colors) +
    # Facet by the chosen grouping variable.
    facet_grid(as.formula(paste0("~ `", group_col, "`")), scales = "free_x", space = "free_x", 
               labeller = as_labeller(x_labels)) +
    labs(
      x = "Counts",
      y = "Age Group",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 11),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.major.y = element_blank(),
      axis.text.y.left = element_text(size = 10, face = "italic"),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.box.spacing = unit(0.2, "lines")
    )
  
  return(p)
}


####################################
##########Scatterplot###############
####################################
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


# make_interactive_scatterplot <- function(comparison_df, group1, group2, N,
#                                          highlight_col   = NULL,
#                                          scatter_colors = c(
#                                            "not significant"                 = "dodgerblue",
#                                            "significant prior correction"    = "forestgreen",
#                                            "significant post FDR correction" = "firebrick"
#                                          )){
#   # 2) Decide on the aesthetic mapping for color:
#   if (!is.null(highlight_col) && highlight_col %in% names(comparison_df)) {
#     comparison_df <-  comparison_df %>%
#       mutate(
#         # Build a tooltip text only for the “significant” categories
#         tooltip_txt = ifelse(
#           !!sym(highlight_col) == FALSE,
#           NA_character_,
#           paste0(
#             "Peptide: ",     Peptide,      "<br>",
#             "Desc: ",        Description,  "<br>",
#             "Organism: ",    Organism_complete_name, "<br>",
#             group1, ": ",    !!sym(group1), " / ",
#             group2, ": ",    !!sym(group2), "<br>"
#           )
#         )
#       ) %>%
#       filter(is.finite(log2(ratio)), is.finite(-log10(pvals_not_adj)))  # drop any Inf or NaN
#     # We color by highlight_col (TRUE/FALSE).
#     color_aes   <- aes(color = !!sym(highlight_col), text = tooltip_txt)
#     color_scale <- scale_color_manual(values = scatter_colors, breaks = "TRUE",
#                                       labels = highlight_col) #name = highlight_col)
#     # Legend location: top‐right inside
#     legend_theme <- theme(
#       legend.position      = c(0.97, 0.97),
#       legend.justification = c(1, 1),
#       legend.background    = element_rect(fill = alpha("white", 0.8), color = "gray80"),
#       legend.key.size      = unit(10, "pt"),
#       legend.text          = element_text(size = 8),
#       legend.title         = element_text(size = 9, face = "bold")
#     )
#     show_legend <- TRUE
#   } else {
#     comparison_df <-  comparison_df %>%
#       mutate(
#         # Build a tooltip text only for the “significant” categories
#         tooltip_txt = ifelse(
#           categories == "not significant",
#           NA_character_,
#           paste0(
#             "Peptide: ",     Peptide,      "<br>",
#             "Desc: ",        Description,  "<br>",
#             "Organism: ",    Organism_complete_name, "<br>",
#             group1, ": ",    !!sym(group1), " / ",
#             group2, ": ",    !!sym(group2), "<br>"
#           )
#         )
#       ) %>%
#       filter(is.finite(log2(ratio)), is.finite(-log10(pvals_not_adj)))  # drop any Inf or NaN
#     
#     color_aes   <- aes(color = categories, text = tooltip_txt)
#     color_scale <- scale_color_manual(values = scatter_colors, name = NULL)
#     
#     legend_theme <- theme(legend.position = "none")
#     show_legend <- FALSE
#   }
#   
# 
#   
#   # Generate scatter plot for the comparison
#   p <- ggplot(comparison_df, aes(x = !!sym(group1), y = !!sym(group2)) ) +
#     
#     geom_point(color_aes, alpha = 0.65) +
#     color_scale +
#     labs(x = paste0("% ", group1, " in whom\na peptide is significantly bound\n(n = ",
#                     N[1], ")"),
#          y = paste0("% ", group2, " in whom\na peptide is significantly bound\n(n = ", 
#                     N[2], ")")
#     ) +
#     theme_bw(base_size = 12) +  # Use a minimal theme for elegance and set a base font size
#     theme(
#       #legend.position = "none",
#       #aspect.ratio = 1,
#       panel.grid.major = element_blank(),  # Remove major grid lines
#       panel.grid.minor = element_blank(),  # Remove minor grid lines
#       panel.border = element_rect(colour = "black", fill = NA),  # Keep border if desired
#       plot.margin = margin(t = 10, r = 15, b = 15, l = 10, unit = "pt"),  # Adjust margins in point
#       axis.text.y.left = element_text(size = 10, face = "italic"),  # Style y-axis labels
#       #axis.ticks.y.right = element_blank(),  # Remove right-side y-axis ticks
#       axis.title.x = element_text(face = "bold"),  # Add margin and bold to x-axis title
#       axis.title.y = element_text(face = "bold")  # Bold y-axis title
#     ) +
#     legend_theme
#   
#   # Convert to plotly for interactivity
#   interactive_plot <- ggplotly(p, tooltip = c("text")) %>% 
#     layout(
#       showlegend = show_legend,
#       margin = list(l = 80, r = 80, b = 80, t = 80, pad = 0),
#       height = 550,  # Define the plot height in pixels
#       width = 550,
#       hoverlabel = list(font = list(size = 10)),
#       xaxis = list(#range = c(-2, 102),
#         scaleratio = 1,
#         scaleanchor = "y"
#       ),
#       yaxis = list(#range = c(-2, 102), 
#         scaleratio = 1,
#         scaleanchor = "x"),
#       legend        = if (show_legend) list(
#         x           = 0.02,
#         y           = 0.98,
#         xanchor     = "left",
#         yanchor     = "top",
#         bgcolor     = "rgba(255,255,255,0.8)",
#         bordercolor = "gray80",
#         borderwidth = 0,
#         font        = list(size = 9)
#       ) else NULL
#     )    
#   return(interactive_plot)
# }

group_test_analysis <- function(percentage_group_test_list, num_samples_per_group,
                                         prevalence_threshold = 5) {
  results_list <- list()
  
  # Loop through each group_test column in the list
  for (group_col in names(percentage_group_test_list)) {
    # Get the dataframe for the current group_test
    df <- percentage_group_test_list[[group_col]]
    
    # Identify unique groups within this group_test column
    groups <- colnames(df)[grepl("_count$", colnames(df))] %>%
      sub("_count$", "", .)
    
    # Ensure there are at least two groups for pairwise comparison
    if (length(groups) < 2) next
    
    
    # Generate pairwise comparisons for each group_test column
    comparison_results <- list()
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        group1 <- groups[i]
        group2 <- groups[j]
        
        N1 <- num_samples_per_group[[group_col]][[group1]]
        N2 <- num_samples_per_group[[group_col]][[group2]]

        # ——— apply pairwise prevalence filter ——————————————
        df_pair <- df %>%
          # keep only peptides with >= threshold in at least one of the two
          filter(
            .data[[ group1 ]] >= prevalence_threshold |
              .data[[ group2 ]] >= prevalence_threshold
          )
        if (nrow(df_pair) == 0) {
          warning("No peptides pass the prevalence filter for ", 
                  group1, " vs ", group2)
          next
        }
        
        # Calculate p-values and ratios for each pair of groups
        delta_ratio_vals <- numeric(nrow(df_pair))
        ratio_vals       <- numeric(nrow(df_pair))
        pvals_chisq      <- numeric(nrow(df_pair))
        epsilon_prop      <- 0.5
        epsilon_delta    <- 1
        
        # Loop over each row in the filtered dataframe
        for (k in 1:nrow(df_pair)) {
          # Extract the values for the two groups being compared
          val1 <- df_pair[[paste0(group1, "_count")]][k]
          val2 <- df_pair[[paste0(group2, "_count")]][k]
          
          # Construct the contingency table using counts from num_samples_per_group
          chitable <- matrix(c(val1 + 1, N1 - val1 + 1, 
                               val2 + 1, N2 - val2 + 1), 
                             nrow = 2, byrow = TRUE)
          
          # Chi-squared test 
          #pvals_chisq[k] < chisq.test(chitable)$p.value

          # Fisher test 
          test_result <- fisher.test(chitable)
          pvals_chisq[k] <- test_result$p.value
          
    
          # ratio (with 0.5 pseudocount only for zeros)
          val1_prop <- (val1 + ifelse(val1 == 0, epsilon_prop, 0)) / N1
          val2_prop <- (val2 + ifelse(val2 == 0, epsilon_prop, 0)) / N2
          ratio_vals[k] <- (val1_prop / val2_prop)
          

          # Delta-ratio (with 1 pseudocount only for zeros)
          val1_delta <- (val1 + ifelse(val1 == 0, epsilon_delta, 0)) / N1
          val2_delta <- (val2 + ifelse(val2 == 0, epsilon_delta, 0)) / N2
          if (val1_delta >= val2_delta) {
            delta_ratio_vals[k] <- val1_delta / val2_delta - 1
          } else {
            delta_ratio_vals[k] <- -(val2_delta / val1_delta - 1)
          }
        }
        
        # Add p-values and ratios to the dataframe
        comparison_df <- df_pair %>%
          mutate(
            Delta_ratio = delta_ratio_vals,
            ratio = ratio_vals,
            pvals_not_adj = pvals_chisq
          ) %>%
          mutate(
            passed_not_adj = ifelse(pvals_not_adj < 0.05, TRUE, FALSE),
            pvals_bh = p.adjust(pvals_not_adj, method = "BH"),
            passed_bh = p.adjust(pvals_not_adj, method = "BH") < 0.05,#,
            categories = case_when(
              passed_bh                     ~ "significant post FDR correction",
              passed_not_adj & !passed_bh   ~ "significant prior correction",
              TRUE                          ~ "not significant"
            ),
            
            #Significant_bh = ifelse(passed_bh, "Yes", "No")
          )  %>%
          dplyr::arrange(
            Delta_ratio,            # then by direction of change
            !passed_bh,            # TRUE (significant) comes first
            pvals_bh,              # then smallest BH-adjusted p-values
            pvals_not_adj,        # then smallest raw p-values
          ) %>%
          dplyr::select(Peptide, Description, `full name`, everything())
        

        # Generate scatter plot for the comparison
        p <- make_interactive_scatterplot(comparison_df = comparison_df, group1 = group1, group2 = group2, N = c(N1,N2))
        
        # Store the plot and results in the list
        comparison_results[[paste0(group1, "_vs_", group2)]] <- list(
          plot = p,
          comparison_df = comparison_df
        )
      }
    }
    
    # Store all comparison results for this group_test
    results_list[[group_col]] <- comparison_results
  }
  
  return(results_list)
}


#########################################
############ Volcano plot ###############
#########################################
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
    plotly_width <- 500
    margins <- list(l = 0, r = 0, t = 0, b = 0, pad = 2)
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

#########################################
###############MDS#######################
#########################################


plot_mds <- function(features_target, group_col, 
                     custom_colors,  group_cols =NA,
                     method = "jaccard",
                     permutations = 999) {
  
  # Full matrix sames as exist for the corresponding group and transposed
  binary_data_all <- features_target %>%
    filter(!is.na(.data[[group_col]])) %>%       # Only keep rows with non-NA group
    dplyr::select(-any_of(group_cols)) %>%
    tibble::column_to_rownames("SampleName")
  
  # Distance matrices
  dist_all <- vegan::vegdist(binary_data_all, method = method)
  mds_all <- cmdscale(dist_all, k = 2, eig = TRUE)
  var_exp_all <- round(100 * mds_all$eig[1:2] / sum(mds_all$eig), 2)
  
  
  # MDS dataframe
  mds_df_all <- as.data.frame(mds_all$points) %>%
    tibble::rownames_to_column("SampleName") %>%
    left_join(features_target %>% select(any_of(c("SampleName", group_col))), by = "SampleName")
  
  # Permanova test
  permanova <-   vegan::adonis2(dist_all ~ mds_df_all[[group_col]], permutations = permutations)
  p_value <- permanova$`Pr(>F)`[1]
  
  # Centroids for plotting
  centroids <- mds_df_all %>%
    group_by(.data[[group_col]]) %>%
    summarise(V1 = mean(V1), V2 = mean(V2), .groups = "drop")
  
  # Main plot
  p <- ggplot(mds_df_all, aes(x = V1, y = V2, fill = .data[[group_col]])) +
    geom_point(size = 3, alpha = 0.5, shape = 21, color = "black",  stroke = 0, aes(text = paste("Sample:", SampleName))) +
    geom_point(data = centroids, aes(x = V1, y = V2, fill = .data[[group_col]]),
               show.legend = FALSE, size = 5, shape = 21, color = "black") +
    scale_fill_manual(values = custom_colors) +
    labs(
      x = paste0("MDS 1 (", var_exp_all[1], "%)"),
      y = paste0("MDS 2 (", var_exp_all[2], "%)"),
      title = paste0("MDS (PCoA) with Jaccard Distance\nPERMANOVA p = ", round(p_value, 4),
                     " | permutations = ", permutations)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  p
}

####################################
#############PCA ###################
####################################
plot_pca <- function(fold_df, group_col, custom_colors,
                     group_cols = NULL, prevalence_cutoff = 0.05,
                     log_transform = TRUE) {
  
  # 1) Identify peptide columns
  peptide_cols <- setdiff(
    names(fold_df),
    c("SampleName", group_col, group_cols)
  )
  
  # 2) Find peptides with ≥ cutoff prevalence in any group
  keep_peptides <- fold_df %>%
    filter(!is.na(.data[[group_col]])) %>%
    tidyr::pivot_longer(all_of(peptide_cols),
                        names_to  = "peptide",
                        values_to = "fc") %>%
    mutate(present = fc != 0) %>%
    group_by(.data[[group_col]], peptide) %>%
    summarise(prevalence = mean(present), .groups = "drop") %>%
    group_by(peptide) %>%
    filter(any(prevalence >= prevalence_cutoff)) %>%
    pull(peptide) %>%
    unique()
  
  # 3) Build filtered matrix (samples × kept peptides)
  mat <- fold_df %>%
    filter(!is.na(.data[[group_col]])) %>%
    select(SampleName, all_of(keep_peptides)) %>%
    tibble::column_to_rownames("SampleName") %>%
    as.matrix()
  
  # 4) Optionally log2-transform fold-changes
  if (log_transform) {
    mat <- log2(mat + 1)  
    # +1 avoids log(0); use another pseudocount if you prefer
  }
  
  # 5) PCA
  pca_res <- prcomp(mat, center = TRUE, scale. = FALSE)
  var_exp  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  pca_df <- as.data.frame(pca_res$x[,1:2]) %>%
    tibble::rownames_to_column("SampleName") %>%
    left_join(
      fold_df %>% select(SampleName, any_of(group_col)),
      by = "SampleName"
    )
  
  # 6) Plot PC1 vs PC2
  ggplot(pca_df, aes(x = PC1, y = PC2, fill = .data[[group_col]])) +
    geom_point(shape = 21, color = "black", size = 3, alpha = 0.6) +
    scale_fill_manual(values = custom_colors) +
    labs(
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)"),
      fill = group_col
    ) +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

####################################
############Similarity pot##########
####################################
compute_patient_correlation <- function(data_matrix, metadata, 
                                        ind_id_col = "ind_id", 
                                        samples_t1_col = "SampleName_t1", 
                                        samples_t2_col = "SampleName_t2",
                                        method = "pearson") {
  
  # Order the metadata by patient_id to align pairs consistently.
  metadata <- metadata[order(metadata[[ind_id_col]]), ]
  
  # Extract baseline and t2 sample names from metadata as character vectors.
  t1_samples <- as.character(metadata[[samples_t1_col]])
  t2_samples       <- as.character(metadata[[samples_t2_col]])
  
  # Create an empty correlation matrix with rows = baseline samples and columns = t2 samples.
  # Note: The dimnames come directly from the metadata values.
  corr_mat <- matrix(NA, nrow = length(t1_samples), ncol = length(t2_samples),
                     dimnames = list(t1_samples, t2_samples))
  
  # Loop over each pair of baseline and t2 sample names.
  for (i in seq_along(t1_samples)) {
    for (j in seq_along(t2_samples)) {
      b_sample <- t1_samples[i]
      t_sample <- t2_samples[j]
      
      # If one of the sample names is missing (or is NA/empty), leave that cell as NA.
      if (is.na(b_sample) || is.na(t_sample) || b_sample == "" || t_sample == "") {
        corr_mat[i, j] <- NA
      } else if (!(b_sample %in% rownames(data_matrix)) || !(t_sample %in% rownames(data_matrix))) {
        # If the sample is not found in the data matrix, set correlation to NA.
        corr_mat[i, j] <- NA
      } else {
        # Otherwise, extract the peptide profiles (each is a numeric vector)
        vec_b <- as.numeric(data_matrix[b_sample, ])
        vec_t <- as.numeric(data_matrix[t_sample, ])
        
        # If either profile has zero variance, correlation is undefined (set to NA)
        if (sd(vec_b, na.rm = TRUE) == 0 || sd(vec_t, na.rm = TRUE) == 0) {
          corr_mat[i, j] <- NA
        } else {
          # Compute the correlation between the two samples.
          corr_mat[i, j] <- cor(vec_b, vec_t, use = "pairwise.complete.obs", method = method)
        }
      }
    }
  }
  
  return(corr_mat)
}

plot_correlation <- function(phiseq_df, metadata,
                             ind_id_col = "ind_id", samples_t1_col = "SampleName_t1", samples_t2_col = "SampleName_t2",
                             label_x = "t1", label_y = "t2",
                             sort_by_status = FALSE, pre_status_col = "status_t1", post_status_col = "status_t2",
                             method = "pearson", require_both_timepoints = FALSE) {
  
  # Optionally filter metadata to only include individuals with both timepoints
  use_fallback <- FALSE
  if (require_both_timepoints) {
    tmp_meta <- metadata %>%
      filter(g1_exists & g2_exists)
    if (nrow(tmp_meta) == 0) {
      print("No individuals with both timepoints. Falling back to those with either timepoint.")
      tmp_meta <- metadata
      use_fallback <- TRUE
    }
    metadata <- tmp_meta
  }
  
  # Compute correlation matrix for this subgroup
  corr_matrix <- compute_patient_correlation(
    data_matrix = t(phiseq_df),
    metadata = metadata,
    ind_id_col = ind_id_col,
    samples_t1_col = samples_t1_col,
    samples_t2_col = samples_t2_col,
    method = method
  )
  
  if (sort_by_status == TRUE){
    metadata <- metadata %>% 
      arrange((!!sym(pre_status_col)), (!!sym(post_status_col)), (!!sym(ind_id_col))) %>% 
      slice(rev(row_number()))
  }
  
  # 2. Extract the sample names in the desired order.
  ordered_ids <- metadata[[ind_id_col]]
  
  
  # Format for plotting
  rownames(corr_matrix) <- as.character(ordered_ids)
  colnames(corr_matrix) <- as.character(ordered_ids)
  
  cor_df <- reshape2::melt(corr_matrix)
  
  if (sort_by_status == TRUE){
    cor_df <- cor_df %>%
      left_join(metadata %>%  select(patient_id = !!sym(ind_id_col), pre_status = !!sym(pre_status_col)),
                by = c("Var1" = ind_id_col)) %>%
      left_join(metadata %>% select(patient_id = !!sym(ind_id_col), post_status = !!sym(post_status_col)),
                by = c("Var2" = ind_id_col))
  }
  
  cor_df$Var1 <- factor(cor_df$Var1, levels = ordered_ids)
  cor_df$Var2 <- factor(cor_df$Var2, levels = ordered_ids)
  
  if (use_fallback) {
    cor_df <- cor_df %>% filter(!is.na(value))
  }
  
  # Plot
  p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(colour = "gray60", linewidth = 0.25) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(0, 1),
                         na.value = "gray", name = "Pearson\nCorrelation") +
    coord_fixed() +
    theme_bw() +
    theme(
      axis.text.x = element_text(size=6.5, angle = 45, vjust = 1, hjust = 0.5),
      axis.text.y = element_text(size = 6.5),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    #scale_x_discrete(breaks = col_breaks) +
    #scale_y_discrete(breaks = row_breaks) 
    labs(
      x = paste0(label_x," Ig epitope repertoire ", "\n(individual ID)"),
      y = paste0(label_y," Ig epitope repertoire ", "\n(individual ID)")) 
  
  return(p)
}


plot_correlation_distribution <- function(phiseq_df, metadata,
                                          ind_id_col = "ind_id",
                                          samples_t1_col = "SampleName_t1",
                                          samples_t2_col = "SampleName_t2",
                                          label_x = "t1", label_y = "t2",
                                          method = "pearson",
                                          bins = 20,
                                          require_both_timepoints = FALSE) {
  use_fallback <- FALSE
  if (require_both_timepoints) {
    tmp_meta <- metadata %>%
      filter(g1_exists & g2_exists)
    if (nrow(tmp_meta) == 0) {
      print("No individuals with both timepoints. Falling back to those with either timepoint.")
      tmp_meta <- metadata
      use_fallback <- TRUE
    }
    metadata <- tmp_meta
  }
  
  corr_matrix <- compute_patient_correlation(
    data_matrix = t(phiseq_df),
    metadata = metadata,
    ind_id_col = ind_id_col,
    samples_t1_col = samples_t1_col,
    samples_t2_col = samples_t2_col,
    method = method
  )
  
  ordered_ids <- metadata[[ind_id_col]]
  rownames(corr_matrix) <- as.character(ordered_ids)
  colnames(corr_matrix) <- as.character(ordered_ids)
  
  cor_df <- reshape2::melt(corr_matrix) %>%
    mutate(pair_type = ifelse(Var1 == Var2, "matched", "random"))
  
  
  if (use_fallback) {
    cor_df <- cor_df %>% filter(!is.na(value))
  }
  
  # Separate matched and random
  matched_df <- cor_df %>% filter(pair_type == "matched")
  random_df  <- cor_df %>% filter(pair_type == "random")
  
  # Initialize an empty plot
  p <- ggplot()
  
  has_matched <- nrow(matched_df) > 0
  has_random  <- nrow(random_df) > 0
  
  if (has_random) {
    random_hist <- ggplot(random_df, aes(x = value)) + geom_histogram(bins = bins)
    random_data <- ggplot_build(random_hist)$data[[1]]
    p <- p + geom_col(data = random_data, aes(x = x, y = y),
                      width = random_data$width, fill = "palegreen4", alpha = 0.4)
  }
  
  if (has_matched) {
    matched_hist <- ggplot(matched_df, aes(x = value)) + geom_histogram(bins = bins)
    matched_data <- ggplot_build(matched_hist)$data[[1]]
    
    # Fallback to default if no matched histogram data was produced
    has_matched_hist <- !is.null(matched_data) && nrow(matched_data) > 0 && max(matched_data$y) > 0
    
    if (has_random && has_matched_hist) {
      scale_factor <- max(random_data$y) / max(matched_data$y)
      matched_data <- matched_data %>%
        mutate(y_scaled = y * scale_factor)
      
      p <- p + geom_col(data = matched_data, aes(x = x, y = y_scaled),
                        width = matched_data$width, fill = "dodgerblue3", alpha = 0.4) +
        scale_y_continuous(
          name = "Random Pair Count",
          sec.axis = sec_axis(~ . / scale_factor, name = "Matched Pair Count")
        )
    } else if (has_matched_hist) {
      # Only matched → plot as primary y-axis
      p <- p + geom_col(data = matched_data, aes(x = x, y = y),
                        width = matched_data$width, fill = "dodgerblue3", alpha = 0.4) +
        scale_y_continuous(name = "Matched Pair Count")
    }
  }
  
  p <- p +
    scale_x_continuous(name = paste0("Pearson correlation of\n", label_x, " vs ", label_y)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(p)
}


###################################################
############ ratios subgroups######################
###################################################

# plot_ratios_by_subgroup <- function(comparison_df, group1, group2, subgroup_lib_df, prevalence_threshold = 5) {
#   # Join peptide subgroup flags
#   long_ratios <- comparison_df %>%
#     mutate(ratio = log2(ratio)) %>% 
#     left_join(subgroup_lib_df, by = "Peptide") %>%
#     tidyr::pivot_longer(cols = all_of(SUBGROUPS_TO_INCLUDE),
#                         names_to = "subgroup_flag",
#                         values_to = "belongs_to_group") %>%
#     dplyr::filter(belongs_to_group) %>%
#     mutate(subgroup = factor(SUBGROUPS_TO_NAME[subgroup_flag],
#                              levels = SUBGROUPS_ORDER)) %>%
#     dplyr::filter(.data[[group1]] >= prevalence_threshold | .data[[group2]] >= prevalence_threshold)
#   
#   
#   # Skip empty plots
#   if (nrow(long_ratios) < 10) {
#     warning("Not enough data for ", group1, " vs ", group2)
#     return(NULL)
#   }
#   
#   # Get pairwise subgroup comparisons
#   subgroup_labels <- levels(long_ratios$subgroup)
#   pairwise_combos <- combn(subgroup_labels, 2, simplify = FALSE)
#   
#   # Get significant pairs
#   sig_comparisons <- ggpubr::compare_means(
#     formula = ratio ~ subgroup,
#     data = long_ratios,
#     method = "wilcox.test",
#     comparisons = pairwise_combos,
#     p.adjust.method = "bonferroni"
#   ) %>% filter(p.adj < 0.01)
#   
#   # Format for stat_compare_means
#   sig_pairs <- lapply(seq_len(nrow(sig_comparisons)), function(i) {
#     c(sig_comparisons$group1[i], sig_comparisons$group2[i])
#   })
#   
#   # Plot
#   p <- ggplot(long_ratios, aes(x = subgroup, y = ratio, fill = subgroup)) +
#     geom_boxplot(outlier.shape = 21, outlier.size = 1, width = 0.6) +
#     scale_fill_brewer(palette = "Paired") +
#     ggpubr::stat_compare_means(
#       method = "wilcox.test",
#       comparisons = sig_pairs,
#       p.adjust.method = "BH",
#       label = "p.signif",
#       hide.ns = TRUE,
#       size = 4
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.position = "none"
#     ) +
#     labs(
#       x = "Subgroups of the antigen library",
#       y = paste("log-ratio of antibody responses\nin", group1, "and", group2, sept=" ")
#     )
#   
#   return(p)
# }

plot_ratios_by_subgroup <- function(comparison_df,
                              group1, group2,
                              subgroup_lib_df, custom_colors,
                              subgroup_colors = NULL,
                              prevalence_threshold = 5,
                              #min_peptides = 5, 
                              add_subgroups = NULL) {
  

  subgroups_to_include = SUBGROUPS_TO_INCLUDE
  subgroups_order = SUBGROUPS_ORDER
  subgroups_to_name = SUBGROUPS_TO_NAME

  if (!is.null(add_subgroups)) {
    for (flag in add_subgroups) {
      subgroups_to_include = c(subgroups_to_include, flag)
      subgroups_order = c(subgroups_order, flag)
      subgroups_to_name[flag] = flag 
    }
  }
  
  # 1) build the long table
  long_ratios <- comparison_df %>%
    mutate(ratio = log2(ratio)) %>% 
    left_join(subgroup_lib_df, by = "Peptide") %>%
    tidyr::pivot_longer(
      cols      = all_of(subgroups_to_include),
      names_to  = "subgroup_flag",
      values_to = "in_subgroup"
    ) %>%
    dplyr::filter(in_subgroup) %>%
    mutate(
      subgroup = factor(
        subgroups_to_name[subgroup_flag],
        levels = subgroups_order
      )
    ) %>%
    dplyr::filter(.data[[group1]] >= prevalence_threshold |
                    .data[[group2]] >= prevalence_threshold) %>%
    select(Peptide, subgroup, ratio)
  
  if(nrow(long_ratios) < 10) {
    warning("Too few points for ", group1, " vs ", group2)
    return(NULL)
  }
  
  posthoc <- long_ratios %>%
    #dunn_test(ratio ~ subgroup, p.adjust.method = "BH") %>%
    pairwise_wilcox_test(ratio ~ subgroup, p.adjust.method = "BH") %>%  
    #dplyr::filter(n1 >= min_peptides, n2 >= min_peptides) %>%
    add_significance("p.adj") #%>%
    # add_xy_position(
    #   data    = long_ratios,
    #   formula = ratio ~ subgroup,
    #   step.increase = diff(range(long_ratios$ratio)) * 0.1
    # )
  
  # Manually compute y positions so they stack
  # if (nrow(posthoc) > 0) {
  #   top   <- max(long_ratios$ratio, na.rm = TRUE)
  #   # small offset proportional to the data range
  #   offset <- diff(range(long_ratios$ratio, na.rm = TRUE)) * 0.1
  #   posthoc_filter <- posthoc %>%
  #     dplyr::filter(p.adj < 0.01) %>% 
  #     arrange(p.adj) %>%
  #     mutate(y.position = top + seq_len(n()) * offset)
  # }
  
  if (is.null(subgroup_colors)){
    color_scale <- scale_fill_brewer(palette = "Paired") 
  } else{
    color_scale <- scale_color_manual(values = subgroup_colors)
  }
  
  

  # grab the first up to 9 colours from Paired
  base_cols <- RColorBrewer::brewer.pal(9, "Paired")
  names(base_cols) <- subgroups_order[1:9]
  
  if (!is.null(subgroup_colors)) {
    # only keep colours for subgroups we actually plot
    extra <- subgroup_colors[names(subgroup_colors) %in% subgroups_order]
    palette_vals <- c(base_cols, extra)[subgroups_order]
  } else {
    palette_vals <- base_cols[subgroups_order]
  }
  
  p1 <- ggplot(long_ratios, aes(subgroup, ratio, fill = subgroup)) +
    geom_boxplot(outlier.shape = 21, width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray1") +
    scale_fill_manual(values = palette_vals, limits = subgroups_order) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
      ,    plot.margin =  margin(2, 2, 2, 2),
    ) +
    labs(
      x = "Subgroups of the antigen library",
      y = paste("log₂-ratio of antibody responses\nin", group1, "and", group2, sept=" ")
    )   #+
    # # global KW p-value
    # stat_compare_means(
    #   method  = "kruskal.test",
    #   label   = "p.format",
    #   label.x = 4, #min(long$ratio) * 1.02,
    #   #label.x.npc = "centre",
    #   label.y.npc = "bottom",
    #   size        = 3
    # )
  
  #  Add pairwise stars if any
  # if (nrow(posthoc_filter) > 0) {
  #   p1 <- p1 +
  #     stat_pvalue_manual(
  #       posthoc_filter,
  #       label       = "p.adj.signif",
  #       xmin        = "group1",
  #       xmax        = "group2",
  #       y.position  = "y.position",
  #       tip.length  = 0.01,
  #       inherit.aes = FALSE
  #     )
  #}
  
  p2 <- ggplot(posthoc, aes(x = group1, y = group2, fill = p.adj.signif)) +
    geom_tile(color = "grey90", size = 0.2) +
    geom_text(aes(label = sprintf("%.1g", p.adj)),
              color = "black", size = 1.5) +
    scale_fill_manual(
      values = c(
        "****"  = "#67001F",
        "***"   = "#B2182B",
        "**"    = "#D6604D",
        "*"     = "#F4A582",
        "ns"    = "dodgerblue1"
      ),
      na.value = "white",
      name     = "Significance"
    ) +
    # Put x-axis on the bottom:
    scale_x_discrete(
      limits  = subgroups_order,
      position = "bottom"
    ) +
    # Reverse y so first factor is at the top:
    scale_y_discrete(
      limits   = rev(subgroups_order),
      position = "left"
    ) +
    labs(x = NULL, y = NULL) +
    #coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.margin =  margin(2, 2, 2, 2),
      legend.position  = "right",
      legend.margin    = margin(0, 0, 0, -25),
      legend.box.margin = margin(0, 0, 0, -25),  # pull legend 10px left
      legend.key.size  = unit(10, "pt"),
      #legend.title     = element_text(size = 8),
      #legend.text      = element_text(size = 7)
    )
  
  
  # comparison_df has columns Peptide, <group1>_count, <group2>_count
  g1_cnt <- paste0(group1, "_count")
  g2_cnt <- paste0(group2, "_count")
  
  comparison_df %>%
    dplyr::filter(.data[[group1]] >= prevalence_threshold |
                    .data[[group2]] >= prevalence_threshold) 
  
  set1 <- comparison_df %>%
    filter(!!sym(g1_cnt) >  0) %>%   # present in group1
    pull(Peptide) %>% 
    unique()
  
  set2 <- comparison_df %>%
    filter(!!sym(g2_cnt) >  0) %>%   # present in group2
    pull(Peptide) %>% 
    unique()
  
  
  # plot Venn with ggvenn
  venn_input <- setNames(list(set1, set2), c(group1, group2))
  p3 <- ggvenn(venn_input,
         fill_color   = c(custom_colors[group1][[1]], custom_colors[group2][[1]]),
         stroke_size  = 0.1,
         set_name_size= 3,
         text_size    = 3, 
         auto_scale = F,
         show_outside = "none") +
    theme(
      #aspect.ratio = 1,                   # keep the Venn itself square
      plot.margin  = margin(0,0,0,0)   # top = 15px, right/bottom/left = 2px
    )
  
  #combined <- (p1  | ( p3 / p2 + plot_layout(heights = c(0.95)) ))
  combined <- (p3 + p1 + p2) +
    plot_layout(
      ncol   = 3,
      widths = c(1.4, 1.7, 1.7),
      #align  = "v"
    )
  return(combined)
}


####################################################
################## MSA plot ########################
####################################################

make_flagged_msa_plot <- function(df,
                                  flag,
                                  log2_ratio_cut   = 1,
                                  pval_cut         = 0.05,
                                  seq_length       = 64,
                                  coords           = c(NULL, NULL),   # e.g. c(18, 71)
                                  msa_method       = "Muscle",
                                  msa_type         = "protein",
                                  msa_font         = "helvetic",
                                  msa_color_scheme = "Chemistry_AA",
                                  msa_char_width   = 0.5
) {
  
  # Filter
  df_hits <- df %>%
    mutate(log2ratio = log2(ratio)) %>%
    filter(
      !!sym(flag) == TRUE,
      (log2_ratio_cut >= 0 & log2ratio >  log2_ratio_cut) |
        (log2_ratio_cut <  0 & log2ratio <  log2_ratio_cut),
      pvals_not_adj <= pval_cut,
      nchar(aa_seq) >= seq_length
    )
  
  if (nrow(df_hits) < 2) {
    stop("Not enough peptides passed the filters.")
  }
  
  # Build AAStringSet
  seqs <- setNames(df_hits$aa_seq, df_hits$Peptide)
  aa_set <- Biostrings::AAStringSet(seqs, use.names = TRUE)
  
  # MSA
  aln <- msa::msa(aa_set,
                  method = msa_method,
                  type   = msa_type,
                  order  = "input")
  class(aln) <- "AAMultipleAlignment"
  
  # Plot
  p <- ggmsa(
    aln,
    start = coords[1], end = coords[2],
    color    = msa_color_scheme,
    font = msa_font,
    char_width = msa_char_width, 
    seq_name = TRUE,
    consensus_views = TRUE,
    disagreement = FALSE,
    ignore_gaps = FALSE
    ) +
    geom_seqlogo(color = msa_color_scheme, font = msa_font, adaptive = TRUE) +
    geom_msaBar() 
  
  return (p)
}
