source("../../phipseq_report-generator/template/utils/report_utils.R")


# Function to generate count tables per taxonomic rank
generate_taxa_count_tables <- function(meta, exist_df) {
  meta <- meta %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "peptide")
  
  exist_df <- exist_df %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "peptide")
  
  ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  common_peptides <- intersect(meta$peptide, exist_df$peptide)
  
  meta <- meta %>% filter(peptide %in% common_peptides)
  exist_df <- exist_df %>% filter(peptide %in% common_peptides)
  
  taxa_counts <- list()
  for (rank in ranks) {
    # filter out missing annotations
    df <- meta %>%
      select(peptide, !!sym(rank)) %>%
      filter(!is.na(!!sym(rank))) %>%
      inner_join(exist_df, by = "peptide") %>%
      group_by_at(rank) %>%
      summarise(across(-peptide, sum, .names = "{.col}"), .groups = "drop") %>%
      tibble::column_to_rownames(var = rank)
    
    taxa_counts[[rank]] <- df
  }
  
  return(taxa_counts)
}


compute_alpha_diversity <- function(ct, sample_meta, group_col = "Group") {
  # Ensure samples x taxa
  mat <- if (all(colnames(ct) %in% sample_meta$SampleName)) t(ct) else ct
  # Compute metrics
  metrics <- data.frame(
    SampleName   = rownames(mat),
    Richness = vegan::specnumber(mat),
    "Shannon Diversity"  = vegan::diversity(mat, index = "shannon"),
    "Simpson Diversity"  = vegan::diversity(mat, index = "simpson"),
    check.names       = FALSE
  )
  # Merge metadata
  metrics <- merge(metrics, sample_meta, by.x = "SampleName", by.y = "SampleName",     check.names       = FALSE)
  # Rename grouping column to 'Group' for plotting consistency
  colnames(metrics)[which(names(metrics) == group_col)] <- "Group"
  return(metrics)
}

normalize_taxa_counts <- function(ct, norm_method = c("relative","hellinger","log","none")){
  norm_method <- match.arg(norm_method)
  if (norm_method == "relative") {
    ct <- sweep(ct, 1, rowSums(ct), "/")
  } else if (norm_method == "hellinger") {
    ct <- decostand(ct, "hellinger")
  } else if (norm_method == "log") {
    ct <- log1p(ct)
  } 
  return(ct)
}

prep_abund_data <- function(ct, sample_meta, group_col = "Group", norm_method = c("relative","hellinger","log","none")) {
  
  # ensure samples × taxa and align to metadata
  mat <- if (all(colnames(ct) %in% sample_meta$SampleName)) t(ct) else ct
  
  sample_meta <- sample_meta %>%
    rename(Group = !!sym(group_col)) %>%
    filter(!is.na(Group))
  
  common <- intersect(rownames(mat), sample_meta$SampleName)
  mat         <- mat[common, , drop = FALSE]

  #. normalize
  mat <- normalize_taxa_counts(mat, norm_method)
  # norm_method <- match.arg(norm_method)
  # if (norm_method == "relative") {
  #   mat <- sweep(mat, 1, rowSums(mat), "/")
  # } else if (norm_method == "hellinger") {
  #   mat <- decostand(mat, "hellinger")
  # } else if (norm_method == "log") {
  #   mat <- log1p(mat)
  # }
  # # else norm_method == "none": leave raw counts
  
  list(mat = mat, sample_meta = sample_meta)
}

compute_beta_diversity <- function(ct, sample_meta, group_col = "Group", method = "bray", 
                                   norm_method = c("relative","hellinger","log", "none"), permutations = 1000) {
  
  prep  <- prep_abund_data(ct, sample_meta, group_col, norm_method)
  mat   <- prep$mat
  sample_meta  <- prep$sample_meta
  
  dist_mat <- vegan::vegdist(mat, method = method)
  pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)
  var_exp_all <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig), 2)
  
  perm <- vegan::adonis2(dist_mat ~ Group, data = sample_meta, permutations = permutations)
  p_value <- perm$`Pr(>F)`[1]
  
  disp <- vegan::betadisper(dist_mat, sample_meta$Group)

  scores <- data.frame(
    SampleName = rownames(pcoa$points),
    PC1    = pcoa$points[,1],
    PC2    = pcoa$points[,2]
  )
  
  scores <- merge(scores, sample_meta, by.x = "SampleName", by.y = "SampleName")
  
  return(list(dist = dist_mat, permanova = perm, pcoa = scores, beta_dispersion = disp, exp_variance = var_exp_all))
}


# nice_p <- function(p, cutoff = 0.001, digits = 3) {
#   if (is.na(p)) return(NA_character_)
#   if (p < cutoff) {
#     paste0("<", format(cutoff, scientific = FALSE, trim = TRUE))
#   } else {
#     sprintf(paste0("%.", digits, "f"), p)
#   }
# }


compute_pca <- function(ct, sample_meta,
                                  group_col   = "Group",
                                  norm_method = c("relative","hellinger","log","none")) {
  
  prep  <- prep_abund_data(ct, sample_meta, group_col, norm_method)
  mat   <- prep$mat
  samp  <- prep$sample_meta
  
  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 2)
  
  scores <- data.frame(
    SampleName = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2]
  ) %>% merge(samp, by = "SampleName")
  
  
  list(
       scores = scores,
       exp_variance   = var_exp
       )
}


# Given:
#  • rank            — e.g. "species" (must exist in counts_list)
#  • group_col       — grouping variable name in metadata
#  • counts_list     — named list of count matrices
#  • metadata        — data.frame with SampleName + all group_cols
#  • abundance_cutoff = 0.01
#  • min_samples      = 5
#  • top_n            = how many species to keep
#  • order_by         = which level of Group to use for ordering (or "All" for overall)
# Returns a list with:
#   df_long,            # full long format used for stats/plot
#   pval_df,            # Dunn’s test results with xmin/xmax/y.position
#   species_colors      # named vector for species colors
prep_taxa_stats <- function(
    rank,
    group_col,
    counts_list,
    metadata,
    drop_taxa = T,
    abundance_cutoff = 0.01,
    min_samples      = 5,
    top_n            = 12,
    order_by         = "All"    # "All" or a specific group level
) {
  # pull & clean counts
  counts <- counts_list[[rank]] 
  
  if(drop_taxa & rank  == "species"){
    counts <- counts %>%
      tibble::rownames_to_column("Species") %>%
      dplyr::filter(Species != "Homo sapiens") %>% #let's ignore homo sapiens for now
      tibble::column_to_rownames("Species")
  }
  
  metadata  <- metadata %>% 
    select(SampleName=SampleName, group_col) %>% rename(Group = !!sym(group_col))  
  
  # relative abundance & filter low-prevalence
  rel <- normalize_taxa_counts(t(counts), "relative")
  keep <- colSums(rel > abundance_cutoff) >= min_samples
  rel <- rel[, keep, drop = FALSE]
  
  # pivot into long & join metadata
  df_long <- dplyr::as_tibble(rel, rownames = "SampleName") %>%
    tidyr::pivot_longer(-SampleName, names_to = "Rank", values_to = "Abundance") %>%
    dplyr::left_join(metadata, by = "SampleName") %>%
    dplyr::filter(!is.na(Group)) %>%
    dplyr::mutate(Abundance = tidyr::replace_na(Abundance, 0))
  
  # Kruskal–Wallis per taxa
  # kw_res <- df_long %>%
  #   group_by(Rank) %>%
  #   kruskal_test(Abundance ~ Group) %>%
  #   adjust_pvalue(method = "BH") %>%
  #   filter(p.adj <= 0.1)  # keep only taxa with significant omnibus
  # 
  #  Dunn post-hoc only on those taxa
  
  # Dunn’s test + position
  pval_df <- df_long %>%
    #filter(Species %in% kw_res$Species) %>%
    dplyr::group_by(Rank) %>%
    rstatix::dunn_test(Abundance ~ Group, p.adjust.method = "BH") %>%
    rstatix::add_xy_position(x = "Group") %>%
    dplyr::filter(p.adj < 0.05) %>%
    dplyr::mutate(p.adj.label = signif(p.adj, 2)) %>%
    dplyr::ungroup()
  
  # pick top_n Rank
  unique_sp <- unique(pval_df$Rank)
  n_keep    <- min(top_n, length(unique_sp))
  top_Rank <- unique_sp[seq_len(n_keep)]
  
  # filter df_long to these Rank
  df_long_filtered <- df_long %>%
    dplyr::filter(Rank %in% top_Rank)
  
  # determine ordering
  ordering_data <- if (order_by != "All") {
    df_long_filtered %>% dplyr::filter(Group == order_by)
  } else {
    df_long_filtered
  }
  
  ordering <- ordering_data %>%
    dplyr::group_by(Rank) %>%
    dplyr::summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(mean_ab)) %>%
    dplyr::pull(Rank)
  
  df_long_filtered <- df_long_filtered %>%
    dplyr::mutate(Rank = factor(Rank, levels = ordering))
  
  pval_df <- pval_df %>%
    dplyr::filter(Rank %in% ordering) %>% 
    dplyr::mutate(Rank = factor(Rank, levels = ordering))
  
  # assign colors
  
  Rank_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(df_long_filtered$Rank)))
  names(Rank_colors) <- unique(df_long_filtered$Rank)  # Map names to colors
  
  
  
  list(
    rel_abundance        = df_long, 
    rel_abundance_filtered = df_long_filtered,
    pval_df              = pval_df,
    rank_colors       = Rank_colors
  )
}


## plots

plot_pcoa <- function(beta_diversity, custom_colors, 
                      ellipse = T, permutations = 999, show_legend = T) {
  
  mds_scores <- beta_diversity$pcoa
  permanova  <- beta_diversity$perm
  disp       <- beta_diversity$beta_dispersion
  var_exp_all <- beta_diversity$exp_variance
  
  disp_perm  <- vegan::permutest(disp, permutations = permutations)
  centroids <- disp$centroids %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var="Group")
  
  
  # extract the raw p’s
  p_perm_val <- permanova$`Pr(>F)`[1]
  p_disp_val <- disp_perm$tab$`Pr(>F)`[1]
  
  # build your labels
  lab_perm <- paste0("PERMANOVA p = ", format.pval(p_perm_val, digits = 1, eps = 0.01))
  lab_disp <- paste0( "Dispersion   p = ", format.pval(p_disp_val, digits = 1, eps = 0.01))
  
  
  p <- ggplot(mds_scores, aes(x = PC1, y = PC2, fill = Group))
  
  # Conditionally add the ellipse
  if (ellipse) {
    p <- p + stat_ellipse(aes(colour = Group), type = "t", level = 0.95, fill = NA, geom = "path", size = 0.8, alpha = 0.5, show.legend = FALSE)
  }
  
  # Add all the remaining layers, making sure to use + at the end of each line
  p <- p +  geom_point(size = 3, alpha = 0.5, shape = 21, color = "black",  stroke = 0, aes(text = paste("Sample:", SampleName))) +
    geom_point(data = centroids, aes(x = PCoA1, y = PCoA2, fill = Group),
               show.legend = FALSE, size = 5, shape = 21, color = "black") +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    guides(color = "none") +
    labs(
      x = paste0("MDS 1 (", var_exp_all[1], "%)"),
      y = paste0("MDS 2 (", var_exp_all[2], "%)")
      #,title = glue::glue(
      #   "PERMANOVA p={format(permanova$`Pr(>F)`[1], digits=2)}; dispersion p={format(disp_perm$tab$`Pr(>F)`[1],digits=2)}"
      # )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    annotate(
      "text",
      x = Inf,           # right edge
      y = -Inf,          # bottom edge
      label = paste(lab_perm, lab_disp, sep = "\n"),
      hjust = 1.1,       # nudge left a bit
      vjust = -0.1,      # nudge up a bit
      size = 3           # adjust to taste
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}



# Given the output of prep_Rank_stats(), makes the ggplot + annotations
plot_taxa_violin <- function(prep_out, ncol = 3, nudge_y = -1) {
  p <- ggplot(prep_out$rel_abundance_filtered, aes(x = Group, y = Abundance)) +
    geom_violin(fill = "gray80") +
    geom_jitter(aes(color = Rank), width = 0.2, alpha = 0.6, size = 1) +
    stat_summary(fun = median, geom = "point", color = "black", size = 2) +
    scale_color_manual(values = prep_out$rank_colors) +
    #scale_y_log10(expand = expansion(mult = c(0, 0.11))) +
    scale_y_continuous(
      trans  = "log10", 
      #labels = percent_format(accuracy = 1),   # “0.01”→“1%”, “0.1”→“10%”
      expand = expansion(mult = c(0, 0.11))  # 5% below, 20% above
    ) +
    facet_wrap(~ Rank, ncol = ncol)+#, scales = "free_y") +
    theme_pubclean(base_size = 10) +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.margin      = margin(t = 10, r = 5, b = 10, l = 10)
    ) +
    labs(
      x = "Group",
      y = expression("log"[10]~"(Relative Abundance)")
    ) 
  p + stat_pvalue_manual(
    data          = prep_out$pval_df,
    label         = "p.adj.label",
    #xmin         = "xmin",
    #xmax         = "xmax",
    y.position    = "y.position",
    tip.length    = 0.02,
    bracket.size  = 0.25,
    size          = 2.5,
    inherit.aes   = FALSE,
    step.increase = 0.1,
    step.group.by = "Rank",
    bracket.nudge.y = nudge_y
  )
}


plot_pca_abundance <- function(pca, custom_colors, 
                   ellipse = T,  show_legend = T) {
  
  scores <- pca$scores
  var_exp_all <- pca$exp_variance
 
  
  p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Group))
  
  # Conditionally add the ellipse
  if (ellipse) {
    p <- p + stat_ellipse(aes(colour = Group), type = "t", level = 0.95, fill = NA, geom = "path", size = 0.8, alpha = 0.5, show.legend = FALSE)
  }
  
  p <- p +  geom_point(size = 3, alpha = 0.5, shape = 21, color = "black",  stroke = 0, aes(text = paste("Sample:", SampleName))) +
    scale_fill_manual(values = custom_colors) +
    scale_colour_manual(values = custom_colors) +
    guides(color = "none") +
    labs(
      x = paste0("PC 1 (", var_exp_all[1], "%)"),
      y = paste0("PC 2 (", var_exp_all[2], "%)")
      #,title = glue::glue(
      #   "PERMANOVA p={format(permanova$`Pr(>F)`[1], digits=2)}; dispersion p={format(disp_perm$tab$`Pr(>F)`[1],digits=2)}"
      # )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) 
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}



plot_beta_dispersion <- function(disp, sample_meta, custom_colors, group_col = "Group") {
  sample_meta <- sample_meta %>% rename(Group = !!sym(group_col))
  
  df <- data.frame(SampleName = names(disp$distances), Distance = disp$distances)
  df <- merge(df, sample_meta, by.x = "SampleName", by.y = "SampleName")
  
  pairwise_comparisons <- combn(levels(factor(sample_meta[["Group"]])), 2, simplify = FALSE)
  
  p <- plot_groups_boxplots(data = df, 
                            group_col = "Group", 
                            values_col = "Distance",
                            custom_colors = custom_colors, 
                            pairwise_comparisons = pairwise_comparisons,
                            label_axis = c("Group", "Distance to centroid"))
  
  return(p)
}
