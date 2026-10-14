# Functions for plotting data
# Created 12/29/2021 by Jovian Yu and Alan Cooper

##### Packages #####
suppressMessages(library(readxl))
suppressMessages(library(tidyverse))
suppressMessages(library(vroom))

# Plot
suppressMessages(library(ggtext))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(ggfortify))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(rstatix))
library(ggnewscale)
library(ggrepel)
library(patchwork)

# Survival
suppressMessages(library(survival))
suppressMessages(library(survminer))



##### Set up #####

# NOTE: some of the functions below use get_file(), which depends on utils.R

# Save pictures more easily (wrapper for ggplot2::ggsave)
save_pic <- function(name = "pic", plot = last_plot(), width = 6.67, height = 6.67, dpi = 600, prefix = NULL, output_folder = "figures", device = "pdf", ...) {
  # cairo_pdf embeds fonts within images. This makes the files larger, but makes sure fonts like Arial always show regardless of the file viewer
 
  # Make string for device to append at end of file
  if (!is.character(device)) {
    device_name <-  deparse(substitute(device))
    
    if (device_name == "cairo_pdf") {
      device_name <- "pdf"
    }
    
  }else{
    device_name <- device
  }
  
  if (is.null(prefix)) {
    file_name <- name %>% paste(device_name, sep = ".")
  }else if (prefix == "") {
    file_name <- name %>% paste(device_name, sep = ".")
  }else{
    file_name <- paste(prefix, name, sep = "_") %>% paste(device_name, sep = ".")
  }
  
  file_path <- get_file(output_folder, file_name) 
  ggsave(file_path, plot, device = device, width = width, height = height, dpi = dpi, units = "in", ...)
  
}

##### Stats #####


fun_fisher <- function(data, group, var,...){
  
  # Make contingency table
  data_tib <- data %>%
    dplyr::filter(
      !is.na(!!as.symbol(var)),
      !is.na(!!as.symbol(group))) %>% 
    # mutate(across(where(is.factor), ~unfactor(.))) %>% # unfactor everything
    dplyr::count(!!as.symbol(group), !!as.symbol(var)) %>%
    complete(!!as.symbol(group), !!as.symbol(var), fill = list(n = 0)) %>% # fill in 0 for values with no counts
    pivot_wider(names_from = all_of(group), values_from = "n") %>% 
    column_to_rownames(var = var) %>% 
    print()
  
  pval_tib <- data_tib %>% 
    as.matrix() %>% 
    fisher.test(alternative = "two.sided",...) %>% # default to two-sided test
    generics::tidy() 
  
  return(pval_tib)
  
}


fun_chisq <- function(data, group, var){
  
  # Make contingency table
  data_tib <- data %>%
    dplyr::filter(
      !is.na(!!as.symbol(var)),
      !is.na(!!as.symbol(group))) %>% 
    # mutate(across(where(is.factor), ~unfactor(.))) %>% # unfactor everything
    dplyr::count(!!as.symbol(group), !!as.symbol(var)) %>%
    complete(!!as.symbol(group), !!as.symbol(var), fill = list(n = 0)) %>% # fill in 0 for values with no counts
    pivot_wider(names_from = all_of(group), values_from = "n") %>% 
    column_to_rownames(var = var) %>% 
    print()
  
  pval_tib <- data_tib %>% 
    as.matrix() %>% 
    chisq.test() %>%
    generics::tidy() 
  
  return(pval_tib)
  
}

fun_fisher_multi <- function(data, group, var, comparisons = NULL, threshold = 0.05, adjust_method = "BH"){
  
  # Make contingency table
  data_mat <- data %>%
    dplyr::filter(
      !is.na(!!as.symbol(var)),
      !is.na(!!as.symbol(group))) %>% 
    mutate(across(where(is.factor), ~unfactor(.))) %>% # unfactor everything
    dplyr::count(!!as.symbol(group), !!as.symbol(var)) %>%
    complete(!!as.symbol(group), !!as.symbol(var), fill = list(n = 0)) %>%  # fill in 0 for values with no counts
    pivot_wider(names_from = all_of(group), values_from = "n") %>% 
    column_to_rownames(var = var) %>% 
    as.matrix() %>% 
    print()
  
  # Multi-way Fisher's exact test
  multi_tib <- data_mat %>% 
    fisher.test(simulate.p.value = TRUE, B = 1e6) %>% # need simulated p-values for multi-way
    generics::tidy()

  # If multi-way test passeses threshold, do post-hoc tests
  if (multi_tib$p.value >= threshold) {
    
    return(list(multi = multi_tib))
    
  }else{
      # Pairwise tests (rxc = rows x columns)
      # These should only be done if multi-way rxc is significant
      # If 2x2, use rstatix::pairwise_fisher_test()
      # If rx2, use do rstatix::row_wise_fisher_test()
      #   - this should only be done if multi-way rx2 is also significant
      # If rxc,  use another mulit-way Fisher's exact test
  
      if (dim(data_mat)[1] == 2 ) {
        message("2x2 pairwise tests (pairwise_fisher_test)")
        
        pairwise_tib <- data_mat %>% 
          rstatix::pairwise_fisher_test(p.adjust.method = adjust_method, detailed = TRUE) %>% 
          # List comparisons in either order, to filter out comparisons that were not requested
          mutate(
            comparison1 = paste(group1, group2),
            comparison2 = paste(group2, group1))
        
        # Filter to comparisons in input
        if (!is.null(comparisons)) {
          
          comparisons_for_filter <- sapply(comparisons, function(x) paste(x, collapse = " "))
          
          pairwise_tib <- pairwise_tib %>% 
            dplyr::filter(
              comparison1 %in% comparisons_for_filter | comparison2 %in% comparisons_for_filter)
        }
      
      }else if (dim(data_mat)[2] == 2) {
        message("rx2 pairwise tests (row_wise_fisher_test)")
        
        pairwise_tib <- data_mat %>% 
          rstatix::row_wise_fisher_test(p.adjust.method = adjust_method, detailed = TRUE)
        
        # Filter to comparisons in input
        if (!is.null(comparisons)) {
          pairwise_tib <- pairwise_tib %>% 
            dplyr::filter(group == comparisons)
        }
          
      }else{
        message("rxc multi-way pairwise tests")
        
        if (is.null(comparisons)) {
          message("No comparisons gives. Using all possible pairwise comparisons")
          
          comparisons_mat <- t(combn(colnames(data_mat), 2))  # all pairwise groups
          comparisons <- lapply(apply(comparisons_mat, 1, list), unlist) # unlist to get list of length-2 vectors
          
        }
        
        pairwise_tib <- bind_rows(lapply(comparisons, function(x){
          comp_data <- data_mat[,unlist(x)]
          print(comp_data)
          comp_data %>% 
            fisher.test(simulate.p.value = TRUE, B = 1e6) %>% 
            generics::tidy() %>% 
            mutate(
              comparison = paste(x[1], x[2]),
              group1 = x[1],
              group2 = x[2]) %>% 
            return()
        }))
        
        pairwise_tib <- pairwise_tib %>% 
          mutate(p.adj = p.adjust(p.value, method = adjust_method))
      }
    
    pval_l <- list(multi = multi_tib, pairwise = pairwise_tib)
    
    return(pval_l)
  }
}

fun_wilcox <- function(data, group, var, continuity_correct = FALSE){
  
  wilcox_tib <- data %>%
    dplyr::filter(!is.na(!!as.symbol(var)), !is.na(!!as.symbol(group))) %>% 
    wilcox.test(formula = as.formula(paste0("`",var,"`", " ~ ", "`", group, "`")), data = .,
                correct = continuity_correct, conf.int = TRUE, paired = FALSE) %>% 
    generics::tidy()
  
  return(wilcox_tib)
  
}


fun_kw <- function(data, group, var, comparisons = NULL, threshold = 0.05, adjust_method = "BH", post_hoc = "dunn", continuity_correct = FALSE){
  # Can do post-hoc Dunn's Test or pairwise Wilcoxon
  
  kw_tib <- data %>%
    dplyr::filter(!is.na(!!as.symbol(var)), !is.na(!!as.symbol(group))) %>% 
    kruskal_test(formula = as.formula(paste0("`",var,"`", " ~ ", "`", group, "`")))
  
  # If Kruskall-Wallis test passes threshold, do post-hoc tests
  if (kw_tib$p >= threshold) {
    return(list(kw = kw_tib))
    
  }else{
    # Post-hoc Dunn's test with BH-adjusted p-values
    
    if (post_hoc == "dunn") {
      message("Performing post-hoc Dunn's test with adjusted p-values")
      
      posthoc_tib <- data %>%
        dplyr::filter(
          !is.na(!!as.symbol(var)), !is.na(!!as.symbol(group))) %>% 
        dunn_test(formula = as.formula(paste0("`",var,"`", " ~ ", "`", group, "`")), p.adjust.method = adjust_method) %>% 
        mutate(
          comparison1 = paste(group1, group2),
          comparison2 = paste(group2, group1))
    
    }else if (post_hoc == "wilcoxon") {
      message("Performing post-hoc Wilcoxon tests with adjusted p-values")
      
      groups <- data %>% dplyr::pull(!!as.symbol(group)) %>% unique()
      if (is.factor(groups)) groups <- unfactor(groups)
      
      comparisons <- as.list(as.data.frame(combn(groups, 2)))
      
      posthoc_tib <- bind_rows(lapply(comparisons, function(x){
        data %>% 
          dplyr::filter(!!as.symbol(group) %in% unlist(x)) %>% 
          fun_wilcox(group = group, var = var, continuity_correct = continuity_correct) %>% 
          mutate(
            comparison1 = paste(x[1], x[2]),
            comparison2 = paste(x[2], x[1]),
            group1 = x[1],
            group2 = x[2]) %>% 
          return()
      })) %>% 
      dplyr::rename(p = p.value) %>% 
      mutate(p.adj = p.adjust(p, method = adjust_method))
      
      
    }
  
    
    # Filter to comparisons in input
    if (!is.null(comparisons)) {
      
      comparisons_for_filter <- sapply(comparisons, function(x) paste(x, collapse = " "))
      
      posthoc_tib <- posthoc_tib %>% 
        dplyr::filter(comparison1 %in% comparisons_for_filter | comparison2 %in% comparisons_for_filter)
    }
    
    if (post_hoc == "dunn"){
      pval_l <- list(kw = kw_tib, dunn = posthoc_tib)  
    }else if (post_hoc == "wilcoxon"){
      pval_l <- list(kw = kw_tib, wilcoxon = posthoc_tib)  
    }
      
    return(pval_l)
  }
}

fun_cor <- function(data, x, y, method = "pearson"){
  
  cor_tib <- data %>%
    dplyr::filter(!is.na(!!as.symbol(x)), !is.na(!!as.symbol(y))) %>% 
    cor_test(
      !!as.symbol(x), !!as.symbol(y), 
      alternative = "two.sided",
      method = method)
  
  return(cor_tib)
  
}

##### Colors and themes #####

TEXTSIZE_DBL <-  34
LABEL_TEXTSIZE_DBL <- TEXTSIZE_DBL*0.3

plot_theme <- theme_classic() + # previously include base_family = "Arial"
  theme(
    axis.line = element_line(linewidth = 1.25, color = "black"),
    axis.ticks = element_line(linewidth = 1.25, color = "black"),
    axis.ticks.length = unit(2, "mm"),
    axis.title = element_text(size = TEXTSIZE_DBL),
    axis.text.x = element_text(size = TEXTSIZE_DBL, color = "black"),
    axis.text.y = element_text(size = TEXTSIZE_DBL, color = "black"),
    axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)), 
    plot.title = element_text(size = TEXTSIZE_DBL),
    strip.text.x = element_text(size = TEXTSIZE_DBL),
    legend.text = element_text(size = TEXTSIZE_DBL, hjust = 0), #TEXTSIZE_DBL-2
    legend.title = element_text(size = TEXTSIZE_DBL, vjust = -0.5, hjust = 0),
    # strip.background = element_blank()
    )

##### Basic plots #####


get_pvaly <- function(df, y_var){
  # Get y-axis value to put p-value at, based on input dataframe
  # Expects var to be column name in df
  
  var_max <- df %>% pull(!!as.symbol(y_var)) %>% max()
  return(min(var_max*1.32, var_max + 22))
  
}


format_pvals <- function(pvals, label = FALSE, add_stars = FALSE, show_ns = TRUE, point_one = FALSE){
  # Format p-values
  # Rounds to two significant digits
  # Adds scientific notation for p < 0.0001
  # Takes numeric as input, returns character
  # Formats as "p = #" if label = TRUE
  
  # extra - whether to use 0.1
  
  pvals <- as.numeric(pvals)
  
  fix_num <- function(pvals){
    sapply(pvals, function(x) format(x, scientific = 0, digits = 2, nsmall = 2)) %>% return()
  }
  
  if (show_ns) {
    ns_symbol <- "ns"
  } else {
    ns_symbol <- ""
  }
  
  pvals2stars <- function(x, point_one) {
    case_when(
    pvals < 0.0001 ~ "****",
    pvals < 0.001 ~ "***",
    pvals < 0.01 ~ "**",
    pvals < 0.05 ~ "*",
    pvals < 0.1 & point_one ~ "<0.1",
    TRUE ~ ns_symbol)
  }
  
  if (label) {
    pvals <- paste0("p = ", fix_num(pvals))
  }else if (add_stars) {
    pvals <- pvals2stars(pvals, point_one)
  }else{
    pvals <- fix_num(pvals)
  }
  
  return(pvals)
  
}

plot_pca <- function(data, var, new_encircle = FALSE, group_to_emphasize = NULL, 
                     color_vec = pal_nejm()(8), title = NULL, color_encircle = FALSE,
                     pca = list_pca_train$pca, full_data = train_clindata_tib, point_size = 4, point_emphasis = 1.75){
  
  # pca takes output from prcomp
  # - list_pca_train$pca is how PCA is used in processCountsMutations.Rmd
  # color_encircle - whether to color encircling lines by cluster
  
  # Variance explained by each PC
  pca_variance <- round(pca$sdev^2/sum(pca$sdev^2)*100, 2)
  
  # Want to keep original outlines defined on all samples
  # If new_encircle = FALSE and full_data is defined, then will make based on that
  # Whether to color encircle black or by cluster is also updated here (either in aes() or not)
  
  encircle_line_size <- 2
  
  if (!new_encircle) {
    
    if (color_encircle) {
      encircle <- geom_encircle(
        data = full_data,
        inherit.aes = FALSE,
        aes(x = hc_axis, y = coo_axis, group = cluster, color = cluster), 
        alpha = 1, 
        expand = 0, 
        s_shape = 0.6,
        linetype = 2, 
        size = encircle_line_size, 
        show.legend = FALSE)
    }else{
      encircle <- geom_encircle(
        data = full_data,
        inherit.aes = FALSE,
        aes(x = hc_axis, y = coo_axis, group = cluster), 
        color = "black",
        alpha = 1, 
        expand = 0, 
        s_shape = 0.6,
        linetype = 2, 
        size = encircle_line_size, 
        show.legend = FALSE)
    
    }
    
  }else{
    if (color_encircle) {
      encircle <- geom_encircle(
        data = data,
        inherit.aes = FALSE,
        aes(x = hc_axis, y = coo_axis, group = cluster, color = cluster), 
        alpha = 1, 
        expand = 0, 
        s_shape = 0.6,
        linetype = 2, 
        size = encircle_line_size, 
        show.legend = FALSE)
      
    }else{
      encircle <- geom_encircle(
        data = data,
        aes(x = hc_axis, y = coo_axis, group = cluster), 
        color = "black",
        alpha = 1, 
        expand = 0, 
        s_shape = 0.6,
        linetype = 2, 
        size = encircle_line_size, 
        show.legend = FALSE)
    }
  }
  
  plot <- data %>% 
    ggplot(aes(x = hc_axis,  y = coo_axis, color = !!as.symbol(var))) +
    geom_hline(yintercept = 0, linewidth = 0.5) + 
    geom_vline(xintercept = 0, linewidth = 0.5) +
    # geom_point(size = 4, alpha = alpha) +
    ggtitle(title) +
    plot_theme +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = TEXTSIZE_DBL-2), # make label of % variance explained smaller
      axis.text.x = element_text(size = TEXTSIZE_DBL-2)) +
    encircle +
    scale_color_manual(values = color_vec) +
    # scale_x_continuous(paste0("PC1 ","(", pca_variance[1], "%", ")"), breaks = 0) +
    # scale_y_continuous(paste0("PC2 ","(", pca_variance[2], "%", ")"), breaks = 0)
    scale_x_continuous(paste0("PC1 ","(", pca_variance[1], "%", ")"), breaks = NULL) +
    scale_y_continuous(paste0("PC2 ","(", pca_variance[2], "%", ")"), breaks = NULL)

  # Option to make one subset stand out (group_to_emphasize)
  if (!is.null(group_to_emphasize)) {
    plot <- plot + 
      # Add all points with alpha 0.5
      geom_point(size = point_size, alpha = 0.5) + 
      scale_color_manual(values = color_vec) +
      # Add subset of points with alpha of 1
      geom_point(
        data = dplyr::filter(data, !!as.symbol(var) %in% group_to_emphasize),
        aes(hc_axis, coo_axis, color = !!as.symbol(var)),
        size = point_size*point_emphasis, alpha = 1)
    
  }else{
    
    plot <- plot + 
      scale_color_manual(values = color_vec) +
      geom_point(size = point_size, alpha = 1)
  }
  
  # plot <- plot + encircle

  return(plot)
}

plot_bar <- function(data, group, var, ycol,
                     color_vec = pal_nejm()(8), title = NULL,
                     xlab = NULL, ylab = ycol,
                     bar_stat = "identity", bar_position = position_dodge2(preserve = "single"), bar_width = 0.5, 
                     y_breaks = waiver(), ylim = NULL,
                     pvals = NULL, pval_xmin = NULL, pval_xmax = NULL, pval_y = NULL, pval_bar_tip_length = 0.03, pval_bar_size = 0.75,
                     pval_angle = 0, pval_vjust = 0){

  # Choose max. y-values based on bar type whether or not adding p-values
  bar_position_name <- deparse(substitute(bar_position))
  
  if (!is.null(pvals)) {
    
    if (is.null(pval_y)) {
      pval_y <- get_pvaly(data, y_var = ycol)
    }
    
    y_max <- pval_y + max(pull(data, !!as.symbol(ycol)))/5.5 # tried 10 to 5
    
  }else{
    y_max <- max(pull(data, !!as.symbol(ycol)))
    y_max <- min(y_max*1.25, y_max + 15) 
  }
  
  if (grepl("stack", bar_position_name)) {

    y_max <- data %>% 
      group_by(!!as.symbol(group)) %>% 
      summarize(sum = sum(!!as.symbol(ycol))) %>% 
      dplyr::arrange(dplyr::desc(sum)) %>% 
      dplyr::slice(1) %>% 
      pull(sum)
    
    if (is.null(pval_y)) {
      pval_y <- min(y_max*1.25, y_max + 15) 
      
    }else{
      y_max <- pval_y*1.15 
    }
  }
  
  if (is.null(ylim)) {
    ylim <- c(0, y_max)
  }
  

  plot <- data %>%
    ungroup() %>%
    ggplot(aes(x = !!as.symbol(group), y = !!as.symbol(ycol), fill = !!as.symbol(var))) +
    geom_bar(stat = bar_stat, position = bar_position, width = bar_width) +
    plot_theme +
    theme(
      axis.text.x = element_text(margin = ggplot2::margin(2, 0, 0, 0, "mm")),
      axis.ticks.x = element_line(linewidth = 0),
      axis.ticks.length.x = unit(0, "mm")) +
    labs(x = xlab, y = ylab) +
    scale_fill_manual(values = color_vec) +
    scale_y_continuous(expand = c(0,0), limits = ylim, breaks = y_breaks) +
    ggtitle(title)

  if (!is.null(pvals)) {
    p_df <- data.frame(pvals, pval_xmin, pval_xmax, pval_y)
    
    if (is.numeric(pvals)) {
      p_df <- p_df %>% 
        rowwise() %>% 
        mutate(pvals = format_pvals(pvals)) %>% 
        ungroup()
    }
    
    plot <- plot +
      geom_signif(
        data = p_df,
        aes(xmin = pval_xmin, xmax = pval_xmax, y_position = pval_y, annotations = pvals),
        inherit.aes = FALSE, manual = TRUE,
        textsize = LABEL_TEXTSIZE_DBL,
        vjust = pval_vjust,
        angle = pval_angle,
        tip_length = pval_bar_tip_length, size = pval_bar_size)
    
  }
  return(plot)
}

plot_box <- function(data, group, var, color_vec = pal_nejm()(8), title = NULL, xlab = NULL, ylab = var,
                     ylim = NULL,  y_breaks = waiver(), box_color = NULL,
                     extra_group = NULL, extra_color_vec = pal_nejm()(8), violin = FALSE,
                     pvals = NULL, pval_xmin = NULL, pval_xmax = NULL, pval_y = NULL,
                     pval_bar_tip_length = 0.03, pval_bar_size = 0.75, outlier_size = 2.5,
                     extra_group_point_size = 3, extra_group_alpha = 1, box_alpha = 1){


  if (is.null(pval_y)) {
    pval_y <- get_pvaly(data, y_var = var)
  }
  
  y_max <- pval_y + max(pull(data, !!as.symbol(var)))/5 
  y_min <- min(pull(data, !!as.symbol(var)))
  
  if (is.null(ylim)) {
    if (y_min > 0) {
      ylim <- c(0, y_max) 
    }else{
      ylim <- c(y_min, y_max) 
    }
  }
  
  # Don't show outliers if coloring points by extra group
  # Outlier points are just small black dots
  if (!is.null(extra_group)) {
    outlier_size <- -1
  }
    
  if (!violin) {
    
    if (is.null(box_color)) {
      box_color <- "black"
    }
    
    plot <- data %>% 
      ggplot(aes(x = !!as.symbol(group), y = !!as.symbol(var), fill = !!as.symbol(group))) +
      geom_boxplot(outlier.size = outlier_size, lwd = 1, color = box_color, alpha = box_alpha) +
      scale_y_continuous(limits = ylim, breaks = y_breaks)
    
  }else{
    
    if (is.null(box_color)) {
      box_color <- "white"
    }
    
    
    plot <- data %>% 
      ggplot(aes(x = !!as.symbol(group), y = !!as.symbol(var), fill = !!as.symbol(group))) +
      geom_violin(lwd = 0.5) +
      geom_boxplot(width = 0.1, outlier.size = -1, color = box_color, alpha = box_alpha, lwd = 1) +
      scale_y_continuous(expand = c(0,0), limits = ylim, breaks = y_breaks)
  }

  plot <-  plot +
    plot_theme +
    # theme(
    #   axis.text.x = element_text(margin = ggplot2::margin(2, 0, 0, 0, "mm")),
    #   axis.ticks.x = element_line(linewidth = 0),
    #   axis.ticks.length.x = unit(0, "mm")
    #   ) +
    guides(fill = "none", alpha = "none") +
    labs(x = xlab, y = ylab) +
    ggtitle(title)
  
  if (!is.null(color_vec)) {
    plot <- plot + scale_fill_manual(values = color_vec)
  }

  if (!is.null(pvals)) {
    # Group is to fix bug where duplicate p-values show as one
    p_df <- data.frame(pvals, pval_xmin, pval_xmax, pval_y, group = 1:length(pvals))
    
    if (is.numeric(pvals)) {
      p_df <- p_df %>% 
        rowwise() %>% 
        mutate(pvals = format_pvals(pvals)) %>% 
        ungroup()
    }
    
    plot <- plot +
      geom_signif(
        data = p_df,
        aes(xmin = pval_xmin, xmax = pval_xmax, y_position = pval_y, annotations = pvals, group = group),
        inherit.aes = FALSE, manual = TRUE, textsize = LABEL_TEXTSIZE_DBL, pval_bar_tip_length = pval_bar_tip_length, size = pval_bar_size)
  }
    
  if (!is.null(extra_group)) {
    plot <- plot +
      geom_jitter(
        data = data %>% dplyr::filter(!is.na(!!as.symbol(extra_group))),
        aes(color = !!as.symbol(extra_group)), alpha = extra_group_alpha, size = extra_group_point_size) +
      scale_color_manual(values = extra_color_vec) +
      # scale_alpha_manual(values = rep(0, length(group_levels))) +
      guides(fill = "none", alpha = "none")
  }
  
  return(plot)
}


plot_scatter <- function(data, x, y, group = NULL, color_vec = pal_nejm()(8), title = NULL,
                         xlim = NULL, ylim = NULL, point_size = 2, line_size = 2.5, line_color = "red",
                         xlab = x, ylab = y, cor_method = "pearson", line_ci = FALSE, y_breaks = waiver(), x_breaks = waiver(),
                         cor_label = NULL){
  
  # line_ci is for confidence interval for linear fit

  if (is.null(names(color_vec)) & !is.null(group)) {
    group_levels <- data %>% 
      pull(!!as.symbol(group)) %>% 
      unique()
    
    group_levels <- group_levels[!is.na(group_levels)] %>% sort()
  }else{
    group_levels <- names(color_vec)
  }
  
  plot_data <- data %>% 
    dplyr::filter(!is.na(!!as.symbol(x)), !is.na(!!as.symbol(y)))
  
  
  # Correlation
  if (is.null(cor_label)) {
    cor_tib <- fun_cor(data, x, y, method = cor_method)
    if (cor_method == "pearson") {
      cor_label <- bquote(~italic("r")*.(paste0(" = ", cor_tib$cor,", p = ", format_pvals(cor_tib$p))))
    }
  }
  

  # y limits - need to add offset to fit correlation label without hitting data
  if (is.null(ylim)) {
    y_max <- data %>%
      dplyr::arrange(dplyr::desc(!!as.symbol(y))) %>%
      pull(!!as.symbol(y)) %>% 
      max()
    
    y_min <- data %>%
      group_by(!!as.symbol(y)) %>%
      pull(!!as.symbol(y)) %>% 
      min()
    
    ylim <- c(y_min, y_max*1.1)
  }

  
  # x limits
  if (is.null(xlim)) {
    x_max <- data %>%
      dplyr::arrange(dplyr::desc(!!as.symbol(x))) %>%
      pull(!!as.symbol(x)) %>% 
      max()
    
    x_min <- data %>%
      group_by(!!as.symbol(x)) %>%
      pull(!!as.symbol(x)) %>% 
      min()
    
    xlim <- c(x_min, x_max*1.05)
  }

  plot <- plot_data %>% 
    ggplot(aes(x = !!as.symbol(x), y = !!as.symbol(y))) +
    geom_point(size = point_size) +
    # Linear fit
    geom_smooth(method = lm, se = line_ci, linewidth = line_size, color = line_color, formula = "y~x")  +
    # Add correlation + p-value
    annotate(
      geom = "text", size = LABEL_TEXTSIZE_DBL,
      label = cor_label, x = -Inf, y = Inf, hjust = -0.1, vjust = 1,
    ) +
    plot_theme +
    labs(x = xlab, y = ylab) +
    scale_fill_manual(values = color_vec) +
    scale_alpha_manual(values = rep(1, length(group_levels))) +
    scale_y_continuous(limits = ylim, breaks = y_breaks) +
    scale_x_continuous(limits = xlim, breaks = x_breaks) +
    ggtitle(title)

  
  # For coloring dots based on group
  if (!is.null(group)) {
    plot <- plot_data %>% 
      ggplot(aes(x = !!as.symbol(x), y = !!as.symbol(y))) +
      geom_point(aes(color = !!as.symbol(group)), size = point_size) +
      # Linear fit
      geom_smooth(method = lm, se = line_ci, linewidth = line_size, color = line_color, formula = "y~x")  +
      # Add correlation + p-value
      annotate(
        geom = "text", size = LABEL_TEXTSIZE_DBL,
        label = cor_label, x = -Inf, y = Inf, hjust = -0.1, vjust = 1,
      ) +
      plot_theme +
      labs(x = xlab, y = ylab) +
      scale_fill_manual(values = color_vec) +
      scale_alpha_manual(values = rep(1, length(group_levels))) +
      scale_y_continuous(limits = ylim, breaks = y_breaks) +
      scale_x_continuous(limits = xlim, breaks = x_breaks) +
      scale_color_manual(values = color_vec)
      ggtitle(title)
  }
  
  return(plot)
}

plot_alluvial <- function(data, width, color_vec){
  
  n_axes <- ncol(data)-1
  
  alluvial_tib <- to_lodes_form(data, key = "Variable", value = "Group", id = "Cohort", axes = 1:n_axes)
  
  
  plot <- ggplot(alluvial_tib, aes(x = Variable, y = total, stratum = Group, alluvium = Cohort)) +
    geom_alluvium(aes(fill = Group), width = 2/5) +
    geom_stratum(aes(fill = Group), width = 2/5) +
    geom_richtext(
      stat = "stratum", 
      fill = NA, text.color = "white",
      vjust = 0.6,
      aes(label = after_stat(stratum)), # postpone aes mapping to after transformation
      size = LABEL_TEXTSIZE_DBL,
      label.color = NA,
    ) +
    scale_fill_manual(values = color_vec) +
    theme_bw() + # previously had base_family = "Arial"
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()) +
    guides(fill = "none") 
  
  return(plot)
}

plot_confusion_matrix <- function(data, true_var, pred_var, true_lab = "True class", acc_text_x = 1.5, calc_acc = TRUE, pred_lab = "Predicted class", matching_classes = TRUE, norm = "none"){
  # Expects input from dplyr::count
  
  # norm is one of "row" "col" ; default is none
  
  
  
  raw_data <- data # data that is not normalized
  
  if (norm != "none") {
    if (norm == "row") {
      data <- data %>% 
        group_by(!!as.symbol(true_var)) %>% 
        mutate(prop = 100*n/sum(n)) %>% 
        ungroup()
    }else if (norm == "col") {
      data <- data %>% 
        group_by(!!as.symbol(pred_var)) %>% 
        mutate(prop = 100*n/sum(n)) %>% 
        ungroup()
    }else if (norm == "all") {
      data <- data %>% 
        mutate(prop = 100*n/sum(n)) 
    }
    
  }


  # Plot
  # Make label differently depending on normalization
  # - n for no normalization
  # - n (%) for normalization
  if (norm == "none") {
    plot <- data %>% 
      ggplot(aes(x = !!as.symbol(pred_var), y = !!as.symbol(true_var), fill = n)) +
      geom_tile(color = "black", size = 2) +
      geom_text(aes(label = n), color = "white", size = 10) 
    
  }else{
    plot <- data %>% 
      ggplot(aes(x = !!as.symbol(pred_var), y = !!as.symbol(true_var), fill = prop)) +
      geom_tile(color = "black", linewidth = 2) +
      geom_text(aes(label = paste0(n, " (", round(prop), "%)")), color = "white", size = 10)
  }
  
  
  plot <- plot +
    plot_theme +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0, color = "black", size = plot_theme$axis.text.y$size),
      axis.text.y = element_text(color = "black", size = plot_theme$axis.text.x$size),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.key.height = unit(1.5, "cm"),
      plot.margin = ggplot2::margin(1, 1, 20, 1),
    ) +
    scale_x_discrete(position = "top") +
    labs(x = pred_lab, y = true_lab)
  
  # Accuracy
  # Whether to include only overlapping classes in accuracy
  if (calc_acc) {
    if (matching_classes) {
      
      true_classes <- data %>% pull(!!as.symbol(true_var)) %>% unique()
      pred_classes <- data %>% pull(!!as.symbol(pred_var)) %>% unique()
      
      common_classes <- intersect(true_classes, pred_classes)
      
      acc <- raw_data %>% 
        dplyr::filter(!!as.symbol(pred_var) %in% common_classes, !!as.symbol(true_var) %in% common_classes) %>% 
        mutate(
          !pred_var := factor(!!as.symbol(pred_var)),
          !!true_class := factor(!!as.symbol(true_var))) %>% 
        mutate(correct = !!as.symbol(true_var) == !!as.symbol(pred_var)) %>%
        summarize(acc = sum(correct*n)/sum(n)) %>% 
        pull(acc)
      
    }else{
      acc <- raw_data %>% 
        mutate(across(where(is.factor), ~unfactor(.))) %>%
        mutate(correct = !!as.symbol(true_var) == !!as.symbol(pred_var)) %>%
        summarize(acc = sum(correct*n)/sum(n)) %>% 
        pull(acc)
    }
    
    # Add accuracy to plot
    plot +
      annotate(
        "text", x = acc_text_x, y = 0.3,
        size = LABEL_TEXTSIZE_DBL,
        label = paste0("Accuracy = ", round(acc,2)*100, "%"),
      )  +
      coord_cartesian(clip = "off") %>% 
      return()
  }else{
    return(plot)
  }
  
}


plot_pie <- function(data, group, value, label_color = "black", width = 1, label_values = TRUE, color_vec = pal_nejm()(8), manual=FALSE){
  
  # Adding labels for fractions
  # - label_values = TRUE or 1 : add %
  # - label_values = TRUE or 1 : add % and #
  # - manual = FALSE for doing dplyr::count() on input
  
  if (!manual) {
    data <-  data %>% 
      dplyr::count(!!as.symbol(group)) %>% 
      ungroup() %>% 
      mutate(prop = 100*n/sum(n)) %>% 
      ungroup()
  }

  plot <- data %>% 
    ggplot(aes(x = "", y = prop, fill = !!as.symbol(group))) +
    geom_bar(stat = "identity", width = width, color = label_color) +
    coord_polar("y", start = 0) +
    plot_theme +
    theme(
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()) +
    scale_fill_manual(values = color_vec)
  

  if (label_values == TRUE | label_values == 1) {
    plot <- plot +
      geom_text(aes(label = paste0(round(prop, 1),"%"),1.7),
                position = position_stack(vjust = 0.5),
                size = TEXTSIZE_DBL*0.25)
  }else if (label_values == 2){
    plot <- plot +
      geom_text(aes(label = paste0(round(prop, 1),"%", " (",n,")"),1.85),
                position = position_stack(vjust = 0.5),
                size = TEXTSIZE_DBL*0.25)
  }
  
  return(plot)
}


##### Survival Curves ####

fun_surv_fit <- function(data, group, time_col, status_col, end_time = 5, return_summary = FALSE, p_adjust = "none") {
  
  # Assumes 0 is no event (censored) and 1 is an event (not censored). An event includes death, progression, etc.
  # Assumes time is in years but it can be flexible
  survival <- str_remove(status_col, "_.+")
  
  # Make survival data stop at `end_time` for statistical tests
  # This makes overall p-value reflect time for `end_time` rather than all available time
  # StackExchange: https://stats.stackexchange.com/a/377339
  surv_df <- data %>%
    dplyr::filter(!is.na(!!as.symbol(time_col)), !is.na(!!as.symbol(status_col))) %>%
    mutate(new_status_col = !!as.symbol(status_col)) %>%
    # Set any survival greater than the max. time to be that max. value
    mutate(new_time_col = ifelse(!!as.symbol(time_col) > end_time, end_time, !!as.symbol(time_col))) %>% 
    # Set any status after the years to 0, no event
    mutate(new_status_col := ifelse(!!as.symbol(time_col) > end_time, 0, !!as.symbol(status_col)))
  
  
  # Makes formula of Surv(years, status) ~ group
  surv_fomrula <- as.formula(paste0("Surv(new_time_col, new_status_col) ~ ", group))
  
  fit <- surv_fit(formula = surv_fomrula, data = surv_df)
  
  # Overall log-rank test
  message("Overall log-rank test")
  logrank_df <- surv_pvalue(fit, method = "survdiff")
  print(logrank_df)
  
  # Pairwise survival differences
  # Uses log-rank test (rho = 0)
  message("Pairwise comparisons using log-rank test are without p-value adjustment")
  pairwise_diff <- pairwise_survdiff(surv_fomrula, data = surv_df, p.adjust.method = p_adjust, rho = 0)
  
  # Median survival differences
  message("Median survival with 95% CI")
  surv_med <- surv_median(fit)

  
  # Survival summary - returns estimate of survival in column `surv`
  if (return_summary) {
    surv_sum <- surv_summary(fit)
  }else{
    surv_sum <- NULL
  }
  
  # Survival fit for all available years
  # This is used by plot
  surv_df <- data %>%
    dplyr::filter(!is.na(!!as.symbol(time_col)), !is.na(!!as.symbol(status_col)))
  
  surv_fomrula <- as.formula(paste0("Surv(",time_col, ",", status_col, ") ~ ", group))
  
  fit <- surv_fit(formula = surv_fomrula, data = surv_df)
  
  # Output
  surv_l <- list(
    fit = fit, logrank = logrank_df,
    pairwise = pairwise_diff,
    survival_summary = surv_sum,
    survival_median = surv_med)
  
  return(surv_l)
}

plot_survival <- function(
    data, group, time_col, status_col,
    end_time = 10, break_time_by = 1, 
    pval_x = end_time*0.1, pval_y = 0.25,
    title = "", legend_title = "",
    ylab = NULL, xlab = "Time in years",
    color_vec = pal_nejm()(8),
    custom_group_labels = NULL,
    risk_table = TRUE,
    curve_size = 1.2,
    censor_size = 5.5,
    add_num_to_name = FALSE,
    legend_pos = "right", table_height = 0.25,include_legend = TRUE,
    comparisons = NULL, pvals = TRUE, p_adjust = "none", ...){
  
  
  # if `title`  and `legend_title` are NULL then `ggsurvplot` uses "strata"

  if (add_num_to_name) {
    custom_group_labels <- data %>% 
      dplyr::count(!!as.symbol(group)) %>% 
      mutate(!!group := paste0(!!as.symbol(group), " (n = ", n, ")")) %>% 
      pull(group)
  }
  
  surv_l <- fun_surv_fit(data = data, group = group, time_col = time_col, status_col = status_col, end_time = end_time, p_adjust = p_adjust)
  fit <- surv_l$fit
  pairwise_diff <-  surv_l$pairwise
  

  if (pvals) {
    pvals <- format_pvals(surv_l$logrank$pval, label = TRUE)
    
    # Show pairwise comparisons on plot
    if (!is.null(comparisons)) {
      
      # Build comparisons into data frame
      comparisons <- sapply(comparisons, function(x) paste(sort(x), collapse = " ")) # sort and make into text
      comparison_df <- data.frame(comparison = unname(comparisons), name = names(comparisons))
      
      # Filter pairwise tests into those in comparison
      pvals <- pairwise_diff %>% 
        tidy() %>% 
        rowwise() %>% 
        mutate(comparison = paste(sort(c(group1, group2)), collapse = " ")) %>% 
        dplyr::filter(comparison %in% comparisons) %>% 
        left_join(comparison_df, by = "comparison") %>% 
        arrange(factor(name, levels = names(comparisons))) %>% 
        rowwise() %>% 
        mutate(pval_text = paste0(name, ": ", format_pvals(p.value, label = TRUE))) %>% 
        pull(pval_text) %>% 
        paste(collapse = "\n")
      
    }
  }

  
  # Labels and axes
  group_labels <- data %>% 
    pull(!!as.symbol(group)) %>% 
    unique()
  
  group_labels <- group_labels[!is.na(group_labels)] %>% sort()
  
  if (is.null(ylab)) {
    ylab <- str_remove(status_col, "_.+")
  }

  # Change table theme
  # Remove x-axis and make text smaller
  table_theme <- plot_theme + 
    theme(
      plot.title = element_text(size = 0.8*TEXTSIZE_DBL),
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank())
  
  # Plot
  surv_plot <- ggsurvplot(
    fit,
    pval = pvals,
    ylab = ylab,
    pval.method = FALSE,
    # pval.method.coord = c(pval_x, pval_y),
    pval.coord = c(pval_x, pval_y - 0.1),
    size = curve_size,  #1.2
    censor.size = censor_size,
    ggtheme = plot_theme,
    risk.table.y.text = FALSE,
    risk.table = risk_table,
    tables.theme = table_theme,
    tables.height = table_height,
    fontsize = 0.8*LABEL_TEXTSIZE_DBL,
    pval.size = 0.8*LABEL_TEXTSIZE_DBL,
    pval.method.size = 0.8*LABEL_TEXTSIZE_DBL,
    title = title,
    legend.title = legend_title,
    legend = legend_pos,
    palette = color_vec,
    legend.labs = group_labels,
    tables.y.text = FALSE,
    xlim = c(0, end_time),
    break.time.by = break_time_by, 
    xlab = xlab,...)
  
  # Include/exclude legend
  # Make legend icons bigger
  if (include_legend) {
    if (is.character(legend_pos)) {
      surv_plot$plot <- surv_plot$plot +
        plot_theme +
        xlab(xlab) + 
        guides(color = guide_legend(override.aes = list(size = 7), label.position = legend_pos)) +
        theme(legend.text = element_markdown())
    }else{
      surv_plot$plot <- surv_plot$plot +
        plot_theme +
        xlab(xlab) + 
        guides(color = guide_legend(override.aes = list(size = 7))) +
        theme(legend.position = legend_pos, legend.text = element_markdown())
    }
  }else{
      surv_plot$plot <- surv_plot$plot + guides(color = "none")
  }
  
  # Add new labels for groups (e.g. for adding in italics)
  if (!is.null(custom_group_labels)) {
    surv_plot$plot <- surv_plot$plot +
      scale_color_manual(labels = custom_group_labels, values = color_vec) + 
      theme(legend.text = element_markdown())
  }
  
  # Add survival curves and risk table plots together so they are on ggplot object
  layout <- 
    "A
     B"
  
  if (risk_table) { 
    plot <- surv_plot$plot + surv_plot$table + plot_layout(design = layout, heights = c(2,1))   
  }else{
    plot <- surv_plot$plot
  }

  
  return(list(plot = plot, diff = pairwise_diff))
}



##### DGE and GSEA #####

plot_dge <- function(dge_toptable, genes_to_label = NULL,
                     box_label = TRUE, img_file = NULL, title = NULL,force_label = FALSE,
                     xlim = NULL, ylim = NULL, y_breaks = waiver(), x_breaks = waiver(),
                     log2fc_threshold = 1.5, pval_threshold = 0.05, n_genes = 50,
                     nudge_x = 0, nudge_y = 0, type = "limma", pval_col_adjusted = TRUE,
                     label_line_size = 1.25,
                     point_size = 1.25,
                     max_overlaps = 20,
                     extra_genes = NULL, extra_label_color = "red",
                     group_label_x = NULL, group_label_y = NULL, group1 = "High", group2 = "Low",
                     pathway_genes = NULL, left_color = "red", right_color = "red"){
  # Para,s
  # force_label - label genes that are not significant
  # extra_genes - genes to label in a different color by chagning text color
  
  # Expects data from limma::topTable()
  if (type == "limma") {
    logfc_col <- "logFC"
    
    if (pval_col_adjusted) {
      pval_col <- "adj.P.Val"
    }else{
      pval_col <- "P.Value"
    }
    
  # Expects data from DESeq2::results()
  } else if (type == "DESeq2") {
    logfc_col <- "log2FoldChange"
    
    if (pval_col_adjusted) {
      pval_col <- "padj"
    }else{
      pval_col <- "pvalue"
    }
  }
  
  if (pval_col_adjusted) {
    pval_label <- "-log<sub>10</sub>(adj. p-value)"
    # pval_label <- bquote(~-log[10]*("adj. p-value"))
  }else{
    pval_label <- "$-$log<sub>10</sub>(p-value)"
    # pval_label <- bquote(~-log[10]*("p-value"))
  }
  
  # Reformat data for plotting
  # Define genes to label based on logFC and p-value cutoff
  plot_data <- dge_toptable %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    dplyr::rename(logFC = !!as.symbol(logfc_col)) %>% 
    mutate(
      signif = (!!as.symbol(pval_col) <= pval_threshold) & (abs(logFC) >= log2fc_threshold),
      log.adj.P.Val = -log10(!!as.symbol(pval_col)),
      logfc_0 = case_when(
        signif & logFC > 0 ~ "> 0",
        signif & logFC < 0 ~ "< 0")) %>% 
    # Label genes based on there logFC and p-value
    group_by(logfc_0) %>% 
    mutate(
      logfc_rank = dense_rank(dplyr::desc(abs(logFC))),
      pval_rank = dense_rank(dplyr::desc(log.adj.P.Val)),
      label = case_when(
        (logfc_rank < n_genes | pval_rank < n_genes) & signif ~ TRUE,
        gene %in% genes_to_label & (signif | force_label) ~ TRUE,
        TRUE ~ FALSE),
      extra_label = label & (gene %in% extra_genes)) %>% 
    ungroup()
  
  # X-axis breaks
  x_max <- max(abs(plot_data$logFC), na.rm = TRUE)

  if (!is.numeric(x_breaks)) {
    if (x_max > 5) {
      x_breaks <- seq(-1*round(x_max), round(x_max), 2)
    }else{
      x_breaks <- seq(-1*round(x_max), round(x_max))
    }
  }
  
  if (is.null(xlim)) {
    xlim <- c(-1*x_max*1.25, x_max*1.25)
  }
  
  # Y-axis
  y_max <- max(abs(plot_data$log.adj.P.Val), na.rm = TRUE)
  
  if (is.null(ylim)) {
    ylim <- c(0, 1.25*y_max) # this gives padding for group labels
  }
  
  # Coordinates for adding group labels
  if (is.null(group_label_x)) {
    group_label_x <- 2*log2fc_threshold
  }
  
  if (is.null(group_label_y)) {
    group_label_y <- y_max
  }
  
  # Plot
  
  if (is.null(pathway_genes)) {
    vplot <- plot_data %>%
      ggplot(aes(x = logFC, y = log.adj.P.Val, col = logfc_0, alpha = signif)) +
      geom_point(size = point_size)
  
  # For plotting all genes in a given pathway (`pathway_genes`) regardless of significance
  }else{
    plot_data <-  plot_data %>% 
      mutate(
        label = FALSE,
        show = ifelse(gene %in% pathway_genes, TRUE, FALSE))
    
    vplot <- plot_data %>% 
      ggplot(aes(x = logFC, y = log.adj.P.Val, col = show, alpha = show)) +
      # Layer `show == TRUE` points on top
      geom_point(data = plot_data %>% dplyr::filter(!show)) +
      geom_point(data = plot_data %>% dplyr::filter(show), size = point_size)
  }
  
  
  vplot <- vplot +
    geom_vline(xintercept = c(-1*log2fc_threshold, log2fc_threshold), linetype = "dashed", linewidth = 1, color = "#6B6B6B") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", linewidth = 1, color = "#6B6B6B") +
    plot_theme +
    scale_color_manual(values = c("< 0" = left_color, "> 0" = right_color, "FALSE" = "grey")) +
    scale_alpha_manual(values = c("TRUE" = 0.75, "FALSE" = 0.3)) +
    guides(color = "none", alpha = "none") +
    labs(x = "log<sub>2</sub>FC", y = pval_label) +
    scale_x_continuous(limits = xlim, breaks = x_breaks) +
    scale_y_continuous(limits = ylim, breaks = y_breaks) +
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown()) +
    annotate(
      geom = "richtext",
      x = group_label_x, y = group_label_y,
      label = paste0("**",group1, "**"), size = LABEL_TEXTSIZE_DBL, fill = NA, label.color = NA,
      label.size = 0, label.padding = unit(rep(0, 4), "pt")) +
    annotate(
      geom = "richtext",
      x = -1*group_label_x, y = group_label_y, 
      label = paste0("**",group2, "**"), size = LABEL_TEXTSIZE_DBL, fill = NA, label.color = NA,
      label.size = 0, label.padding = unit(rep(0, 4), "pt")) +
    ggtitle(title)
  
  if (box_label) {
    vplot <- vplot +
      new_scale_color() +  
      geom_label_repel(
        data = plot_data %>% dplyr::filter(label), 
        aes(label = gene, color = extra_label), nudge_x = nudge_x, nudge_y = nudge_y, segment.size = label_line_size,
        alpha = 1, fontface = "italic", size = LABEL_TEXTSIZE_DBL*0.75, max.overlaps = max_overlaps) +
      scale_colour_manual(values = c("black", extra_label_color)) +
      guides(color = "none")
    
  }else{
    vplot <- vplot +
      new_scale_color() +  
      geom_text_repel(
        data = plot_data  %>% dplyr::filter(label), 
        aes(label = gene, color = extra_label), nudge_x = nudge_x, nudge_y = nudge_y, segment.size = label_line_size,
        alpha = 1, fontface = "italic", size = LABEL_TEXTSIZE_DBL*0.7, max.overlaps = max_overlaps) +
      scale_colour_manual(values = c("black", extra_label_color)) +
      guides(color = "none")
  }
  
  return(vplot) 
}

plot_gsea <- function(data, title = "GSEA", filter_to_signif = TRUE, color_vec = c("red", "blue"), scale_gs_name = 1, adjpval_threshold = 0.05, bar_width = NULL){
  
  if (filter_to_signif) {
    data <- data %>% 
      dplyr::filter(qvalue < adjpval_threshold)
  }
  
  # If data has column `signif` (which has T/F or "q < 0.1" and "q >= 0.1" for instance)
  # use as alpha
  if (sum(colnames(data) %in% c("signif")) > 0) {
    plot <- data %>% 
      ggplot(aes(x = NES, y = ID, fill = enrich, alpha = signif)) +
      scale_alpha_manual(values = c("signif" = 1, "ns" = 0.25)) +
      labs(x = "Normalized Enrichment Score", y = "", fill = "Enrichment", alpha = "q-value")
    
  } else{
    plot <- data %>% 
      ggplot(aes(x = NES, y = ID, fill = enrich)) +
      labs(x = "Normalized Enrichment Score", y = "", fill = "Enrichment") 
  }
  
  plot <- plot +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = bar_width) +
    plot_theme +
    theme(
      axis.text.y = element_text(size = TEXTSIZE_DBL*scale_gs_name),
      # axis.text.x = element_text(size = TEXTSIZE_DBL*0.75)
      ) +
    scale_fill_manual(values = color_vec) +
    scale_x_continuous(breaks = seq(-3, 3, 1)) +
    guides(alpha = "none", fill = "legend") +
    ggtitle(title)
  
  return(plot)
}

plot_gsea_enrichment <- function(gsea, geneset, title = geneset, show_qval = TRUE,
                           left_enrich_label = "", right_enrich_label = "", left_color = "red", right_color = "blue"){
  
  nes <- gsea@result %>% 
    dplyr::filter(ID == geneset) %>% 
    pull(NES)
  
  pval <- gsea@result %>% 
    dplyr::filter(ID == geneset) %>% 
    pull(pvalue)
  
  qval <- gsea@result %>% 
    dplyr::filter(ID == geneset) %>% 
    pull(qvalue)
  
  anno_y <- gsea@result %>% 
    dplyr::filter(ID == geneset) %>% 
    pull(enrichmentScore)
  
  if(anno_y < 0){
    anno_y <- anno_y*0.6
  }else{
    anno_y <-  anno_y*0.3
  }
  
  if (show_qval) {
    plot_anno <- paste0(
      paste("NES =", format(round(nes, 2), nsmall = 2)), "\n", 
      paste("p =", format_pvals(pval)), "\n",
      paste("q =", format_pvals(qval)))
  }else{
    plot_anno <- paste0(
      paste("NES =", format(round(nes, 2), nsmall = 2)), "\n", 
      paste("p =", format_pvals(pval)), "\n")
  }


  # Plot 1 - running enrichment score
  plot1 <- enrichplot::gseaplot2(gsea, geneSetID = geneset, title = title, subplots = 1) +
    theme(
      plot.title = plot_theme$plot.title,
      axis.text = plot_theme$axis.text.y,
      axis.title = plot_theme$axis.title,
      axis.ticks = element_line(linewidth = 2, color = "black"),
      axis.line = element_line(linewidth = 1.25, color = "black"),
      axis.ticks.length = unit(2, "mm"),
    ) +
    scale_x_discrete(breaks = 0) +
    scale_y_continuous(n.breaks = 4) +
    annotate(
      "label", x = 1000, y = anno_y, size = LABEL_TEXTSIZE_DBL,
      label = plot_anno, hjust = 0, label.size = 0) +
    labs(y = "Running\nEnrichment Score") +
    coord_cartesian(clip = "off") 
  
  # Plot 2 - color bar
  # Code is mainly pulled from source (https://rdrr.io/github/GuangchuangYu/enrichplot/src/R/gseaplot.R)
  gsdata <- enrichplot:::gsInfo(gsea, geneset)
  
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(11) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values = "black") +
    coord_cartesian(clip = "off") 
  
  
  v <- seq(1, sum(gsdata$position), length.out=9)
  inv <- findInterval(rev(cumsum(gsdata$position)), v)
  if (min(inv) == 0) inv <- inv + 1
  
  col <- colorRampPalette(c(right_color, "white", left_color))(10)
  
  ymin <- min(p2$data$ymin)
  yy <- ymin*1.25
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  d <- data.frame(ymin = ymin, ymax = yy,
                  xmin = xmin,
                  xmax = xmax,
                  col = col[unique(inv)])
  plot2 <- p2 + geom_rect(
    aes_(xmin=~xmin,
         xmax=~xmax,
         ymin=~ymin,
         ymax=~ymax,
         fill=~I(col)),
    data=d,
    alpha=.9,
    inherit.aes=FALSE)
  
  # Plot 3 - genes in ranked list
  plot3 <- enrichplot::gseaplot2(gsea, geneSetID = geneset, title = title, subplots = 3) +
    theme(
      plot.title = plot_theme$plot.title,
      axis.text = plot_theme$axis.text.y,
      axis.title = plot_theme$axis.title,
      axis.ticks = element_line(linewidth = 2, color = "black"),
      axis.line = element_line(linewidth = 1.25, color = "black"),
      axis.ticks.length = unit(2, "mm")) +
    annotate("text", x = 1500, y = Inf, size = LABEL_TEXTSIZE_DBL,
             label = left_enrich_label) +
    annotate("text", x = 10000, y = Inf, size = LABEL_TEXTSIZE_DBL,
             label = right_enrich_label) +
    coord_cartesian(clip = "off")
  
  
  # Combine
  layout <- 
    "A
     B
     C"
  gsea_plot <- plot1 + wrap_elements(plot2, clip = FALSE) + plot3 + plot_layout(design = layout)
  
  return(gsea_plot)
}


# Function from https://rdrr.io/github/PeeperLab/Rtoolbox/src/R/ReplotGSEA.R
replot_gsea <- function(path, gene.set, class.name = "", metric.range,
                       enrichment.score.range,
                       left_color = "blue", right_color = "red",
                       left_enrich_label = NULL, right_enrich_label = NULL) {
  
  ## Script by Thomas Kuilman
  ## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
  ## gene.set argument: name of the gene set (e.g. V$AP1_Q2).
  ## It is used in a grep command, so multiple matching is possible.
  ## Also, R regular expressions can be handled, e.g. "IL2[0-9]$"
  ## Leading "V$" from gene set names are stripped to allow using the grep command.
  ## In case of multiple grep matches a warning is given and the first option is plotted.
  ## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)
  ## metric.range: the range of the metric; defaults to [min(DEFINED RANGE), max(DEFINED RANGE)]
  ## enrichment.score.range: the range of the enrichment score; defaults to [min(ENRICHMENT SCORE), max(ENRICHMENT SCORE)]
  
  # Modfied 
  # `path` points to the path of edb folder directly
  # `class.name` default is "" and that label on plot no longer pulls from `gsea.template`
  # Added custom colors and labels
  
  
  if (missing(path)) {
    stop("Path argument is required")
  }
  if (!file.exists(path)) {
    stop("The path folder could not be found. Please change the path")
  }
  if (missing(gene.set)) {
    stop("Gene set argument is required")
  }
  
  ## Load .rnk data
  path.rnk <- list.files(path = path,
                         pattern = ".rnk$", full.names = TRUE)
  gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
  if (missing(metric.range)) {
    metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
  }  
  
  ## Load .edb data
  path.edb <- list.files(path = path,
                         pattern = ".edb$", full.names = TRUE)
  gsea.edb <- read.delim(file = path.edb,
                         header = FALSE, stringsAsFactors = FALSE)
  gsea.edb <- unlist(gsea.edb)
  gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
  gsea.metric <- unlist(strsplit(gsea.metric, " "))
  gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
  gsea.metric <- gsub("METRIC=", "", gsea.metric)
  gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
  
  # Select the right gene set
  if (length(gsea.edb) == 0) {
    stop(paste("The gene set name was not found, please provide",
               "a correct name"))
  }
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
    warning(paste("More than 1 gene set matched the gene.set",
                  "argument; the first match is plotted"))
  }
  gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]
  
  # Get template name
  gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
  gsea.edb <- unlist(strsplit(gsea.edb, " "))
  gsea.template <- gsea.edb[1]
  
  # Get gene set name
  gsea.gene.set <- gsea.edb[2]
  gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
  
  # Get enrichment score
  gsea.enrichment.score <- gsea.edb[3]
  gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
  
  # Get gene set name
  gsea.normalized.enrichment.score <- gsea.edb[4]
  gsea.normalized.enrichment.score <- gsub("NES=", "",
                                           gsea.normalized.enrichment.score)
  
  # Get nominal p-value
  gsea.p.value <- gsea.edb[5]
  gsea.p.value <- gsub("NP=", "", gsea.p.value)
  gsea.p.value <- as.numeric(gsea.p.value)
  
  # Get FDR
  gsea.fdr <- gsea.edb[6]
  gsea.fdr <- gsub("FDR=", "", gsea.fdr)
  gsea.fdr <- as.numeric(gsea.fdr)
  
  # Get hit indices
  gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
  gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
  gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
  gsea.hit.indices <- as.integer(gsea.hit.indices)
  
  # Get ES profile
  gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
  gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
  gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
  gsea.es.profile <- as.numeric(gsea.es.profile)
  
  # Set enrichment score range
  if (missing(enrichment.score.range)) {
    enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
  }

  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
  layout.show(gsea.layout)
  
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
       c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       ylim = enrichment.score.range,
       main = list(gsea.gene.set, font = 1, cex = 1),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score), adj = c(0, 0))
  } else {
    text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))
  abline(v = gsea.hit.indices, lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))
  rank.colors <- gsea.rnk$metric - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c(right_color, "white", left_color))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea.rnk$metric)))
  box()
  text(length(gsea.rnk$metric) / 2, 0.7,
       labels = class.name) #ifelse(!missing(class.name), class.name, gsea.template)
  text(length(gsea.rnk$metric) * 0.01, 0.7, left_enrich_label, adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.7, right_enrich_label, adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk$metric, digits = 2))
  plot(gsea.rnk$metric, type = "n", xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
       ylim = metric.range, yaxs = "i",
       ylab = if(gsea.metric == "None") {"Ranking metric"} else {gsea.metric},
       panel.first = abline(h = seq(metric.range[1] / 2,
                                    metric.range[2] - metric.range[1] / 4,
                                    metric.range[2] / 2), col = "gray95", lty = 2))
  
  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
          xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
          ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
          ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), space = 0, add = TRUE)
  box()
  
  # Reset to default
  par(def.par)

  

  

  
}
