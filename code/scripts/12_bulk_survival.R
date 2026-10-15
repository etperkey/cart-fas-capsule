# 12_bulk_survival.R
# ------------------------------------------------------------------
# Reproduces Fig 1A, 1B, 1C, 1D, 1E and Supplementary Fig 1A, 1B, 1C.
#
# Inputs:
#   data/derived/panels/cart_df.rds
#   data/derived/panels/cutpoints.rds
#
# Outputs:
#   results/fig1a_drs_efs_km.pdf
#   results/fig1b_drs_gene_maxstat.pdf
#   results/fig1c_drs_nofas_efs_km.pdf
#   results/fig1d_cd19xfas_efs_km.pdf
#   results/fig1e_cd19xfas_os_km.pdf
#   results/sup1a_drs_gene_corr_heatmap.pdf
#   results/sup1b_cd19low_fas_efs_km.pdf
#   results/sup1c_cd19low_fas_os_km.pdf
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

cart_df   <- readRDS(panels_path("cart_df.rds"))
cutpoints <- readRDS(panels_path("cutpoints.rds"))

fas_cut      <- cutpoints$fas_efs$cutpoint[["cutpoint"]]
cd19_cut     <- cutpoints$cd19_efs$cutpoint[["cutpoint"]]
drs_cut      <- cutpoints$drs_efs$cutpoint[["cutpoint"]]
drsNoFAS_cut <- cutpoints$drsNoFAS_efs$cutpoint[["cutpoint"]]

dlbcl_hgbl_pre <- cart_df %>%
  dplyr::filter(biopsy_time == "Pretreatment",
                diagnosis %in% c("DLBCL", "HGBL"))

### Fig 1A — DRS maxstat EFS -----------------------------------------
message("[12_bulk_survival] Fig 1A — DRS maxstat KM EFS.")
fig1a <- dlbcl_hgbl_pre %>%
  dplyr::filter(!is.na(efs_months), !is.na(efs_status)) %>%
  mutate(drs_group = factor(ifelse(
    percent_rank(DRS_gsva) > drs_cut, "DRS high", "DRS low"))) %>%
  plot_survival(
    time_col = "efs_months", status_col = "efs_status",
    ylab = "Probability of EFS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "drs_group", color_vec = drs_colors,
    curve_size = 1.8, censor_size = 6,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Event-free survival", p_adjust = "none",
    pvals = TRUE, include_legend = TRUE)
save_panel("fig1a_drs_efs_km",
           fig1a$plot + theme(legend.text = element_markdown()),
           width = 10, height = 7)

### Fig 1B — per-DRS-gene maxstat chi-square + BH --------------------
message("[12_bulk_survival] Fig 1B — per-DRS-gene maxstat FDR.")
drs_genes <- drs_gene_sets$DRS_genes
max_stat_df <- purrr::map_dfr(drs_genes, function(g) {
  gcol <- paste0(g, "_gene")
  if (!gcol %in% names(cart_df)) return(NULL)
  cut_res <- maxstat_cut(dlbcl_hgbl_pre, "efs_months", "efs_status", gcol)
  cut_val <- cut_res$cutpoint[["cutpoint"]]
  data <- dlbcl_hgbl_pre %>%
    dplyr::filter(!is.na(efs_months), !is.na(efs_status)) %>%
    mutate(gene_group = factor(ifelse(
      percent_rank(.data[[gcol]]) > cut_val, "high", "low")))
  diff <- survdiff(Surv(efs_months, efs_status) ~ gene_group, data, rho = 0)
  tibble(gene = g,
         cutpoint = cut_val,
         pval  = 1 - pchisq(diff$chisq, 1),
         chisq = diff$chisq,
         stat  = sqrt(diff$chisq))
})

fig1b <- max_stat_df %>%
  mutate(
    p_adj      = p.adjust(pval, method = "BH"),
    adj_signif = ifelse(p_adj < 0.05, "p adj. < 0.05", "p adj. ≥ 0.05"),
    gene       = paste0("*", gene, "*"),
    gene       = fct_reorder(gene, stat)
  ) %>%
  ggplot(aes(x = gene, y = stat)) +
  geom_point(size = 8, aes(color = adj_signif)) +
  plot_theme +
  scale_color_manual(values = c("red", "black")) +
  labs(y = "Maximal standardized\nlog-rank statistic",
       x = "DRS genes", color = NULL) +
  theme(axis.text.x = element_markdown(angle = 90)) +
  geom_hline(yintercept = abs(sqrt(qchisq(0.95, 1))),
             linetype = 2, linewidth = 1)
save_panel("fig1b_drs_gene_maxstat", fig1b, width = 15, height = 8)

### Sup 1A — DRS gene correlation heatmap ----------------------------
message("[12_bulk_survival] Sup 1A — DRS gene correlation heatmap.")
ordered_genes <- max_stat_df %>% arrange(desc(stat)) %>% pull(gene)
present <- ordered_genes[paste0(ordered_genes, "_gene") %in% names(cart_df)]
mat <- cart_df %>%
  dplyr::select(sample, all_of(paste0(present, "_gene"))) %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()
hmap_mat <- cor(mat, method = "pearson")
gene_labels <- str_remove(rownames(hmap_mat), "_gene")
color_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                        c(pal_gsea()(12)[c(2, 4)], "white",
                          pal_gsea()(12)[c(9, 11)]))
hmap <- Heatmap(
  hmap_mat, col = color_fun, name = "Pearson correlation",
  cluster_columns = FALSE, cluster_rows = FALSE,
  row_labels = gt_render(gene_labels),
  column_labels = gt_render(gene_labels),
  row_names_gp = gpar(fontface = "italic"),
  column_names_gp = gpar(fontface = "italic"),
  heatmap_legend_param = list(direction = "horizontal")
)
# Open a null PDF device so ComplexHeatmap's internal `draw()` has an
# active device to measure against. Without this, grid::grid.grabExpr
# triggers an implicit `pdf("Rplots.pdf")` in the CWD, which is read-only
# in-container.
grDevices::pdf(NULL)
hmap_gb <- grid::grid.grabExpr(draw(
  hmap,
  column_title = "Pairwise correlation of genes in the Death Receptor Score",
  heatmap_legend_side = "bottom"))
grDevices::dev.off()
saveRDS(hmap_gb, panels_path("sup1a_drs_gene_corr_heatmap.rds"))
grDevices::cairo_pdf(results_path("sup1a_drs_gene_corr_heatmap.pdf"),
                     width = 7, height = 7, bg = "transparent")
grid::grid.draw(hmap_gb); dev.off()

### Fig 1C — DRS without FAS, maxstat EFS ----------------------------
message("[12_bulk_survival] Fig 1C — DRS (no FAS) maxstat KM EFS.")
fig1c <- dlbcl_hgbl_pre %>%
  dplyr::filter(!is.na(efs_months), !is.na(efs_status)) %>%
  mutate(drs_group = factor(ifelse(
    percent_rank(DRS_noFAS_gsva) > drsNoFAS_cut, "DRS high", "DRS low"))) %>%
  plot_survival(
    time_col = "efs_months", status_col = "efs_status",
    ylab = "Probability of EFS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "drs_group", color_vec = drs_colors,
    curve_size = 1.8, censor_size = 6,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Event-free survival", p_adjust = "none",
    pvals = TRUE, include_legend = TRUE)
save_panel("fig1c_drs_nofas_efs_km",
           fig1c$plot + theme(legend.text = element_markdown()),
           width = 10, height = 7)

### Fig 1D / 1E — CD19 × FAS maxstat, EFS and OS ---------------------
message("[12_bulk_survival] Fig 1D/1E — CD19 × FAS maxstat KM.")
fas_cd19_df <- dlbcl_hgbl_pre %>%
  dplyr::filter(!is.na(efs_months), !is.na(efs_status)) %>%
  mutate(
    CD19_group = factor(ifelse(percent_rank(CD19_gene) > cd19_cut, "high", "low")),
    FAS_group  = factor(ifelse(percent_rank(FAS_gene)  > fas_cut,  "high", "low")),
    group      = paste("*FAS*", FAS_group, "and", "*CD19*", CD19_group),
    group      = factor(group, levels = names(fas_cd19_colors))
  )

fig1d <- fas_cd19_df %>%
  plot_survival(
    time_col = "efs_months", status_col = "efs_status",
    ylab = "Probability of EFS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "group", color_vec = fas_cd19_colors,
    curve_size = 1.8, censor_size = 7,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Event-free survival", p_adjust = "BH",
    pvals = TRUE, pval_y = 0.2, include_legend = TRUE)
save_panel("fig1d_cd19xfas_efs_km",
           fig1d$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

fig1e <- fas_cd19_df %>%
  dplyr::filter(!is.na(os_months), !is.na(os_status)) %>%
  plot_survival(
    time_col = "os_months", status_col = "os_status",
    ylab = "Probability of OS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "group", color_vec = fas_cd19_colors,
    curve_size = 1.8, censor_size = 6,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Overall survival", p_adjust = "BH", pval_y = 0.2,
    pvals = TRUE, include_legend = TRUE)
save_panel("fig1e_cd19xfas_os_km",
           fig1e$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

### Sup 1B / 1C — CD19low × FAS 1-vs-3 EFS and OS --------------------
message("[12_bulk_survival] Sup 1B/1C — CD19low × FAS 1-vs-3 KM.")
fas_cd19_1vs3_df <- fas_cd19_df %>%
  mutate(
    group = as.character(group),
    group = ifelse(group == "*FAS* low and *CD19* low",
                   "*FAS* low and *CD19* low", "All other cases"),
    group = factor(group, levels = names(fas_cd19_1vs3_colors))
  )

sup1b <- fas_cd19_1vs3_df %>%
  plot_survival(
    time_col = "efs_months", status_col = "efs_status",
    ylab = "Probability of EFS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "group", color_vec = fas_cd19_1vs3_colors,
    curve_size = 1.8, censor_size = 7,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Event-free survival", p_adjust = "BH", pval_y = 0.2,
    pvals = TRUE, include_legend = TRUE)
save_panel("sup1b_cd19low_fas_efs_km",
           sup1b$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

sup1c <- fas_cd19_1vs3_df %>%
  dplyr::filter(!is.na(os_months), !is.na(os_status)) %>%
  plot_survival(
    time_col = "os_months", status_col = "os_status",
    ylab = "Probability of OS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "group", color_vec = fas_cd19_1vs3_colors,
    curve_size = 1.8, censor_size = 6,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Overall survival", p_adjust = "BH", pval_y = 0.2,
    pvals = TRUE, include_legend = TRUE)
save_panel("sup1c_cd19low_fas_os_km",
           sup1c$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

message("[12_bulk_survival] Done.")
