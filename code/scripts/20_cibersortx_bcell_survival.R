# 20_cibersortx_bcell_survival.R
# ------------------------------------------------------------------
# Reproduces Fig 1G, 1H, 1I and Supplementary Fig 1D, 1E.
#
# CIBERSORTx High-Resolution runs on the CIBERSORTx web portal (not as an
# R package), so the capsule does not perform the deconvolution itself.
# This stage consumes the saved CIBERSORTx outputs and produces the
# downstream figures.
#
# Inputs (REQUIRED — the CIBERSORTx run output must be deposited here):
#   data/derived/cibersortx/stanford_cibersortx_cell_fractions.csv
#   data/derived/cibersortx/stanford_cibersortx_bcell_gep.csv  (gene-expression profile,
#                                                              rows = genes, cols = samples)
#   data/derived/panels/cart_df.rds        (from stage 10 — bulk CD19 for scatter)
#   data/derived/panels/cutpoints.rds      (from stage 11)
#
# Outputs:
#   results/fig1g_cibersortx_cd19_scatter.pdf    # bulk CD19 vs CIBERSORTx B-cell fraction
#   results/fig1h_cibersortx_cd19xfas_efs_km.pdf # B-cell GEP CD19 x FAS EFS
#   results/fig1i_cibersortx_cd19xfas_os_km.pdf  # B-cell GEP CD19 x FAS OS
#   results/sup1d_cibersortx_cd19_fas_scatter.pdf
#   results/sup1e_cibersortx_variants.pdf
#
# CIBERSORTx settings documented in the paper's Methods:
#   signature matrix: LM22 merged into 10 B-cell–resolved classes
#   permutations: 500
#   batch correction: B-mode
#
# Note on provenance: the specific CSV layouts above are what Cooper's
# downstream analysis expects. If your CIBERSORTx export uses a different
# layout (e.g., long format, extra header rows), convert before placing
# into `data/derived/cibersortx/`.
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

cibersort_dir <- derived_path("cibersortx")
frac_path <- file.path(cibersort_dir, "stanford_cibersortx_cell_fractions.csv")
gep_path  <- file.path(cibersort_dir, "stanford_cibersortx_bcell_gep.csv")

if (!file.exists(frac_path) || !file.exists(gep_path)) {
  stop(sprintf(paste0(
    "[20_cibersortx_bcell_survival] Missing CIBERSORTx inputs.\n",
    "Expected:\n  %s\n  %s\n\n",
    "These must be exported from the CIBERSORTx web portal using:\n",
    "  - LM22 signature merged into 10 B-cell–resolved classes\n",
    "  - 500 permutations\n",
    "  - B-mode batch correction\n",
    "and placed under data/derived/cibersortx/.\n",
    "If the capsule's bundled data asset supplies them, confirm it was ",
    "mounted at /data."),
    frac_path, gep_path))
}

cart_df   <- readRDS(panels_path("cart_df.rds"))
cutpoints <- readRDS(panels_path("cutpoints.rds"))

frac_df <- read_csv(frac_path, show_col_types = FALSE)
gep_mat <- read_csv(gep_path,  show_col_types = FALSE) %>%
  { stopifnot("gene" %in% names(.)); . } %>%
  column_to_rownames("gene") %>%
  as.matrix()

### Fig 1G — bulk CD19 vs CIBERSORTx B-cell fraction -----------------
message("[20_cibersortx_bcell_survival] Fig 1G — bulk CD19 vs B-cell fraction.")
bcell_cols <- grep("B.cell|Bcell|B_cell|Memory.B|Naive.B|Plasma",
                   names(frac_df), value = TRUE, ignore.case = TRUE)
if (length(bcell_cols) == 0) {
  stop("[20_cibersortx_bcell_survival] Could not locate B-cell columns in ",
       "cibersortx cell-fraction CSV. Expected columns matching ",
       "'B.cell|Bcell|B_cell|Memory.B|Naive.B|Plasma'.")
}
frac_bcell <- frac_df %>%
  mutate(bcell_frac = rowSums(across(all_of(bcell_cols)))) %>%
  dplyr::select(sample = any_of(c("Mixture", "sample", "Sample")), bcell_frac)

fig1g <- cart_df %>%
  dplyr::select(sample, CD19_gene) %>%
  inner_join(frac_bcell, by = "sample") %>%
  ggplot(aes(x = bcell_frac, y = CD19_gene)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  ggpubr::stat_cor() +
  labs(x = "CIBERSORTx B-cell fraction",
       y = "Bulk CD19 (VST)") +
  plot_theme
save_panel("fig1g_cibersortx_cd19_scatter", fig1g, width = 7, height = 7)

### B-cell GEP CD19 & FAS → per-sample B-cell-resolved values --------
message("[20_cibersortx_bcell_survival] Building B-cell GEP CD19/FAS table.")
bcell_df <- tibble(sample = colnames(gep_mat)) %>%
  mutate(
    CD19_bcell = if ("CD19" %in% rownames(gep_mat)) gep_mat["CD19", sample] else NA_real_,
    FAS_bcell  = if ("FAS"  %in% rownames(gep_mat)) gep_mat["FAS",  sample] else NA_real_
  )

merged_df <- cart_df %>%
  dplyr::select(sample, biopsy_time, diagnosis,
                efs_months, efs_status, os_months, os_status) %>%
  inner_join(bcell_df, by = "sample") %>%
  dplyr::filter(biopsy_time == "Pretreatment",
                diagnosis %in% c("DLBCL", "HGBL"),
                !is.na(CD19_bcell), !is.na(FAS_bcell))

cd19b_cut <- maxstat_cut(merged_df, "efs_months", "efs_status", "CD19_bcell")$cutpoint[["cutpoint"]]
fasb_cut  <- maxstat_cut(merged_df, "efs_months", "efs_status", "FAS_bcell")$cutpoint[["cutpoint"]]

merged_df <- merged_df %>%
  mutate(
    CD19_group = factor(ifelse(percent_rank(CD19_bcell) > cd19b_cut, "high", "low")),
    FAS_group  = factor(ifelse(percent_rank(FAS_bcell)  > fasb_cut,  "high", "low")),
    group      = paste("*FAS*", FAS_group, "and", "*CD19*", CD19_group),
    group      = factor(group, levels = names(fas_cd19_colors))
  )

### Sup 1D — B-cell CD19 vs FAS scatter -------------------------------
message("[20_cibersortx_bcell_survival] Sup 1D — B-cell CD19 vs FAS scatter.")
sup1d <- merged_df %>%
  ggplot(aes(x = CD19_bcell, y = FAS_bcell)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  ggpubr::stat_cor() +
  labs(x = "B-cell CD19 (CIBERSORTx HiRes)",
       y = "B-cell FAS (CIBERSORTx HiRes)") +
  plot_theme
save_panel("sup1d_cibersortx_cd19_fas_scatter", sup1d, width = 7, height = 7)

### Fig 1H / 1I — B-cell CD19 × FAS KM EFS and OS --------------------
message("[20_cibersortx_bcell_survival] Fig 1H/1I — B-cell CD19 × FAS KM.")
fig1h <- merged_df %>%
  dplyr::filter(!is.na(efs_months), !is.na(efs_status)) %>%
  plot_survival(
    time_col = "efs_months", status_col = "efs_status",
    ylab = "Probability of EFS",
    end_time = 36, xlab = "Time in months", break_time_by = 12,
    group = "group", color_vec = fas_cd19_colors,
    curve_size = 1.8, censor_size = 7,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Event-free survival", p_adjust = "BH", pval_y = 0.2,
    pvals = TRUE, include_legend = TRUE)
save_panel("fig1h_cibersortx_cd19xfas_efs_km",
           fig1h$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

fig1i <- merged_df %>%
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
save_panel("fig1i_cibersortx_cd19xfas_os_km",
           fig1i$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

### Sup 1E — CIBERSORTx variants assembly ----------------------------
message("[20_cibersortx_bcell_survival] Sup 1E — variants panel.")
sup1e <- (fig1h$plot | fig1i$plot) +
  patchwork::plot_annotation(title = "CIBERSORTx HiRes B-cell survival variants")
save_panel("sup1e_cibersortx_variants", sup1e, width = 20, height = 7)

message("[20_cibersortx_bcell_survival] Done.")
