# 13_nci_bulk_survival.R
# ------------------------------------------------------------------
# Reproduces Fig 1F and Supplementary Fig 1F — NCI R-CHOP cohort KM
# using the CD19 and FAS maxstat cutpoints learned on Stanford+ZUMA-1.
#
# Inputs:
#   data/derived/CAR-T_FAS/output/nci_vsd.rds       (VST matrix, NCI)
#   data/derived/CAR-T_FAS/output/nci_metadata.csv  (OS/PFS clindata)
#   data/derived/panels/cutpoints.rds               (fas/cd19 cutpoints)
#
# Outputs:
#   results/fig1f_nci_cd19xfas_os_km.pdf
#   results/sup1f_nci_pfs_km.pdf
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

nci_vsd_path  <- get_file("output", "nci_vsd.rds")
nci_meta_path <- get_file("output", "nci_metadata.csv")

if (!file.exists(nci_vsd_path) || !file.exists(nci_meta_path)) {
  stop(sprintf(paste0(
    "[13_nci_bulk_survival] Missing NCI inputs.\n",
    "Expected:\n  %s\n  %s\n",
    "Run scripts/00_prepare_data.R first (requires dbGaP phs001444 counts), ",
    "or drop in the bundled derived data asset."),
    nci_vsd_path, nci_meta_path))
}

cutpoints <- readRDS(panels_path("cutpoints.rds"))
fas_cut  <- cutpoints$fas_efs$cutpoint[["cutpoint"]]
cd19_cut <- cutpoints$cd19_efs$cutpoint[["cutpoint"]]

message("[13_nci_bulk_survival] Reading NCI VST + clindata.")
nci_vsd      <- readRDS(nci_vsd_path)
nci_meta     <- read_csv(nci_meta_path, show_col_types = FALSE)
nci_cart_df  <- build_cohort_df(nci_vsd, nci_meta)

### Sup 1F — NCI PFS at median FAS ------------------------------------
message("[13_nci_bulk_survival] Sup 1F — NCI PFS, median FAS cut.")
sup1f <- nci_cart_df %>%
  dplyr::filter(!is.na(pfs_years), !is.na(pfs_status)) %>%
  mutate(fas_group = factor(ifelse(
    FAS_gene > median(FAS_gene), "*FAS* high", "*FAS* low"))) %>%
  plot_survival(
    time_col = "pfs_years", status_col = "pfs_status",
    ylab = "Probability of PFS",
    end_time = 5, xlab = "Time in years", break_time_by = 1,
    group = "fas_group", color_vec = fas_colors,
    curve_size = 1.8, censor_size = 7,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Progression-free survival", p_adjust = "none",
    pvals = TRUE, include_legend = TRUE)
save_panel("sup1f_nci_pfs_km",
           sup1f$plot + theme(legend.text = element_markdown()),
           width = 10, height = 7)

### Fig 1F — NCI OS with Stanford+ZUMA maxstat cutpoints -------------
message("[13_nci_bulk_survival] Fig 1F — NCI OS with cross-cohort cuts.")
tmp_df <- nci_cart_df %>%
  mutate(
    CD19_group = factor(ifelse(percent_rank(CD19_gene) > cd19_cut, "high", "low")),
    FAS_group  = factor(ifelse(percent_rank(FAS_gene)  > fas_cut,  "high", "low")),
    group      = paste("*FAS*", FAS_group, "and", "*CD19*", CD19_group),
    group      = factor(group, levels = names(fas_cd19_colors))
  )

fig1f <- tmp_df %>%
  dplyr::filter(!is.na(os_years), !is.na(os_status)) %>%
  plot_survival(
    time_col = "os_years", status_col = "os_status",
    ylab = "Probability of OS",
    end_time = 5, xlab = "Time in years", break_time_by = 1,
    group = "group", color_vec = fas_cd19_colors,
    curve_size = 1.8, censor_size = 7,
    add_num_to_name = FALSE, risk_table = TRUE,
    title = "Overall survival", p_adjust = "BH", pval_y = 0.2,
    pvals = TRUE, include_legend = TRUE)
save_panel("fig1f_nci_cd19xfas_os_km",
           fig1f$plot + theme(legend.text = element_markdown()),
           width = 12, height = 7)

message("[13_nci_bulk_survival] Done.")
