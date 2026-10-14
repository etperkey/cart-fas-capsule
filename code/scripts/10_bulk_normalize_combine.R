# 10_bulk_normalize_combine.R
# ------------------------------------------------------------------
# Loads the combined Stanford + ZUMA-1 VST matrix and metadata, computes
# GSVA scores and per-sample single-gene expression, and writes the
# per-sample analysis frame used by stages 11, 12.
#
# Panels produced: none directly (prep for Fig 1A–E, Sup 1A–C).
#
# Inputs:
#   data/derived/DLBCL_CART_tumorBulkRNAseq/output/zuma1_stanford_vsd.rds
#   data/derived/DLBCL_CART_tumorBulkRNAseq/output/zuma1_stanford_metadata.csv
#
# Outputs:
#   data/derived/panels/cart_df.rds        # per-sample tibble
#   data/derived/panels/exprs_mat.rds      # vst-expression matrix (genes × samples)
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

vsd_path <- get_cart_rnaseq_file("output", "zuma1_stanford_vsd.rds")
meta_path <- get_cart_rnaseq_file("output", "zuma1_stanford_metadata.csv")

if (!file.exists(vsd_path) || !file.exists(meta_path)) {
  stop(sprintf(paste0(
    "[10_bulk_normalize_combine] Missing combined Stanford+ZUMA-1 inputs.\n",
    "Expected:\n  %s\n  %s\n",
    "Place the bundled derived data asset under /data/derived/ or run ",
    "scripts/00_prepare_data.R with raw dbGaP counts first."),
    vsd_path, meta_path))
}

message("[10_bulk_normalize_combine] Reading combined VST + clindata.")
zuma1_stanford_vsd <- readRDS(vsd_path)
clindata_tib <- read_csv(meta_path, show_col_types = FALSE)

exprs_mat <- SummarizedExperiment::assay(zuma1_stanford_vsd)
clindata_tib <- clindata_tib %>% dplyr::filter(sample %in% colnames(exprs_mat))

message(sprintf("[10_bulk_normalize_combine] %d samples, %d genes.",
                ncol(exprs_mat), nrow(exprs_mat)))

cart_df <- build_cohort_df(zuma1_stanford_vsd, clindata_tib)

saveRDS(cart_df,   derived_path("panels", "cart_df.rds"))
saveRDS(exprs_mat, derived_path("panels", "exprs_mat.rds"))

message("[10_bulk_normalize_combine] Wrote cart_df + exprs_mat to derived/panels/.")
