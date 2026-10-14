# Shared setup for all stage scripts.
# Each stage script sources this first; libraries + helpers + palettes.

suppressMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(maxstat)
  library(ggsci)
  library(ggtext)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

SEED <- 123
set.seed(SEED)

# Resolve /code. The `run` entrypoint cd's into /code, and the Dockerfile
# sets WORKDIR /code; for native use the reader is expected to run from the
# capsule's code/ directory. Fall back to walking up for the paths.R case
# where it's used from a non-standard shell.
.CODE_DIR <- (function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  # If CWD already is /code (contains paths.R + scripts/), use it.
  if (file.exists(file.path(cwd, "paths.R")) &&
      dir.exists(file.path(cwd, "scripts"))) return(cwd)
  # If CWD is /code/scripts, step up one.
  if (file.exists(file.path(cwd, "_common.R")) &&
      file.exists(file.path(dirname(cwd), "paths.R"))) return(dirname(cwd))
  # Fall back to walking up.
  here <- cwd
  while (dirname(here) != here) {
    if (file.exists(file.path(here, "paths.R"))) return(here)
    here <- dirname(here)
  }
  stop("Could not locate the capsule's code/ directory from CWD=", cwd)
})()

suppressMessages({
  source(file.path(.CODE_DIR, "paths.R"))
  source(file.path(.CODE_DIR, "utils.R"))
  source(file.path(.CODE_DIR, "analysis_func.R"))
  source(file.path(.CODE_DIR, "plotting_func.R"))
})

ensure_dir(results_path())
ensure_dir(derived_path("panels"))

# Palettes used across stages (lifted verbatim from Cooper's Rmd setup).
fas_colors <- c("*FAS* low" = "purple", "*FAS* high" = "green3")
drs_colors <- c("DRS low" = "purple3", "DRS high" = "green4")
cd19_colors <- fas_colors
names(cd19_colors) <- str_replace(names(cd19_colors), "FAS", "CD19")
lymphoma_colors <- pal_nejm()(4)[c(1, 2, 4, 3)]
names(lymphoma_colors) <- c("HGBL", "DLBCL", "tFL", "PMBL")

# 2x2 FAS × CD19 group palette used by Fig 1D–F.
fas_cd19_colors <- c(
  "*FAS* high and *CD19* high" = "green3",
  "*FAS* high and *CD19* low"  = "skyblue3",
  "*FAS* low and *CD19* high"  = "purple",
  "*FAS* low and *CD19* low"   = "red3"
)
# 1-vs-3 collapse used by Sup 1B/C.
fas_cd19_1vs3_colors <- c(
  "*FAS* low and *CD19* low" = "red3",
  "All other cases"          = "skyblue4"
)

# Death Receptor Signature gene sets (Singh et al. 2020). Three variants of
# the DRS appear in Fig 1A, 1C, and a sensitivity analysis.
drs_gene_sets <- list(
  DRS_genes = c("RIPK1", "TNFRSF10A", "TRADD", "TNFRSF1B", "TNFRSF25", "FAF1",
                "FAS", "CASP10", "CASP8", "RIPK3", "FADD", "TNFRSF9", "BID",
                "CASP3", "CASP7", "TNFRSF10B", "TNFRSF1A"),
  DRS_noFAS_genes = c("RIPK1", "TNFRSF10A", "TRADD", "TNFRSF1B", "TNFRSF25",
                      "FAF1", "CASP10", "CASP8", "RIPK3", "FADD", "TNFRSF9",
                      "BID", "CASP3", "CASP7", "TNFRSF10B", "TNFRSF1A"),
  DRS_noFAS_noTNFRSF9_genes = c("RIPK1", "TNFRSF10A", "TRADD", "TNFRSF1B",
                                "TNFRSF25", "FAF1", "CASP10", "CASP8", "RIPK3",
                                "FADD", "BID", "CASP3", "CASP7", "TNFRSF10B",
                                "TNFRSF1A"),
  immunosign21_genes = c("CD3D", "CD3E", "CD3G", "CD8A", "CD69", "ICOS", "IL5",
                         "CCL2", "CCR2", "IRF1", "GZMA", "GZMB", "GZMK", "GZMM",
                         "PRF1", "CXCR3", "CXCL10", "CXCL11", "STAT1", "STAT4",
                         "TBX21")
)

# Single genes carried through every cohort join.
single_genes_of_interest <- c("FAS", "FASLG", "CD19", "BIRC2", "BIRC3", "XIAP",
                              "RIPK1", "TNFRSF10A", "TRADD", "TNFRSF1B",
                              "TNFRSF25", "FAF1", "CASP10", "CASP8", "RIPK3",
                              "FADD", "TNFRSF9", "BID", "CASP3", "CASP7",
                              "TNFRSF10B", "TNFRSF1A")

#' Build the per-sample analysis frame (gene expression + GSVA scores + clindata).
#'
#' @param vsd A `DESeqTransform` object (output of `varianceStabilizingTransformation`).
#' @param clindata Clinical metadata tibble, one row per sample column of `vsd`.
#' @return tibble joining clindata with `<GENE>_gene` columns and `<SET>_gsva` columns.
build_cohort_df <- function(vsd, clindata) {
  exprs_mat <- SummarizedExperiment::assay(vsd)

  gset_l <- lapply(drs_gene_sets, limma::alias2SymbolTable, species = "Hs")
  gset_l <- lapply(gset_l, function(x) x[!is.na(x)])
  gset_tib <- bind_rows(lapply(
    names(gset_l),
    function(x) tibble(gene = gset_l[[x]],
                       gs = paste0(str_remove(x, "_genes"), "_gsva"))
  ))

  gsva_mat <- fun_gsva(lcpm = exprs_mat, genesets = gset_tib)

  gset_scores_df <- gsva_mat %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    dplyr::rename(Immunosign_score = immunosign21_gsva) %>%
    mutate(Immunosign_group = ifelse(
      Immunosign_score > quantile(Immunosign_score, 0.25), "High", "Low"))

  present_genes <- intersect(single_genes_of_interest, rownames(exprs_mat))
  gene_df <- exprs_mat[present_genes, , drop = FALSE] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    rename_with(.cols = -c("sample"), ~paste0(., "_gene"))

  clindata %>%
    dplyr::filter(sample %in% colnames(exprs_mat)) %>%
    left_join(gset_scores_df, by = "sample") %>%
    left_join(gene_df, by = "sample")
}

#' maxstat cutpoint on a percent-rank-transformed continuous covariate.
#'
#' @return a `surv_cutpoint` object.
maxstat_cut <- function(df, time_col, status_col, value_col, minprop = 0.1) {
  df %>%
    dplyr::filter(!is.na(.data[[time_col]]), !is.na(.data[[status_col]])) %>%
    mutate(.rank = percent_rank(.data[[value_col]])) %>%
    surv_cutpoint(time = time_col, event = status_col,
                  ".rank", minprop = minprop, progressbar = FALSE)
}

#' Save a figure's underlying ggplot/grob to `derived/panels/<id>.rds`
#' AND a standalone PDF to `/results/<id>.pdf`.
save_panel <- function(id, obj, width = 10, height = 7) {
  rds_path <- derived_path("panels", paste0(id, ".rds"))
  pdf_path <- results_path(paste0(id, ".pdf"))
  saveRDS(obj, rds_path)
  ggsave(pdf_path, obj, width = width, height = height, device = cairo_pdf)
  invisible(c(rds = rds_path, pdf = pdf_path))
}
