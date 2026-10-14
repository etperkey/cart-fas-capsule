# Capsule shim for Alan Cooper's `get_file()` / `get_cart_rnaseq_file()` helpers.
#
# In Cooper's original Kline-lab project (CAR-T_FAS/utils.R), these functions
# resolve paths against the Box drive via the `boxrdrive` package. In this
# Code Ocean capsule we instead resolve against /data (raw inputs) and
# /data/derived (precomputed intermediates), with a laptop fallback via
# paths.R. Signatures are preserved so the ported analysis scripts in
# `scripts/` and the vendored notebook can keep calling these verbatim.

if (!exists("CAPSULE_ROOT")) {
  source(file.path(dirname(sys.frame(1)$ofile %||% "."), "paths.R"))
}

# Original CAR-T_FAS project resources live under /data (and /code where
# code was co-located with the project).
get_file <- function(...) {
  rel <- file.path(...)
  candidates <- c(
    raw_path("CAR-T_FAS", rel),
    derived_path("CAR-T_FAS", rel),
    data_path(rel),
    code_path(rel)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    # Fall back to the first candidate so write paths still work.
    hit <- candidates[1]
  }
  hit
}

# The companion DLBCL_CART_tumorBulkRNAseq project supplies the Stanford
# vsd / metadata that bulk_rnaseq_analysis.Rmd reads via this helper.
get_cart_rnaseq_file <- function(...) {
  rel <- file.path(...)
  candidates <- c(
    raw_path("DLBCL_CART_tumorBulkRNAseq", rel),
    derived_path("DLBCL_CART_tumorBulkRNAseq", rel),
    data_path(rel)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) hit <- candidates[1]
  hit
}
