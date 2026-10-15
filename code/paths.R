# Capsule-aware path resolution.
# Code Ocean sets working-tree mounts at /code, /data, /results. On a laptop
# we fall back to the project directory relative to this file.

detect_capsule_root <- function() {
  if (dir.exists("/code") && dir.exists("/data") && dir.exists("/results")) {
    return("/")
  }
  here <- normalizePath(dirname(sys.frame(1)$ofile %||% getwd()),
                        winslash = "/", mustWork = FALSE)
  repeat {
    if (all(c("code", "data", "results") %in% list.files(here))) return(here)
    parent <- dirname(here)
    if (parent == here) break
    here <- parent
  }
  getwd()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

CAPSULE_ROOT <- detect_capsule_root()

code_path    <- function(...) file.path(CAPSULE_ROOT, "code",    ...)
data_path    <- function(...) file.path(CAPSULE_ROOT, "data",    ...)
raw_path     <- function(...) file.path(CAPSULE_ROOT, "data", "raw",     ...)
# derived_path remains /data/derived for READ access to the bundled derived
# data asset (zuma1_stanford_vsd.rds, nci_vsd.rds, cibersortx CSVs). On Code
# Ocean /data is mounted read-only, so panel intermediates (.rds) and figure
# outputs (.pdf) both live under results_path() instead — see panels_path().
derived_path <- function(...) file.path(CAPSULE_ROOT, "data", "derived", ...)
results_path <- function(...) file.path(CAPSULE_ROOT, "results", ...)
panels_path  <- function(...) file.path(CAPSULE_ROOT, "results", "panels", ...)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
