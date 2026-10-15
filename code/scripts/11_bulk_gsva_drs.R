# 11_bulk_gsva_drs.R
# ------------------------------------------------------------------
# Compute the two sets of cutpoints used throughout the paper:
#   - fas_efs_cut_res, cd19_efs_cut_res  (DLBCL+HGBL pretreatment, EFS)
#   - drs_efs_cut_res (full DRS)
#   - drsNoFAS_efs_cut_res (DRS without FAS)
# These are shared across Fig 1A/C/D/E and the NCI validation Fig 1F.
#
# Panels produced: none directly (prep for stages 12, 13).
#
# Inputs:  data/derived/panels/cart_df.rds
# Outputs: data/derived/panels/cutpoints.rds
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

cart_df <- readRDS(panels_path("cart_df.rds"))

base_pretreatment <- cart_df %>%
  dplyr::filter(biopsy_time == "Pretreatment",
                diagnosis %in% c("DLBCL", "HGBL"))

message("[11_bulk_gsva_drs] Fitting maxstat cutpoints.")
fas_efs_cut_res  <- maxstat_cut(base_pretreatment, "efs_months", "efs_status", "FAS_gene")
cd19_efs_cut_res <- maxstat_cut(base_pretreatment, "efs_months", "efs_status", "CD19_gene")
drs_efs_cut_res  <- maxstat_cut(base_pretreatment, "efs_months", "efs_status", "DRS_gsva")
drsNoFAS_efs_cut_res <- maxstat_cut(base_pretreatment, "efs_months", "efs_status", "DRS_noFAS_gsva")

cutpoints <- list(
  fas_efs      = fas_efs_cut_res,
  cd19_efs     = cd19_efs_cut_res,
  drs_efs      = drs_efs_cut_res,
  drsNoFAS_efs = drsNoFAS_efs_cut_res
)
saveRDS(cutpoints, panels_path("cutpoints.rds"))

for (nm in names(cutpoints)) {
  message(sprintf("[11_bulk_gsva_drs]   %-13s cutpoint = %.4f",
                  nm, cutpoints[[nm]]$cutpoint[["cutpoint"]]))
}
