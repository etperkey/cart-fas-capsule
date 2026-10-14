# 30_render_figures.R
# ------------------------------------------------------------------
# Assembles per-panel .rds objects from stages 12, 13, 20 into the
# publication-ready composite Figure 1 and Supplementary Figure 1.
#
# Inputs:  data/derived/panels/*.rds
# Outputs:
#   results/figure1_composite.pdf
#   results/supfig1_composite.pdf
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

load_panel <- function(id) {
  p <- derived_path("panels", paste0(id, ".rds"))
  if (!file.exists(p)) {
    warning(sprintf("[30_render_figures] Missing panel: %s", id))
    return(NULL)
  }
  readRDS(p)
}

fig_panels <- list(
  "1A" = load_panel("fig1a_drs_efs_km"),
  "1B" = load_panel("fig1b_drs_gene_maxstat"),
  "1C" = load_panel("fig1c_drs_nofas_efs_km"),
  "1D" = load_panel("fig1d_cd19xfas_efs_km"),
  "1E" = load_panel("fig1e_cd19xfas_os_km"),
  "1F" = load_panel("fig1f_nci_cd19xfas_os_km"),
  "1G" = load_panel("fig1g_cibersortx_cd19_scatter"),
  "1H" = load_panel("fig1h_cibersortx_cd19xfas_efs_km"),
  "1I" = load_panel("fig1i_cibersortx_cd19xfas_os_km")
)
fig_panels <- Filter(Negate(is.null), fig_panels)

sup_panels <- list(
  "S1B" = load_panel("sup1b_cd19low_fas_efs_km"),
  "S1C" = load_panel("sup1c_cd19low_fas_os_km"),
  "S1D" = load_panel("sup1d_cibersortx_cd19_fas_scatter"),
  "S1E" = load_panel("sup1e_cibersortx_variants"),
  "S1F" = load_panel("sup1f_nci_pfs_km")
)
sup_panels <- Filter(Negate(is.null), sup_panels)

### Figure 1 composite ----------------------------------------------
if (length(fig_panels) > 0) {
  message(sprintf("[30_render_figures] Figure 1 composite from %d panels.",
                  length(fig_panels)))
  tagged <- lapply(seq_along(fig_panels),
                   function(i) fig_panels[[i]] +
                     patchwork::plot_annotation(title = names(fig_panels)[i]))
  fig1 <- patchwork::wrap_plots(fig_panels, ncol = 3) +
    patchwork::plot_annotation(
      title = "Figure 1 — FAS/FASLG and CD19 CAR-T outcomes",
      tag_levels = "A"
    )
  ggsave(results_path("figure1_composite.pdf"),
         fig1, width = 24, height = 24, device = cairo_pdf, limitsize = FALSE)
}

### Supplementary Figure 1 composite --------------------------------
if (length(sup_panels) > 0) {
  message(sprintf("[30_render_figures] Sup Fig 1 composite from %d panels.",
                  length(sup_panels)))
  sup1 <- patchwork::wrap_plots(sup_panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "Supplementary Figure 1 — Sensitivity analyses",
      tag_levels = list(c("B", "C", "D", "E", "F"))
    )
  ggsave(results_path("supfig1_composite.pdf"),
         sup1, width = 18, height = 24, device = cairo_pdf, limitsize = FALSE)
}

# Sup 1A (correlation heatmap) is already a standalone PDF produced by
# scripts/12_bulk_survival.R because its native container is a grob, not
# a ggplot. Keep it separate rather than wrapping it with patchwork.

message("[30_render_figures] Done.")
