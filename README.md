# FAS/FASLG CAR-T bystander killing — Code Ocean capsule

Reproducibility capsule for the bulk RNA-seq analyses in:

> Lin MJ\*, **Chorazeczewski JK**\*, Phillips NA, Cooper A, …, Kline J†, Brody JD†.
> **Harnessing Fas/FasL-dependent bystander killing to mitigate antigen escape in CAR-T therapy.**
> *Blood Cancer Discovery* (under revision, BCD-26-0329).
> \* co-first authors · † co-senior authors (JDB corresponding).

Prepared in response to the Blood Cancer Discovery editorial request that
computationally intensive revisions ship a self-contained Code Ocean capsule
with the underlying code and a Code Availability Statement after Data
Availability.

The capsule reproduces **Figure 1A–I** and **Supplementary Figure 1A–F**.
Figures 2–7 and Supplementary Figures 2–6 are wet-lab or in vivo experiments
and are not in capsule scope. Supplementary Figure 1G uses an internal
Mount Sinai scRNA cohort (Brody lab, Scanpy pipeline) that is not deposited
and is also not in scope.

---

## Scope and data provenance

| Dataset | Access | Role in paper |
|---|---|---|
| Stanford/Sworder bulk RNA-seq, pre-treatment DLBCL/HGBL, n=51 | **dbGaP phs003145.v1.p1** (controlled) | axi-cel 3L+ cohort, Fig 1A–E, 1G–I |
| ZUMA-1 bulk RNA-seq, n=67 | **Kite Pharma** (private; Upadhyay 2021) | axi-cel pivotal cohort, Fig 1A–E |
| NCI/Schmitz bulk RNA-seq, n=230 DLBCL | **dbGaP phs001444.v2.p1** / GDC STAR | R-CHOP comparator, Fig 1F / Sup 1F |
| LM22 signature matrix | public (CIBERSORTx) | B-cell deconvolution upstream of Fig 1G–I |

**None of the raw bulk count matrices are redistributed in this capsule.**
The capsule ships, as a private Code Ocean data asset, the **normalized
variance-stabilizing-transform (VST) matrices and clinical metadata** that
are the downstream inputs for the paper's figures. Readers who hold current
dbGaP approvals for phs001444 and phs003145 can regenerate the VST
intermediates from raw counts via `scripts/00_prepare_data.R`.

---

## Capsule layout

```
16-cart-fas-capsule/
├── code/
│   ├── run                       # Code Ocean entry point (bash)
│   ├── paths.R                   # /code, /data, /results path shim
│   ├── utils.R                   # replaces Cooper's boxrdrive helpers
│   ├── analysis_func.R           # Cooper's analysis helpers (verbatim)
│   ├── plotting_func.R           # Cooper's plotting helpers (verbatim)
│   └── scripts/
│       ├── 00_prepare_data.R           # (opt) raw counts → VST + metadata
│       ├── 10_bulk_normalize_combine.R # DESeq2 VST + limma batch-combine
│       ├── 11_bulk_gsva_drs.R          # GSVA Death Receptor Signature
│       ├── 12_bulk_survival.R          # Fig 1A–E, Sup 1A–C
│       ├── 13_nci_bulk_survival.R      # Fig 1F, Sup 1F
│       ├── 20_cibersortx_bcell_survival.R  # Fig 1G–I, Sup 1D–E
│       └── 30_render_figures.R         # assemble panels into final figures
├── data/
│   ├── raw/                      # dbGaP STAR / Salmon counts (NOT shipped)
│   └── derived/                  # VST .rds, clinical CSVs, CIBERSORTx GEP
├── environment/
│   └── Dockerfile                # R 4.2.3 / Bioc 3.16, pinned
├── metadata/
│   └── metadata.yml              # Code Ocean capsule metadata
└── results/                      # figure outputs (PDF + PNG)
```

---

## Software versions (match Methods, Lin & Chorazeczewski et al.)

- R 4.2.3 (rocker/r-ver base)
- Bioconductor 3.16
- DESeq2 1.38.3 · limma 3.54.2 · tximport 1.26.1 · edgeR
- GSVA 1.46.0
- survival 3.4.0 · survminer 0.4.9 · maxstat
- tidyverse / ggsci / patchwork / cowplot
- clusterProfiler + fgsea (loaded for compatibility; not used in final figures)

Seurat / BPCells / SeuratDisk are **not** installed; the Haradhvala scRNA
thread that Cooper's Box project included was dropped from the final
manuscript, and Supplementary Figure 1G (internal Mount Sinai scRNA) uses a
separate Scanpy pipeline maintained by the Brody lab.

---

## Reproduction

**On Code Ocean (one-click reproduce).** The capsule is pre-wired to the
private data asset of derived intermediates; the "Reproducible Run" button
re-runs `code/run all` and writes figure PDFs under `/results`.

**On a local workstation** (laptop or lab server with Docker):

```bash
docker build -t cart-fas-capsule environment/
docker run --rm \
    -v "$(pwd)/data:/data:ro" \
    -v "$(pwd)/results:/results" \
    -v "$(pwd)/code:/code:ro" \
    -w /code \
    cart-fas-capsule bash run all
```

### Staged runs

```bash
# Downstream figures only (uses /data/derived VST matrices)
bash code/run bulk
bash code/run cibersortx
bash code/run figures

# Rebuild intermediates from raw dbGaP counts
# (requires dbGaP phs001444 + phs003145 STAR / Salmon outputs placed under data/raw/)
bash code/run prepare
```

Each stage short-circuits if its output already exists; set `FORCE=1`
to force recomputation.

---

## Figure ↔ script map

| Panel | Script | Output file |
|---|---|---|
| Fig 1A — Stanford+ZUMA DRS KM EFS | `12_bulk_survival.R` | `results/fig1a_drs_efs_km.pdf` |
| Fig 1B — per-DRS-gene maxstat FDR | `12_bulk_survival.R` | `results/fig1b_drs_gene_maxstat.pdf` |
| Fig 1C — DRS-noFAS KM | `12_bulk_survival.R` | `results/fig1c_drs_nofas_efs_km.pdf` |
| Fig 1D — CD19×FAS maxstat EFS | `12_bulk_survival.R` | `results/fig1d_cd19xfas_efs_km.pdf` |
| Fig 1E — CD19×FAS maxstat OS | `12_bulk_survival.R` | `results/fig1e_cd19xfas_os_km.pdf` |
| Fig 1F — NCI CD19×FAS OS (R-CHOP) | `13_nci_bulk_survival.R` | `results/fig1f_nci_cd19xfas_os_km.pdf` |
| Fig 1G — CIBERSORTx fractions vs CD19 bulk | `20_cibersortx_bcell_survival.R` | `results/fig1g_cibersortx_cd19_scatter.pdf` |
| Fig 1H — B-cell GEP CD19×FAS EFS | `20_cibersortx_bcell_survival.R` | `results/fig1h_cibersortx_cd19xfas_efs_km.pdf` |
| Fig 1I — B-cell GEP CD19×FAS OS | `20_cibersortx_bcell_survival.R` | `results/fig1i_cibersortx_cd19xfas_os_km.pdf` |
| Sup 1A — DRS gene correlation heatmap | `12_bulk_survival.R` | `results/sup1a_drs_gene_corr_heatmap.pdf` |
| Sup 1B — CD19low × FAS EFS | `12_bulk_survival.R` | `results/sup1b_cd19low_fas_efs_km.pdf` |
| Sup 1C — CD19low × FAS OS | `12_bulk_survival.R` | `results/sup1c_cd19low_fas_os_km.pdf` |
| Sup 1D — CIBERSORTx B-cell CD19 vs FAS scatter | `20_cibersortx_bcell_survival.R` | `results/sup1d_cibersortx_cd19_fas_scatter.pdf` |
| Sup 1E — CIBERSORTx B-cell EFS/OS variants | `20_cibersortx_bcell_survival.R` | `results/sup1e_cibersortx_variants.pdf` |
| Sup 1F — NCI PFS under R-CHOP | `13_nci_bulk_survival.R` | `results/sup1f_nci_pfs_km.pdf` |

Supplementary Figure 1G (internal scRNA) is **not** in this capsule — see
the Brody lab's Scanpy repository (to be linked in the final Data Availability
Statement).

---

## Draft Code Availability Statement (for Methods)

> **Code Availability.** The analytical pipeline that produces Figure 1A–I
> and Supplementary Figure 1A–F — DESeq2 variance-stabilizing normalization,
> limma `removeBatchEffect` batch combination, GSVA Death Receptor Signature
> scoring, maxstat survival cut-point selection, and Kaplan–Meier
> estimation — is distributed as a Code Ocean compute capsule (DOI pending).
> The capsule bundles the normalized VST matrices and de-identified clinical
> metadata as a private data asset; raw count matrices are available only
> through dbGaP under accessions phs001444.v2.p1 (NCI) and phs003145.v1.p1
> (Stanford/Sworder), and through Kite Pharma for ZUMA-1. The CIBERSORTx
> High-Resolution deconvolution was performed on the CIBERSORTx web portal
> (https://cibersortx.stanford.edu/) using the LM22 signature matrix merged
> into 10 B-cell–resolved classes, 500 permutations, and B-mode batch
> correction; the resulting fraction and gene-expression-profile matrices
> are included in the capsule's data asset.

---

## Provenance

The analysis scripts in `code/analysis_func.R` and `code/plotting_func.R`
are Alan Cooper's work from the Kline-lab project directory `CAR-T_FAS/`
(originally stored on UChicago Box), lightly adapted to use the Code Ocean
`/data` and `/results` path convention. The companion project
`DLBCL_CART_tumorBulkRNAseq/` supplied the Stanford VST intermediates.
Capsule packaging and staged scripts were prepared by Eric T. Perkey.

## Open items before publication

- [ ] **Locate the CIBERSORTx run outputs.** Cooper's Box project
      `CAR-T_FAS/` does **not** contain the CIBERSORTx deconvolution
      pipeline or its output CSVs — the stage-20 inputs
      (`stanford_cibersortx_cell_fractions.csv`,
      `stanford_cibersortx_bcell_gep.csv`) need to be sourced from Joanna,
      Nicholas Phillips, or a Brody-lab machine. The stage script's
      expected column layout (`Mixture`/`sample`/`Sample` + `B.cell|Memory.B|
      Naive.B|Plasma` fraction columns, and a `gene`×`sample` GEP CSV) is a
      best-inference layout and must be verified against the real export.
- [ ] **Confirm the Sup 1E panel content.** The capsule currently assembles
      Sup 1E as a side-by-side wrap of Fig 1H/1I; the real panel content
      (per the paper figure legend) may differ.
- [ ] **Verify `fun_process_rnaseq_deseq2` output filenames** against the
      Cooper helper's actual behavior if a reader regenerates intermediates
      from raw dbGaP counts (stage `00_prepare_data.R`).
- [ ] Brody-lab sign-off on capsule DOI and author order
- [ ] Final author list and ORCID iDs from submitted manuscript
- [ ] Replace "DOI pending" in the Code Availability Statement with the
      Code Ocean DOI once the capsule is published
- [ ] Confirm CIBERSORTx GEP matrices can be redistributed under the
      LM22 non-commercial research license
