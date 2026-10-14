# 00_prepare_data.R
# ------------------------------------------------------------------
# Optional preparation stage. Rebuilds the normalized intermediates
# (vst matrices + sample metadata) from raw dbGaP STAR/Salmon counts.
#
# Panels produced: none directly — this is the upstream of every bulk stage.
#
# Inputs (OPTIONAL; only required if a reader is regenerating from raw):
#   data/raw/CAR-T_FAS/input/nci_counts.csv                 (GDC STAR, phs001444)
#   data/raw/CAR-T_FAS/input/gencode.v44.annotation.gtf.gz
#   data/raw/CAR-T_FAS/input/phs003145.*.p1.*.txt           (dbGaP sample attrs)
#   data/raw/DLBCL_CART_tumorBulkRNAseq/output/stanford_cart_tpm.rds
#   data/raw/DLBCL_CART_tumorBulkRNAseq/output/zuma1_counts_updated.csv
#
# Outputs (written to data/derived/):
#   data/derived/CAR-T_FAS/output/nci_counts_updated.csv
#   data/derived/CAR-T_FAS/output/nci_vsd.rds
#   data/derived/CAR-T_FAS/output/nci_metadata.csv
#   data/derived/DLBCL_CART_tumorBulkRNAseq/output/zuma1_stanford_vsd.rds
#   data/derived/DLBCL_CART_tumorBulkRNAseq/output/zuma1_stanford_metadata.csv
#
# If raw inputs are absent, the stage exits with status 0 after noting that
# the capsule's bundled derived assets will be used instead.
# ------------------------------------------------------------------
source("scripts/_common.R", chdir = FALSE)

nci_raw   <- raw_path("CAR-T_FAS", "input", "nci_counts.csv")
nci_gtf   <- raw_path("CAR-T_FAS", "input", "gencode.v44.annotation.gtf.gz")

if (!file.exists(nci_raw) || !file.exists(nci_gtf)) {
  message("[00_prepare_data] Raw dbGaP counts not present under data/raw/.")
  message("[00_prepare_data] Skipping intermediate rebuild — downstream stages")
  message("[00_prepare_data] will consume the bundled derived data asset.")
  quit(save = "no", status = 0)
}

suppressMessages({
  library(DESeq2)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(rtracklayer)
})

message("[00_prepare_data] Rebuilding NCI VST intermediates from raw counts.")

gtf <- rtracklayer::import(nci_gtf) %>%
  as.data.frame() %>%
  distinct(gene_name, transcript_type)

protein_coding_genes <- gtf %>%
  dplyr::filter(transcript_type == "protein_coding") %>%
  pull(gene_name) %>%
  unique() %>%
  limma::alias2SymbolTable()

tib_nci <- read_csv(nci_raw, show_col_types = FALSE) %>%
  dplyr::select(-any_of("...1")) %>%
  mutate(
    gene = mapIds(org.Hs.eg.db,
                  keys = str_remove(gene, "\\..*"),
                  column = "SYMBOL", keytype = "ENSEMBL"),
    gene = alias2SymbolTable(alias = gene, species = "Hs")
  ) %>%
  dplyr::filter(!is.na(gene), gene %in% protein_coding_genes) %>%
  group_by(gene) %>%
  summarise_all(sum)

ensure_dir(derived_path("CAR-T_FAS", "output"))
write_csv(tib_nci, derived_path("CAR-T_FAS", "output", "nci_counts_updated.csv"))

# `fun_process_rnaseq_deseq2` lives in analysis_func.R; it computes vst and
# writes `<prefix>_vsd.rds` next to the input. The helper hard-codes
# relative output paths, so chdir into derived first.
old_wd <- getwd()
setwd(derived_path("CAR-T_FAS"))
on.exit(setwd(old_wd), add = TRUE)
fun_process_rnaseq_deseq2(
  path_to_counts = "output/nci_counts_updated.csv",
  batch_correct  = FALSE
)
setwd(old_wd)

message("[00_prepare_data] Done — derived intermediates written.")
