# Functions for executing analysis.Rmd
# Created 3/5/24 by Alan Cooper


##### Packages #####
suppressMessages(library(readxl))
suppressMessages(library(tidyverse))
suppressMessages(library(vroom))

# Reference genomes and working with gene IDs
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))

# Process RNA-seq
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))
suppressMessages(library(GEOquery))
suppressMessages(library(GenomicFeatures))

# GSEA and GSVA
suppressMessages(library(GSVA))
# clusterProfiler and fgsea are only used by Cooper's GSEA wrapper (not by any
# Fig 1 / Sup 1 panel this capsule reproduces). Loaded only if available so
# the capsule doesn't require them.
suppressWarnings(suppressMessages(require(clusterProfiler, quietly = TRUE)))
suppressWarnings(suppressMessages(require(fgsea, quietly = TRUE)))

# Statistics (other than base stats)
suppressMessages(library(rstatix))


# For converting probes to genes
map_probes2genes <- function(exprs_mat, gene_ids = NULL, data_type = NULL, probe_db = NULL){
  # Remove duplicate probes and map to HUGO IDs
  # Choosing probe with highest median expression
  # Mapping IDs based on corresponding database
  # Also removes probes without genes 
  
  gene_filter <- function(exprs_mat, genes) {
    ranks <- apply(exprs_mat, 1, median)
    df_temp <- data.frame(
      genes = genes,
      row = seq_along(ranks),
      ranks = ranks
    )

    df_temp <- df_temp[order(df_temp$ranks, decreasing = TRUE), ]
    keep <- df_temp[!duplicated(df_temp$genes) & !is.na(df_temp$genes),]$row # !is.na() removes probes without genes

    exprs_mat_return <- exprs_mat[keep, ]
    rownames(exprs_mat_return) <- genes[keep]

    return(exprs_mat_return)
  }

  # Use input vector of genes
  if (is.null(gene_ids)) {
    
    # if (!is.null(data_type)) {
    #   # Map gene symbols to probes based on platform database
    #   if (data_type == "hgu133plus") {
    #     library(hgu133plus2.db)
    #     probe_db <- "hgu133plus2.db"
    #     
    #   } else if (data_type == "illumina") {
    #     library(illuminaHumanv4.db)
    #     probe_db <- "illuminaHumanv4.db"
    #     
    #   } else if (data_type == "hta2") {
    #     library(hta20transcriptcluster.db)
    #     probe_db <- "hta20transcriptcluster.db"
    #     
    #   } else if (!is.null(probe_db)) {
    #   
    #   library(probe_db, character.only = TRUE)
    #   
    #   }else{
    #   
    #   message("Must specificy probe_db")
    #   }
    # }
    
    gene_ids <- mapIds(eval(parse(text = probe_db)), keys = rownames(exprs_mat), column = "SYMBOL", keytype = "PROBEID")
    
  }
  
  gene_ids <- alias2SymbolTable(gene_ids, species = "Hs") # update gene name
  exprs_mat <- gene_filter(exprs_mat, gene_ids)

  return(exprs_mat)

}

download_from_geo <- function(gse, data_title = NULL, download_dir, download_supp = FALSE) {
  # Title for data set
  if (is.null(data_title)) {
    data_title <- gse
  }
  
  geo_dir <- file.path(download_dir, data_title)
  
  # Make directory if does not exist
  if (sum(grepl(geo_dir, list.dirs(path = download_dir, full.names = TRUE, recursive = FALSE))) == 0) {
    dir.create(geo_dir)
  }
  
  if (download_supp) {
    # Download raw data
    # This makes a new directory that getGEO() will deposit into 
    getGEOSuppFiles(gse, makeDirectory = TRUE, baseDir = download_dir)
    
    message(
      paste("Raw data is in:", paste(file.path(geo_dir, list.files(geo_dir)), collapse = " ")))
    
    # Download GEO data set
    # Only using for metadata
    gse_es <- getGEO(gse, AnnotGPL = FALSE, destdir = geo_dir, getGPL = FALSE)[[1]] 
  }else{
    gse_es <- getGEO(gse, AnnotGPL = TRUE, destdir = geo_dir)[[1]]
  }
  
  # Write out metadata
  pData(gse_es) %>% 
    write_csv(get_file("input", sprintf("%s_metadata.csv", data_title)))
  
  return(gse_es)
}

get_geo_microarray <- function(gse, data_title = NULL, download_dir, download_supp = FALSE, process_gse = TRUE, data_type = NULL, probe_db = NULL, use_gpl = FALSE, probe_col = "ID", gene_col = "Gene Symbol") {
  # data_type - one of "hgu133plus",  "illumina", or "hta2". Used for choosing hgu133plus2.db, illuminaHumanv4.db, or hta20transcriptcluster.db for annotations
  
  # Title for data set
  if (is.null(data_title)) {
    data_title <- gse
  }
  geo_dir <- file.path(download_dir, data_title)
  
  # Make directory if does not exist
  if (sum(grepl(geo_dir, list.dirs(path = download_dir, full.names = TRUE, recursive = FALSE))) == 0) {
    dir.create(geo_dir)
  }
  
  if (download_supp) {
    # Download raw data
    # This makes a new directory that getGEO() will deposit into 
    getGEOSuppFiles(gse, makeDirectory = TRUE, baseDir = download_dir)
    
    message(
      paste("Raw data is in:", paste(file.path("\n", geo_dir, list.files(geo_dir)), collapse = " ")))
    
    # Download GEO data set
    # Only using for metadata
    gse_es <- getGEO(gse, AnnotGPL = FALSE, destdir = geo_dir, getGPL = FALSE)[[1]] 
  }

  if (!download_supp & process_gse) {
    
    gse_es <- getGEO(gse, AnnotGPL = TRUE, destdir = geo_dir)[[1]]
  
    
    # Remove duplicate probes and map to HUGO IDs
    # Choosing probe with highest median expression
    # Mapping IDs based on corresponding database
    gene_filter <- function(es, genes) {
      ranks <- apply(exprs(es), 1, median)
      df_temp <- data.frame(
        genes = genes,
        row = seq_along(ranks),
        ranks = ranks
      )

      df_temp <- df_temp[order(df_temp$ranks, decreasing = TRUE), ]
      keep <- df_temp[!duplicated(df_temp$genes) & !is.na(df_temp$genes),]$row

      es_return <- es[keep, ]
      rownames(es_return) <- genes[keep]

      return(exprs(es_return))
    }
    
    
    exprs_mat <- exprs(gse_es)
  
    # Annotation for probes to genes
    if (use_gpl) {
      message(paste("GPL columns:", paste(colnames(fData(gse_es)), collapse = " ")))
      
      gpl_keys <- fData(gse_es) %>% 
        as.data.frame() %>% 
        dplyr::rename(gene = !!as.symbol(gene_col), probe = !!as.symbol(probe_col)) %>% 
        pull(gene, probe)
      
      gene_ids <- unname(gpl_keys)[match(names(gpl_keys), rownames(gse_es))]
      # gene_ids <- gene_ids[gene_ids != ""]
      gene_ids <- alias2SymbolTable(gene_ids, species = "Hs")
      
    }else{
      if (is.null(data_type)) {
        
        if (!is.null(probe_db)) {
          library(probe_db, character.only = TRUE)
          
        } else {
          message("Must specificy probe_db")
        }
        
      }else{
        
        if (data_type == "hgu133plus") {
          library(hgu133plus2.db)
          probe_db <- "hgu133plus2.db"
          
        } else if (data_type == "illumina") {
          library(illuminaHumanv4.db)
          probe_db <- "illuminaHumanv4.db"
          
        } else if (data_type == "hta2") {
          library(hta20transcriptcluster.db)
          probe_db <- "hta20transcriptcluster.db"
          
        }
      }
      gene_ids <- mapIds(eval(parse(text = probe_db)), keys = rownames(exprs_mat), column = "SYMBOL", keytype = "PROBEID")
      gene_ids <- alias2SymbolTable(gene_ids, species = "Hs")
    }

    exprs_mat <- gene_filter(gse_es, gene_ids)
    
    # Avereps - an alternative method for condensing probes to genes
    # exprs_mat <- avereps(exprs(gse_es), gene_ids)
  
    # Write out processed expression matrix
    exprs_mat %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "gene") %>% 
      write_csv(get_file("output", sprintf("%s_microarray.csv", data_title)))
    
    # Write out unprocessed expression matrix
    exprs_mat %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "gene") %>% 
      write_csv(get_file("input", sprintf("%s_microarray_unprocessed.csv", data_title)))
    
    message(
      paste("Processed data is in:", paste(get_file("output", sprintf("%s_microarray.csv", data_title)))))
    
  }
  
  # Write out metadata
  pData(gse_es) %>% 
    write_csv(get_file("input", sprintf("%s_metadata.csv", data_title)))
  
}

# Process RNA-seq
fun_dgelist <- function(path_to_counts, is_ensembl = TRUE, batch_correct = FALSE, output = NULL){
  # Process RNA-seq count matrix by:
  # 1. Converting from ENSEMBL to HUGO gene names (if necessary)
  # 2. Filtering out genes with low counts
  # 3. TMM-normalizing
  # 4. Converting to LCPM (saves this as separate .csv file to output folder)
  # 5. Removing batch effect (if necessary)
  # Saves final DGEList to output folder
  
  # Parameters:
  #   path_to_counts - path name to count matrix .csv file (e.g. "input/counts.csv"). 
  #                - expects column in count matrix for gene IDs to be named "gene"
  #   is_ensembl - if TRUE, converts genes to HUGO IDs using mapIds() with org.Hs.eg.db
  #   batch_correct - if TRUE, performs batch correction using removeBatchEffect() from limma
  #               - batches are extracted from sample names (e.g. sample B1-DLBCL001 is from batch B1)
  #   output - name for final DGEList saved as .rds file in output folder
  #          - defaults to base name of path_to_counts but removes "_counts.csv"
  
  # Convert gene count tibble to DGEList
  fun_tib2dge <- function(tib_data, rownames,
                          samples = NULL, group = NULL) {
    mat_data <- tib_data %>%
      dplyr::select(-all_of(rownames)) %>%
      data.matrix()
    rownames(mat_data) <- dplyr::pull(tib_data, all_of(rownames))
    mat_data %>%
      DGEList(samples = samples, group = group) %>%
      return()
  }
  
  # Transposing the tibble
  fun_t_tib <- function(tib_data) {
    tib_data %>%
      gather(sample, value, 2:ncol(.)) %>%
      spread(gene, value) %>%
      return()
  }
  
  # Read in data and convert gene names if necessary
  if (is_ensembl) {
    tib_data <- read_csv(get_file(path_to_counts)) %>%
      mutate(
        gene = mapIds(
          org.Hs.eg.db,
          keys = str_remove(gene, "\\..*"),
          column = "SYMBOL",
          keytype = "ENSEMBL"
        ),
        gene = alias2SymbolTable(alias = gene, species = "Hs") 
      ) %>%
      group_by(gene) %>%
      summarise_all(sum)
  }else{
    tib_data <- read_csv(get_file(path_to_counts)) %>% 
      group_by(gene) %>%
      summarise_all(sum)
  }
  
  # Make DGEList
  tib_samples <- data.frame(
    row.names = colnames(tib_data %>% dplyr::select(-gene))
  )
  
  dge_pre <- tib_data %>%
    fun_tib2dge("gene", tib_samples, tib_samples$experiment)
  
  rm(tib_data)
  
  
  # Run LCPM
  fun_lcpm_csv <- function(lcpm, str_name) {
    bind_cols(
      gene = rownames(lcpm),
      lcpm %>% as_tibble()
    ) %>%
      fun_t_tib() %>%
      write_csv(get_file("output", sprintf("%s_lcpm_normalized.csv", str_name)))
  }
  
  fun_lcpm_density <- function(dge = dge_pre, name,
                               lcpm = NULL, plot = TRUE) {
    if (is.null(lcpm)) {
      lcpm <- dge %>%
        cpm(normalized.lib.sizes = TRUE,
            log = TRUE)
    }
    if (plot) {
      plot_density <- lcpm %>%
        as_tibble() %>%
        gather(key = "Sample", value = "Intensity") %>%
        ggplot(aes(Intensity, group = Sample)) +
        geom_density() +
        theme(legend.position = "none") +
        ggtitle(paste0("Density plot of ", name, " logCPM")) +
        xlim(c(NA, 25))
      plot_density %>%
        print()
    }
    return(lcpm)
  }
  
  ## If bool_plot = TRUE, makes density plots and saves pre- and post-normalized counts
  fun_save_dge <- function(dge_pre, str_name, bool_plot) {
    if (bool_plot) {
      lcpm_pre <- dge_pre %>%
        fun_lcpm_density(paste0(str_name, " samples, pre-normalization"), plot = bool_plot)
      rm(lcpm_pre)
    }
    dge_post <- calcNormFactors(dge_pre, method = "TMM")
    lcpm_post <- dge_post %>%
      fun_lcpm_density(paste0(str_name, " samples, post-normalization"), plot = bool_plot)
    dge_filt <- dge_post[apply(lcpm_post, 1, mean) >= 0.5,] # cutoff from GSVA paper
    rm(lcpm_post)
    dge_filt %>%
      fun_lcpm_density(paste0(str_name, " samples, filtered"), plot = bool_plot) %>%
      fun_lcpm_csv(str_name) 
    return(dge_filt)
  }
  
  fun_lcpm_final <- function(dge, batch_correct) {
    lcpm <- cpm(dge, log = TRUE)
    
    if (batch_correct) {
      batches <- sapply(
        colnames(dge),
        function(str_name) str_split(str_name, "-")[[1]][1]
      )
      
      if (length(unique(batches)) > 1)
        message(
          paste(
            sprintf("Removing batch effect for %s batches:", length(unique(batches))), 
            paste(unique(batches), collapse = " ")
          ))
      lcpm <- removeBatchEffect(lcpm, batches) 
    }
    return(lcpm)
  }
  
  data_name <- basename(path_to_counts) %>% str_remove_all("_counts.+csv")
  dge_post <- fun_save_dge(dge_pre, data_name, FALSE) 
  
  saveRDS(dge_post, get_file("output", sprintf("dge_tmm_%s.rds", data_name))) # this will be used for limma voom
  
  dge_post$counts <- fun_lcpm_final(dge_post, batch_correct)
  
  rm(dge_pre)
  
  if (is.null(output)) {
    output <- get_file("output", paste0("dge_", data_name,".rds")) 
  }else{
    output <- get_file("output", paste0(output,".rds")) 
  }
  saveRDS(dge_post, output)
  
}

fun_process_rnaseq_deseq2 <- function(path_to_counts,path_to_metadata=NULL, batch_correct = FALSE, bool_plot = TRUE, batch_column = "batch"){
  require(vsn)
  require(DESeq2)
  require(tidyverse)
  # Process RNA-seq count matrix using DESeq2
  
  # Parameters:
  #   path_to_counts - path name to count matrix .csv file (e.g. "input/counts.csv"). 
  #                - expects column in count matrix for gene IDs to be named "gene"
  #   batch_correct - if TRUE, performs batch correction using removeBatchEffect() from limma
  #               - batches are extracted from sample names (e.g. sample B1-DLBCL001 is from batch B1)
  #   output - name for final DGEList saved as .rds file in output folder
  #          - defaults to base name of path_to_counts but removes "_counts.csv"
  
  counts_mat <- read_csv(get_file(path_to_counts)) %>%
    column_to_rownames(var = "gene") %>% 
    as.matrix()
  
  data_name <- basename(path_to_counts) %>% str_remove_all("_counts.+csv")
  
  mode(counts_mat) <- "integer"

  if (!is.null(path_to_metadata)){
    metadata <- read_csv(get_file(path_to_metadata)) %>% 
      dplyr::filter(sample %in% colnames(counts_mat)) %>% 
      mutate(
        # batch = "B1",
        batch = !!as.symbol(batch_column)
        ) %>% 
      arrange(factor(sample, colnames(counts_mat))) %>% 
      column_to_rownames(var = "sample")
    print(table(metadata$batch))
  }else{
    metadata <- data.frame(sample = colnames(counts_mat)) %>% 
      arrange(factor(sample, colnames(counts_mat)))
    rownames(metadata) <- metadata$sample
    
    if(batch_correct){
      metadata <- metadata %>% 
        mutate(batch = str_remove(sample,"-.+"))
      
      print(table(metadata$batch))
    }
  }
  
  
  # # Filter out low counts
  # keep <- apply(counts_mat, 1, median)  > 20
  # counts_mat <- counts_mat[keep,]
  # dim(counts_mat)
  
  # browser()
  
  if (batch_correct) {
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_mat, 
      colData = metadata,
      design = ~batch)
    
    pre_vsd <- vst(dds, blind = FALSE) # batch-aware transformation
    
    # Batch correction
    vsd <- pre_vsd
    assay(vsd) <- removeBatchEffect(assay(pre_vsd), batch = colData(dds)$batch)
    
    if (bool_plot){
      
      plot_pre <- plotPCA(pre_vsd, intgroup = c("batch")) +
        ggtitle("Before batch correction")
      
      plot_post <- plotPCA(vsd, intgroup = c("batch")) +
        ggtitle("After batch correction")
      
      print(plot_pre + plot_post)
    }
    
  }else{
    dds <- DESeqDataSetFromMatrix(
      countData = counts_mat, 
      colData = metadata, 
      design = ~1)
    
    vsd <- vst(dds)
  }
  
  if (bool_plot){
    pre_vsd_plot <- (meanSdPlot(assay(dds))$gg) +
      ggtitle("Pre-vst")
    
    post_vsd_plot <- (meanSdPlot(assay(vsd))$gg) +
      ggtitle("Post-vst")
    
    combined <- pre_vsd_plot + post_vsd_plot
    print(combined)
    
    if (!batch_correct) {
      pca_plot <- autoplot(prcomp(t(assay(vsd))))
      print(pca_plot)
    }
    
    # # Density plots pre- and post-filtering
    # tmp <- calcNormFactors(dge_pre, "TMM")
    # pre_mat <- cpm(tmp, log = TRUE)
    # 
    # post_mat <- dge_lcpm$counts
    # 
    # 
    # pre_density <- pre_mat %>%
    #     as_tibble() %>%
    #     pivot_longer(cols = everything(), names_to = "Sample", values_to = "Intensity") %>% 
    #     ggplot(aes(Intensity, group = Sample)) +
    #     geom_density() +
    #     theme(legend.position = "none") +
    #     ggtitle(paste0("Pre-filtered logCPM")) +
    #     xlim(c(NA, 25))
    # 
    # post_density <- post_mat %>%
    #     as_tibble() %>%
    #     pivot_longer(cols = everything(), names_to = "Sample", values_to = "Intensity") %>% 
    #     ggplot(aes(Intensity, group = Sample)) +
    #     geom_density() +
    #     theme(legend.position = "none") +
    #     ggtitle(paste0("Post-filtered logCPM")) +
    #     xlim(c(NA, 25))
    # 
    # print(pre_density + post_density)
    
  }
  
  
  message(paste("Saving:", get_file("output", paste0(data_name,"_vsd", ".rds"))))
  saveRDS(vsd,
          get_file("output", paste0(data_name,"_vsd", ".rds")))
  
}

##### DGE and GSEA ####
compare_dge <- function(v, contrast = NULL, toptable = TRUE, pval_threshold = 0.1, log2fc_threshold = 0.5, decideTests_method = "global", coef = NULL){
  
  # Contrast can be charcacter or natrix
  
  # CONSIDER USING treat() and topTreat
  # https://support.bioconductor.org/p/9139092/
  # https://support.bioconductor.org/p/65635/#65638
  # treat() produces unshrunk fold change
  
  design <- v$design
  
  if (is.character(contrast)) {
    contrasts_mat <- makeContrasts(contrasts = contrast, levels = colnames(design))
  }else{
    contrasts_mat <- contrast
  }
  

  vfit <- lmFit(v, design)
  if (!is.null(coef)) {
    vfit <- contrasts.fit(vfit, coefficients = coef) 
  }else{
    vfit <- contrasts.fit(vfit, contrasts = contrasts_mat) 
  }
  efit <- eBayes(vfit)

  
  if (toptable) {
    out <- topTable(efit, adjust.method = "BH", sort.by = "none", number = Inf)
  }else{
    out <- decideTests(efit, adjust.method = "BH", p.value = pval_threshold, lfc = log2fc_threshold, method = decideTests_method)
  }
  
  
  return(out) 
}

fun_gsea <- function(dge_data, v, contrast = NULL, geneset_df, dge_type = "limma", column_to_rank_on = "signed_pval"){
  # Run pre-ranked GSEA   
  #   column_to_rank_on - one of "signed_pval" or "t_stat"
  
  
  # Run contrast for results from limma::topTable()
  if (dge_type == "limma") {
    
    # Run contrast if don't supply output from topTable()
    if (!is.null(contrast)) {
      dge_data <- compare_dge(v, contrast)
    }
  
    logfc_col <- "logFC"
    tstat_col <-  "t"
    pval_col <- "P.Value"
    adjpval_col <- "adj.P.Val"
    
  # Expects data from DESeq2::results()
  } else if (dge_type == "DESeq2") {
    
    logfc_col <- "log2FoldChange"
    tstat_col <- "stat"
    pval_col <- "pvalue"
    adjpval_col <- "padj"
  }
  
  if (column_to_rank_on == "t_stat") {
    column_to_rank_on <- tstat_col
    
  } else if (column_to_rank_on == "logfc") {
    logfc_col <- logfc_col
    
  }
  
  message(paste("Ranking genes on", column_to_rank_on))
  
  # Rank genes based on sign(log2FC) * -log10(p-value)
  gene_rank_df <- dge_data  %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    mutate(signed_pval = sign(!!as.symbol(logfc_col))*-log10(!!as.symbol(pval_col)))

  gene_rank_vec <- gene_rank_df %>% pull(!!as.symbol(column_to_rank_on))
  names(gene_rank_vec) <- gene_rank_df$gene
  gene_rank_vec <-  sort(gene_rank_vec, decreasing = T) # sort

  
  gsea <- GSEA(
    geneList = gene_rank_vec, 
    by = "fgsea",
    nPermSimple = 10000, # this prevents an error which results in NA for some gene sets; an arg to fgseaMultilevel()
    TERM2GENE = geneset_df %>% dplyr::select(gs_name, geneID = gene_symbol),
    pAdjustMethod = "BH", 
    pvalueCutoff = 1 # return all results
  )
  
  return(gsea) 
}

fun_files_for_gsea <- function(exprs_mat_filepath, pheno_df, pheno, out_file_name = NULL, batch_correct = FALSE, rnaseq = TRUE, method = "lcpm", deseq2_design = "~batch", file_tag = "", out_dir = NULL){
  # Prepare files for standard GSEA to run on Java client from https://www.gsea-msigdb.org/gsea/downloads.jsp 
  # Make .txt file for gene expression
  # Make .cls file for phenotype
  
  # For `rnaseq` = TRUE, `method` is one of "lcpm" (TMM-normalization + log2-CPM) or "vst" (`DESeq2::vst()`)
  
  if (is.null(out_dir)) {
    out_dir <- dirname(exprs_mat_filepath)
  }
  
  # Read in expression data
  exprs_mat <- read_csv(exprs_mat_filepath) %>% 
    tibble::column_to_rownames(var = "gene") %>% 
    as.matrix()
  
  pheno_vec <- pheno_df %>%
    pull(!!as.symbol(pheno), name = "sample") %>% 
    factor() %>% 
    sort() # makes constructing phenotype file easier
  
  print(
    sprintf("The expression matrix has %s genes and %s samples", nrow(exprs_mat), ncol(exprs_mat))
  )
  
  # Make file for gene expression
  # RNA-seq (from count matrix)
  #  - Counts are either
  #       - TMM-normalized and log2CPM-transformed, or
  #       - normalized with VST (will use `deseq2_design` for design)
  #  - Genes with low counts are not filtered out, per GSEA recommendations
  #  - Batch correction is performed if `batch_correct = TRUE`
  # Microarray
  # - assume it is already normalized
  
  if (rnaseq) {
    
    if (method == "lcpm") {
      dge <- DGEList(exprs_mat)
      dge <- calcNormFactors(dge, method = "TMM")
      lcpm <- cpm(dge, log = TRUE)
      
      if (batch_correct) {
        batches <- sapply(
          colnames(dge),
          function(str_name) str_split(str_name, "-")[[1]][1]
        )
        
        if (length(unique(batches)) > 1) {
          message(
            paste(
              sprintf("Removing batch effect for %s batches:", length(unique(batches))),
              paste(unique(batches), collapse = " ")
            ))
          lcpm <- removeBatchEffect(lcpm, batches)
        }
      }
      
      exprs_mat <- lcpm
      
    } else if (method == "vst") {
      
      mode(exprs_mat) <- "integer"
      
      dds_metadata_df <- pheno_df %>% 
        column_to_rownames(var = "sample")
      
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_mat[,rownames(dds_metadata_df)], 
        colData = dds_metadata_df,
        design = as.formula(deseq2_design))
      
      vsd <- vst(dds)
      
      # Batch correction
      if (batch_correct) {
        batches <- sapply(
          colnames(vsd),
          function(str_name) str_split(str_name, "-")[[1]][1]
        )
        
        if (length(unique(batches)) > 1) {
          message(
            paste(
              sprintf("Removing batch effect for %s batches:", length(unique(batches))),
              paste(unique(batches), collapse = " ")
            ))
          assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch = batches)
        }
      }
      
      exprs_mat <- assay(vsd)
    }
    
  }
  
  # File names for write out
  if (is.null(out_file_name)) {
    
    # Extract name of data and build new file path
    exprs_out_basename <- basename(exprs_mat_filepath) %>% 
      str_remove_all("_counts.+csv") %>% 
      str_remove_all("_microarray.+csv") %>% 
      str_c("_exprs_for_standard_gsea.txt")
    
    exprs_out_file_name <- paste0(out_dir, "/", file_tag, "_", exprs_out_basename)
    
    pheno_out_basename <- basename(exprs_mat_filepath) %>% 
      str_remove_all("_counts.+csv") %>% 
      str_c("_pheno_for_standard_gsea.cls")
    
    pheno_out_file_name <- paste0(out_dir, "/", file_tag, "_", pheno_out_basename)
  }
  
  # Write out expression file
  message(paste("Writing out to", exprs_out_file_name))
  
  
  exprs_mat <- exprs_mat[, match(names(pheno_vec), colnames(exprs_mat))] # order to match phenotype list
  
  exprs_mat %>%
    as.data.frame() %>%
    rownames_to_column(var = "NAME") %>%
    mutate(DESCRIPTION = NA_character_) %>%
    dplyr::select(NAME, DESCRIPTION, everything()) %>%
    write_delim(file = exprs_out_file_name, delim = "\t")

  # Make phenotype file
  # Inspired by https://rdrr.io/bioc/ArrayTools/src/R/output.cls.R
  
  message(paste("Writing out to", pheno_out_file_name))
  
  line1 <- paste(length(pheno_vec), length(unique(pheno_vec)), 1) # N N-classes 1
  line2 <- c("# ", paste(unique(pheno_vec),  collapse = " ")) # class labels
  line3 <- paste(pheno_vec, collapse = " ") # class labels for each sample
  
  cat(line1, "\n", sep = "", file = pheno_out_file_name)
  cat(line2, "\n", sep = "", file = pheno_out_file_name, append = TRUE)
  cat(line3, "\n", sep = "", file = pheno_out_file_name, append = TRUE)
  
  # Print to console
  message("Preview of .cls file:")

  print(line1)
  print(line2)
  print(substr(line3, 1, 20))
  
  
}

##### GSVA Pipeline #####

fun_gsva <- function(lcpm, genesets) {
  # Run GSVA
  
  # Parameters:
  #   lcpm - an LCPM-transformed count matrix
  #   genesets - a dataframe of genesets to be used for GSVA
  #            - expects "gene" column with HUGO gene IDs and "gs" column with the name of the corresponding geneset 
  
  # Returns:
  #   gsva_dge - a matrix of GSVA scores by sample
  
  vec_gs <- unique(genesets$gs)
  gene_sets <- lapply(vec_gs,
                      function(gs1) genesets %>%
                        filter(gs == gs1) %>%
                        pull(gene))
  names(gene_sets) <- vec_gs

  # GSVA 1.50+ (Bioc 3.17+) replaced direct-argument gsva() with parameter
  # objects. Branch on version so this code runs on both GSVA 1.46 (R 4.2.x,
  # Bioc 3.16 — the manuscript's pin) and GSVA 2.0+ (R 4.4.2, Bioc 3.20 —
  # the Code Ocean environment).
  if (utils::packageVersion("GSVA") >= "1.50") {
    gp <- GSVA::gsvaParam(exprData = lcpm, geneSets = gene_sets, kcdf = "Gaussian")
    gsva_dge <- GSVA::gsva(gp, verbose = TRUE)
  } else {
    gsva_dge <- gsva(
      expr = lcpm,
      kcdf = "Gaussian",
      gset.idx.list = gene_sets,
      annotation = NULL,
      method = "gsva",
      verbose = TRUE
    )
  }
  
  # Rename matrix of GSVA output to include gene set names
  if (nrow(gsva_dge) == length(vec_gs)) {
    rownames(gsva_dge) <- vec_gs
  }else{
    
    # Gene sets not included in GSVA output because their genes are not include in the data
    gsets_excluded <-  genesets %>% 
      group_by(gs) %>%
      summarize(gene = list(gene)) %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(sum = sum(c(unlist(gene) %in% rownames(lcpm)))/length(unlist(gene))) %>% 
      dplyr::filter(sum == 0) %>% 
      pull(gs)
    
    # Rename with appropriate columns
    rownames(gsva_dge) <- vec_gs[!vec_gs %in% gsets_excluded]
  }
  
  
  return(gsva_dge)
  
} 

fun_pca <- function(gsva_dge, genesets = tib_genesets){
  # Run PCA on GSVA scores
  # Generally, PC1 is the "immune axis," and PC2 is the "cell of origin axis"
  # Rotates the output such that samples with higher immune scores are on the right and samples with higher ABC-DLBCL scores are on the top
  
  # Parameters:
  #   gsva_dge - output from fun_gsva()
  #   genesets - a dataframe of genesets as used in fun_gsva()
  
  # Returns:
  #   pca_dge - an object of class "prcomp" where PCA scores are located in the "x" matrix
  
  vec_gs <- unique(genesets$gs)
  
  pca_dge <- prcomp(t(gsva_dge))
  
  rownames(pca_dge$rotation) <- rownames(gsva_dge)
  # rownames(pca_dge$rotation) <- vec_gs
  
  
  if (pca_dge$rotation["ifng","PC1"] < 0) {
    pca_dge$x[,"PC1"] <- -pca_dge$x[,"PC1"]
    pca_dge$rotation[,"PC1"] <- -pca_dge$rotation[,"PC1"]
  }
  
  if (pca_dge$rotation["ABCDLBCL-1","PC2"] < 0) {
    pca_dge$x[,"PC2"] <- -pca_dge$x[,"PC2"]
    pca_dge$rotation[,"PC2"] <- -pca_dge$rotation[,"PC2"]
  }
  
  return(pca_dge)
}  


fun_kmeans <- function(pca_dge, k = 4, bool_plot = FALSE){  
  # Perform k-means clustering on PCA scores
  # By default, tries to define 4 clusters:
  #   1. ABC Cold - cluster with low PC1 scores and high PC2 scores
  #   2. ABC Hot - cluster with high PC1 scores and high PC2 scores
  #   3. GCB Cold - cluster with low PC1 scores and low PC2 scores
  #   4. GCB Hot - cluster with high PC1 scores and low PC2 scores
  
  # Parameters:
  #   pca_dge - output from fun_pca()
  #   k - number of clusters to partition data into
  #   seed - seed for reproducible results
  #   bool_plot - if TRUE, makes a PCA plot and a PCA biplot using autoplot()
  
  # Returns a list with:
  #   kmeans - a dataframe with the cluster assignments of each sample
  #   biplot - PCA biplot (optional)
  #   pcaplot - PCA plot (optional)
  
  set.seed(900)
  kmeans_dge <- kmeans(pca_dge$x, k, nstart = 50)
  
  axis_1 <- c("Cold", "Hot")
  axis_2 <- c("GCB", "ABC")
  
  fun_label_cluster <- function(clust) {
    axis_1_mean <- mean(pca_dge$x[kmeans_dge$cluster == clust,"PC1"])
    axis_2_mean <- mean(pca_dge$x[kmeans_dge$cluster == clust,"PC2"])
    return(paste(axis_2[(axis_2_mean > 0) + 1],
                 axis_1[(axis_1_mean > 0) + 1],
                 sep = " "))
  }
  
  df_kmeans <- data.frame(
    row.names = names(kmeans_dge$cluster),
    cluster = sapply(1:k, fun_label_cluster)[kmeans_dge$cluster]
  )
  
  if (bool_plot) {
    plot_biplot <- autoplot(
      pca_dge, data = df_kmeans, colour = "cluster",
      loadings = TRUE, loadings.label = TRUE, frame = FALSE, scale = 0,
    )
    plot_pca <- autoplot(
      pca_dge, data = df_kmeans, colour = "cluster", frame = FALSE, scale = 0,
    )
  }else{
    plot_biplot <- NULL
    plot_pca <- NULL
  }
  
  return(
    list(
      kmeans = df_kmeans,
      biplot = plot_biplot,
      pcaplot = plot_pca
    )
  )
}

