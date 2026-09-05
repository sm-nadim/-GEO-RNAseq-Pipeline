#!/usr/bin/env Rscript
# ==============================================================================
# Modified RNA-seq Analysis Pipeline
# Modifications: 
# 1. Replaced DESeq2 with Limma (to handle pre-normalized logCPM/FPKM data).
# 2. Added dynamic group assignment (Control vs Treatment).
# 3. Handled multiple file formats (Single Matrix, Cuffdiff, HTSeq-style counts).
# ==============================================================================

options(timeout=300)
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(tidyverse)
  library(pheatmap)
  library(Biobase)
  library(pROC)
})

# List of all assigned 11 GEO Datasets
geo_ids <- c("GSE158266", "GSE246030", "GSE140323", "GSE154388",
             "GSE107839", "GSE119785", "GSE124412", "GSE141134", 
             "GSE164073", "GSE21642", "GSE158473")

for (geo_id in geo_ids) {
  cat("\n================================================================\n")
  cat("  ANALYZING:", geo_id, "\n")
  cat("================================================================\n")
  
  out_dir <- file.path(getwd(), geo_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    # Step 1: Download GEO dataset and Phenotype data
    gse <- getGEO(geo_id, GSEMatrix = TRUE)
    pheno <- pData(gse[[1]])
    expr_matrix <- exprs(gse[[1]])
    
    # --- Custom Data Parsers for specific complex data formats ---
    if (geo_id == "GSE158266") {
      file_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158266/suppl/GSE158266_results_TreatedvsControl_ovary_SVremoved.txt.gz"
      dest <- file.path(out_dir, basename(file_url))
      if(!file.exists(dest)) download.file(file_url, dest, mode="wb")
      df <- read.delim(gzfile(dest), check.names=FALSE)
      cpm_cols <- grep("logCPM\\.", colnames(df), value=TRUE)
      expr_matrix <- as.matrix(df[, cpm_cols])
      rownames(expr_matrix) <- df$SYMBOL
      colnames(expr_matrix) <- rownames(pheno)[1:length(cpm_cols)]
    } 
    else if (geo_id == "GSE107839") {
      file_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107839/suppl/GSE107839_genes.count_tracking.txt.gz"
      dest <- file.path(out_dir, basename(file_url))
      if(!file.exists(dest)) download.file(file_url, dest, mode="wb")
      df <- read.delim(gzfile(dest))
      count_cols <- grep("_count$", colnames(df), value=TRUE)
      expr_matrix <- as.matrix(df[, count_cols])
      rownames(expr_matrix) <- df$tracking_id
      pheno <- data.frame(title=count_cols, row.names=count_cols)
      pheno$condition <- factor(ifelse(grepl("CTL", count_cols), "Control", "Treatment"))
      pheno$binary_condition <- pheno$condition
      colnames(expr_matrix) <- rownames(pheno)
    } 
    else if (nrow(expr_matrix) == 0) {
      supp <- getGEOSuppFiles(geo_id, fetch_files = FALSE)
      if(!is.null(supp) && nrow(supp)>0) {
        for(i in 1:nrow(supp)) {
          if (nrow(expr_matrix)>0) break
          if(grepl("RAW\\.tar", supp$fname[i])) next
          dest <- file.path(out_dir, supp$fname[i])
          if(!file.exists(dest)) download.file(supp$url[i], dest, mode="wb")
          df <- tryCatch(read.delim(gzfile(dest), check.names=FALSE), error=function(e) NULL)
          if(!is.null(df) && ncol(df) > 3) {
            is_num <- sapply(df, is.numeric)
            if(sum(is_num) >= 2) {
              expr_matrix <- as.matrix(df[, is_num])
              if(!is_num[1]) rownames(expr_matrix) <- make.unique(as.character(df[[1]]))
            }
          }
        }
      }
    }
    
    if (nrow(expr_matrix) == 0) stop("Expression matrix could not be extracted or is empty.")
    
    # Step 2: Data Cleaning & Normalization check
    expr_matrix[is.na(expr_matrix)] <- 0
    if(max(expr_matrix, na.rm=TRUE) > 100) expr_matrix <- log2(expr_matrix + 1)
    expr_matrix <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
    
    # Step 3: Group Assignment (Control vs Treatment)
    if (geo_id != "GSE107839") {
      condition <- NULL
      char_cols <- grep("characteristics", colnames(pheno), value=TRUE)
      for (col in char_cols) {
        vals <- unique(pheno[[col]])
        if (length(vals) >= 2 && length(vals) <= 10) {
          condition <- factor(make.names(gsub("^[^:]+: ", "", pheno[[col]])))
          break
        }
      }
      if(is.null(condition)) condition <- factor(make.names(pheno$source_name_ch1))
      
      ref_group <- levels(condition)[1]
      for(pat in c("control", "normal", "vehicle", "mock", "wt", "ctl")) {
        m <- grep(pat, levels(condition), ignore.case=TRUE, value=TRUE)
        if(length(m)>0) { ref_group <- m[1]; break }
      }
      condition <- relevel(condition, ref=ref_group)
      pheno$binary_condition <- factor(ifelse(condition==ref_group, "Control", "Treatment"))
      pheno$analysis_condition <- condition
    }

    # Step 4: Run Differential Expression Analysis using Limma
    design <- model.matrix(~ 0 + pheno$analysis_condition)
    colnames(design) <- levels(pheno$analysis_condition)
    fit <- lmFit(expr_matrix, design)
    
    cont_str <- paste0(levels(pheno$analysis_condition)[2], " - ", levels(pheno$analysis_condition)[1])
    fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts=cont_str, levels=design)))
    
    res <- topTable(fit2, number=Inf)
    res$gene_symbol <- rownames(res)
    
    # Step 5: DEG Identification
    deg <- res %>% filter(adj.P.Val < 0.05 & abs(logFC) >= 1)
    if(nrow(deg) < 5) deg <- res %>% filter(P.Value < 0.05 & abs(logFC) >= 1)
    
    # Step 6: Up & Downregulated genes splitting
    upregulated <- deg %>% filter(logFC > 0) %>% arrange(desc(logFC))
    downregulated <- deg %>% filter(logFC < 0) %>% arrange(logFC)
    
    # Step 7: Save CSV Results
    write.csv(res, file.path(out_dir, "All_Results.csv"), row.names=FALSE)
    write.csv(deg, file.path(out_dir, "DEGs_All.csv"), row.names=FALSE)
    write.csv(upregulated, file.path(out_dir, "Upregulated_genes.csv"), row.names=FALSE)
    write.csv(downregulated, file.path(out_dir, "Downregulated_genes.csv"), row.names=FALSE)
    
    # Step 8: Volcano Plot
    res$category <- "Not Significant"
    res$category[res$gene_symbol %in% upregulated$gene_symbol] <- "Upregulated"
    res$category[res$gene_symbol %in% downregulated$gene_symbol] <- "Downregulated"
    
    p <- ggplot(res, aes(x=logFC, y=-log10(P.Value), color=category)) + 
      geom_point(alpha=0.6) + 
      scale_color_manual(values=c("Downregulated"="blue","Not Significant"="grey","Upregulated"="red")) +
      theme_minimal() + 
      labs(title=paste(geo_id, "Volcano Plot"))
    ggsave(file.path(out_dir, "Volcano_plot.png"), p, width=8, height=6)
    
    # Step 9: Heatmap
    if(nrow(deg) >= 2) {
      top_g <- head(deg$gene_symbol, 40)
      png(file.path(out_dir, "Heatmap.png"), width=800, height=800)
      try(pheatmap(expr_matrix[top_g,], scale="row", show_colnames=FALSE, main=paste("Top DEGs -", geo_id)))
      dev.off()
    }
    
    # Step 10: Calculations of AUC values & ROC Curve
    if(nrow(deg) > 0 && length(unique(pheno$binary_condition)) == 2) {
      group_numeric <- ifelse(pheno$binary_condition == "Treatment", 1, 0)
      top_genes <- head(deg$gene_symbol, 50)
      auc_results <- data.frame(Gene = character(), AUC = numeric(), stringsAsFactors=FALSE)
      
      for (gene in top_genes) {
        if(gene %in% rownames(expr_matrix)) {
          expr <- as.numeric(expr_matrix[gene, ])
          roc_obj <- tryCatch(roc(group_numeric, expr, quiet=TRUE, direction="<"), error=function(e) NULL)
          if(!is.null(roc_obj)) {
            auc_val <- as.numeric(auc(roc_obj))
            if(auc_val < 0.5) auc_val <- as.numeric(auc(roc(group_numeric, expr, quiet=TRUE, direction=">")))
            auc_results <- rbind(auc_results, data.frame(Gene = gene, AUC = auc_val))
          }
        }
      }
      
      if(nrow(auc_results) > 0) {
        auc_results <- auc_results %>% arrange(desc(AUC))
        write.csv(auc_results, file.path(out_dir, "Top50_AUC.csv"), row.names=FALSE)
        
        best_gene <- auc_results$Gene[1]
        best_expr <- as.numeric(expr_matrix[best_gene, ])
        best_roc <- roc(group_numeric, best_expr, quiet=TRUE)
        
        png(file.path(out_dir, paste0("ROC_Curve_", best_gene, ".png")), width=800, height=800)
        plot(best_roc, main=paste("ROC Curve -", best_gene), print.auc=TRUE, col="red", lwd=3)
        dev.off()
      }
    }
    
    cat("SUCCESS\n")
    
  }, error = function(e) {
    cat("FAILED:", e$message, "\n")
  })
}
