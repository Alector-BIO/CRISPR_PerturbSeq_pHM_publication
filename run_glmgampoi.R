library('glmGamPoi')
library('Seurat')
library('stringr')
library('dplyr')
library('tidyr')
library('ggplot2')
library('pheatmap')
library('magrittr')
library('SingleCellExperiment')
library('ggrepel')

#Define function
find_log2FC_pvalues <- function(umi_cutoff, seurat_object){
  overall_stats = data.frame()
  sce = as.SingleCellExperiment(seurat_object)
  
  for(gene in genes){
    sce_subset <- sce[rowSums(counts(sce)) > 100, sce$Gene_Targeted %in% c(gene, "NT")]
    assays(sce_subset)[['logcounts']] <- NULL
    
    fit <- glm_gp(sce_subset, design = ~ Gene_Targeted + donor - 1, reference_level="NT", on_disk=FALSE)
    genename = paste0('`Gene_Targeted',gene,"`")
    res <- test_de(fit, eval(parse(text = genename)) - `Gene_TargetedNT`, sort_by = pval)
    
    gene_of_interest = res[res$name == gene,]
    gene_stats = c(umi_cutoff, gene, gene_of_interest$lfc, gene_of_interest$pval, gene_of_interest$adj_pval)
    overall_stats = rbind(overall_stats, gene_stats)
  }
  
  colnames(overall_stats) = c('UMI_cutoff', 'gene', 'log2FC', 'pvalue', 'adj_pval')
  
  overall_stats$pvalue = as.numeric(overall_stats$pvalue)
  overall_stats$adj_pval = as.numeric(overall_stats$adj_pval)
  overall_stats$log2FC = as.numeric(overall_stats$log2FC)
  
  return(overall_stats)
}

#Load RDS files
<<<<<<< Updated upstream
CRISPRi_seurat3 = readRDS("CRISPRi_seurat3.rds")
CRISPRi_seurat5 = readRDS("CRISPRi_seurat5.rds")
CRISPRi_seurat10 = readRDS("CRISPRi_seurat10.rds")
CRISPRi_seurat15 = readRDS("CRISPRi_seurat15.rds")
CRISPRi_seurat20 = readRDS("CRISPRi_seurat20.rds")
CRISPRi_seurat30 = readRDS("CRISPRi_seurat30.rds")
CRISPRi_seurat40 = readRDS("CRISPRi_seurat40.rds")
=======
#CRISPRa_seurat3 = readRDS("CRISPRa_seurat3.rds")
#CRISPRa_seurat5 = readRDS("CRISPRa_seurat5.rds")
#CRISPRa_seurat10 = readRDS("CRISPRa_seurat10.rds")
#CRISPRa_seurat15 = readRDS("CRISPRa_seurat15.rds")
#CRISPRa_seurat20 = readRDS("CRISPRa_seurat20.rds")
#CRISPRa_seurat30 = readRDS("CRISPRa_seurat30.rds")
#CRISPRa_seurat40 = readRDS("CRISPRa_seurat40.rds")
>>>>>>> Stashed changes

#get gene list
genesI <- table(CRISPRi_seurat20[["Gene_Targeted"]]$Gene_Targeted)
genesI <- as.list(genesI)
genesI[setdiff(names(genesI), rownames(CRISPRi_seurat20))] <- NULL
genes <- names(genesI)
#genes <- c('AXL', 'HLA-DRB5', 'ABCA1')

#run functions and save as csv output
<<<<<<< Updated upstream
overall_stats_3umi  <- find_log2FC_pvalues(3, CRISPRi_seurat3)
write.csv(overall_stats_3umi, "overall_stats_3umi_CRISPRi.csv", row.names=FALSE)

overall_stats_5umi  <- find_log2FC_pvalues(5, CRISPRi_seurat5)
write.csv(overall_stats_5umi, "overall_stats_5umi_CRISPRi.csv", row.names=FALSE)

overall_stats_10umi  <- find_log2FC_pvalues(10, CRISPRi_seurat10)
write.csv(overall_stats_10umi, "overall_stats_10umi_CRISPRi.csv", row.names=FALSE)

overall_stats_15umi  <- find_log2FC_pvalues(15, CRISPRi_seurat15)
write.csv(overall_stats_15umi, "overall_stats_15umi_CRISPRi.csv", row.names=FALSE)

overall_stats_20umi  <- find_log2FC_pvalues(20, CRISPRi_seurat20)
write.csv(overall_stats_20umi, "overall_stats_20umi_CRISPRi.csv", row.names=FALSE)

overall_stats_30umi  <- find_log2FC_pvalues(30, CRISPRi_seurat30)
write.csv(overall_stats_30umi, "overall_stats_30umi_CRISPRi.csv", row.names=FALSE)

overall_stats_40umi  <- find_log2FC_pvalues(40, CRISPRi_seurat40)
write.csv(overall_stats_40umi, "overall_stats_40umi_CRISPRi.csv", row.names=FALSE)
=======
#overall_stats_3umi  <- find_log2FC_pvalues(3, CRISPRa_seurat3)
#write.csv(overall_stats_3umi, "overall_stats_3umi_CRISPRa.csv", row.names=FALSE)

#overall_stats_5umi  <- find_log2FC_pvalues(5, CRISPRa_seurat5)
#write.csv(overall_stats_5umi, "overall_stats_5umi_CRISPRa.csv", row.names=FALSE)

#overall_stats_10umi  <- find_log2FC_pvalues(10, CRISPRa_seurat10)
#write.csv(overall_stats_10umi, "overall_stats_10umi_CRISPRa.csv", row.names=FALSE)

#overall_stats_15umi  <- find_log2FC_pvalues(15, CRISPRa_seurat15)
#write.csv(overall_stats_15umi, "overall_stats_15umi_CRISPRa.csv", row.names=FALSE)

#overall_stats_20umi  <- find_log2FC_pvalues(20, CRISPRa_seurat20)
#write.csv(overall_stats_20umi, "overall_stats_20umi_CRISPRa.csv", row.names=FALSE)
>>>>>>> Stashed changes

#overall_stats_30umi  <- find_log2FC_pvalues(30, CRISPRa_seurat30)
#write.csv(overall_stats_30umi, "overall_stats_30umi_CRISPRa.csv", row.names=FALSE)

#overall_stats_40umi  <- find_log2FC_pvalues(40, CRISPRa_seurat40)
#write.csv(overall_stats_40umi, "overall_stats_40umi_CRISPRa.csv", row.names=FALSE)


