## 2026.04.17
## Enrichment Analysis for genes associated with gain/loss of H3K27ac in KOs
## Inspired by: https://sib-swiss.github.io/enrichment-analysis-training/
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/GeneOntologies/ROIs/")

## Perform joint analysis for the six ROI sets
Dir <- c("up", "down")
Enh <- c("actv", "actvM", "prim")
Combinations <- expand.grid("Enhancer" = Enh, "Direction" = Dir)

## -----------
## -- Dnmt TKO

## Create empty list of genes
Gene_list_TKO <- list()
## Iterate over every file
for (i in 1:nrow(Combinations)){
  Name_file <- paste0("WT_D4_vs_TKO_D4_K27ac_",Combinations$Enhancer[i],"_",Combinations$Direction[i],"_AssociatedGenes.txt")
  if (!file.exists(Name_file)){
    warning(paste("File not found:", Name_file))
    next
  }
  Gene_names <- read.table(file = Name_file, header = FALSE, sep = '\t')
  Gene_names <- Gene_names$V1
  ## Extract Entrez IDs
  Gene_EntrezID <- mapIds(x = org.Mm.eg.db,
                          keys = Gene_names,
                          column = "ENTREZID",
                          keytype = "ALIAS",
                          multiVals = "first")
  Gene_EntrezID <- Gene_EntrezID[!is.na(Gene_EntrezID)]
  names(Gene_EntrezID) <- NULL

  ## Add to list
  Gene_list_TKO[[paste0(Combinations$Enhancer[i],"_",Combinations$Direction[i])]] <- Gene_EntrezID
}
rm(i, Gene_EntrezID, Gene_names, Name_file)

## -- Run enrichment; add gene symbols
TKO_enrich <- compareCluster(geneCluster = Gene_list_TKO, fun = enrichGO, OrgDb = org.Mm.eg.db, ont="BP")
TKO_enrich <- setReadable(TKO_enrich, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

## Dotplot (showCategory=5 by default)
TKO_dotplot <- dotplot(TKO_enrich, showCategory=4, font.size = 8) +
  set_enrichplot_color(type = "fill", transform='log10', colors=c("#ae017e", "#feebe2")) +
  scale_size(range=c(2, 8))
TKO_dotplot
ggsave(filename = "WT_D4_vs_TKO_D4_DE-K27ac_clusterProfiler_actv-actvM-prim.pdf", plot = TKO_dotplot, device = "pdf", width = 14, height = 14, units = "cm",
       path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/GeneOntologies/clusterProfiler/")

## Colored by genotype
TKO_dotplot_color <- dotplot(TKO_enrich, showCategory=4, font.size = 8) +
  set_enrichplot_color(type = "fill", transform='log10', colors=c("#FF6A00", "#FFF0E6")) +
  scale_size(range=c(2, 8))
TKO_dotplot_color
ggsave(filename = "WT_D4_vs_TKO_D4_DE-K27ac_clusterProfiler_actv-actvM-prim_color.pdf", plot = TKO_dotplot_color, device = "pdf", width = 14, height = 14, units = "cm",
       path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/GeneOntologies/clusterProfiler/")


## --------------------------

## -----------
## -- Tet DKO

## Create empty list of genes
Gene_list_DKO <- list()
## Iterate over every file
for (i in 1:nrow(Combinations)){
  Name_file <- paste0("WT_D4_vs_DKO_D4_K27ac_",Combinations$Enhancer[i],"_",Combinations$Direction[i],"_AssociatedGenes.txt")
  if (!file.exists(Name_file)){
    warning(paste("File not found:", Name_file))
    next
  }
  Gene_names <- read.table(file = Name_file, header = FALSE, sep = '\t')
  Gene_names <- Gene_names$V1
  ## Extract Entrez IDs
  Gene_EntrezID <- mapIds(x = org.Mm.eg.db,
                          keys = Gene_names,
                          column = "ENTREZID",
                          keytype = "ALIAS",
                          multiVals = "first")
  Gene_EntrezID <- Gene_EntrezID[!is.na(Gene_EntrezID)]
  names(Gene_EntrezID) <- NULL

  ## Add to list
  Gene_list_DKO[[paste0(Combinations$Enhancer[i],"_",Combinations$Direction[i])]] <- Gene_EntrezID
}
rm(i, Gene_EntrezID, Gene_names, Name_file)

## -- Run enrichment; add gene symbols
DKO_enrich <- compareCluster(geneCluster = Gene_list_DKO, fun = enrichGO, OrgDb = org.Mm.eg.db, ont="BP")
DKO_enrich <- setReadable(DKO_enrich, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

## Dotplot (showCategory=5 by default)
DKO_dotplot <- dotplot(DKO_enrich, showCategory=3, font.size = 8) +
  set_enrichplot_color(type = "fill", transform='log10', colors=c("#ae017e", "#feebe2")) +
  scale_size(range=c(2, 8))
DKO_dotplot
ggsave(filename = "WT_D4_vs_DKO_D4_DE-K27ac_clusterProfiler_actv-actvM-prim.pdf", plot = DKO_dotplot, device = "pdf", width = 14, height = 14, units = "cm",
       path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/GeneOntologies/clusterProfiler/")

DKO_dotplot_color <- dotplot(DKO_enrich, showCategory=3, font.size = 8) +
  set_enrichplot_color(type = "fill", transform='log10', colors=c("#33334E", "#E1E1E6")) +
  scale_size(range=c(2, 8))
DKO_dotplot_color
ggsave(filename = "WT_D4_vs_DKO_D4_DE-K27ac_clusterProfiler_actv-actvM-prim_color.pdf", plot = DKO_dotplot_color, device = "pdf", width = 14, height = 14, units = "cm",
       path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/GeneOntologies/clusterProfiler/")
