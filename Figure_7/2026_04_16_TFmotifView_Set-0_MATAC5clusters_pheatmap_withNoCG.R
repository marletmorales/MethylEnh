## 2026.04.16
## TfmotifView with JASPAR2026; WT MATAC 5 clusters
## See "2026_03_31_TFmotifView_MATAC5clusters_pheatmap.Rmd" for reference
## Plot all motifs, including NoCG
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TF_motif_enrichment/TFmotifView/Jaspar2026/Set-0_WT_MATAC_IndividualClusters_StaticBackground/Results/")

## -- Read TfmotifView data
Cluster_all <- read.table(file = "Clust_all_table.txt", header = TRUE, sep = '\t')
cat("Number of motifs:\t",length(unique(Cluster_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(Cluster_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(Cluster_all$cluster_number)),"\n",sep = "")
## 1,019 starting motifs, 77 not found in control peaks

## Keep only motifs that are significant in at least one cluster
pval_cutoff <- (1/4)^6
tmp <- unique(Cluster_all |> filter(p_value < pval_cutoff) |> pull(motif.ID))
Cluster_all <- Cluster_all[Cluster_all$motif.ID %in% tmp,]
rm(tmp)
cat("Number of motifs:\t",length(unique(Cluster_all$motif.ID)),"\n",sep = "")

## -- Add RNA expression that I obtained from Set-1
## Too long to do from scratch here
TF_RNAcounts <- read.table(file = "../../Set-1_DA_TKO_MATAC5clusters/Results/DA_TKO_MATAC5clusters_allDays_wRNA-Cat-Call.txt", header = TRUE, sep = '\t')
TF_RNAcounts <- TF_RNAcounts |> select(motif.ID, Day, WT) |> distinct()
TF_RNAcounts <- TF_RNAcounts |> pivot_wider(names_from = Day, values_from = WT, names_prefix = "WT.")
## Remove those motifs in which the count is 0 for all days
TF_RNAcounts <- TF_RNAcounts |> filter(!(WT.D0 == 0 & WT.D2 == 0 & WT.D4 == 0))

## Filter Cluster_all motifs
Cluster_all <- Cluster_all |> filter(motif.ID %in% TF_RNAcounts$motif.ID)
cat("Number of motifs:\t",length(unique(Cluster_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(Cluster_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(Cluster_all$cluster_number)),"\n",sep = "")


## -------------------------------- ##
## -- Prepare table for pheatmap -- ##
## -------------------------------- ##

## Keep only one motif.ID per cluster_number (the motif with the most variability in enrichment)
Motif_var <- Cluster_all |> group_by(motif.ID, cluster_number) |> summarise("Enrichment_sd" = sd(enrichment_peaks_vs_control))
Main_motif <- Motif_var |> group_by(cluster_number) |> slice_max(Enrichment_sd) |> slice_head(n = 1)
Main_motif <- Cluster_all |> filter(motif.ID %in% Main_motif$motif.ID)
## Keep only relevant columns
Main_motif <- Main_motif |> select(TF, motif.ID, cluster_number, enrichment_peaks_vs_control, p_value, clust, name, degenerate_consensus, Call, Category)

## Make matrix for plotting
pheat_matrix <- Main_motif |> select(motif.ID, enrichment_peaks_vs_control, clust) |> pivot_wider(names_from = clust, values_from = enrichment_peaks_vs_control)
pheat_matrix <- as.matrix(pheat_matrix |> tibble::column_to_rownames(var = "motif.ID"))

## Generate annotation df
Annotation <- Main_motif |> filter(motif.ID %in% rownames(pheat_matrix)) |> select(TF, motif.ID, cluster_number, name, Call, Category) |> distinct()
Annotation <- Annotation |> tibble::column_to_rownames("motif.ID")
Annotation <- Annotation |> mutate(isCG = case_when(
  Category == "NoCG" ~ "NoCG",
  TRUE ~ "CG"))
Annot_colors <- list(Call = c("MethylPlus" = "#CB2260", "Dual" = "orange2", "MethylMinus" = "#305FD2", "Little effect" = "slategray2", "inconclusive" = "gray60", "NoData" = "white"),
                     Category = c("Primary" = "#D43F3AFF", "Secondary" = "#E18727FF", "Tertiary" = "#EEA236FF", "Quaternary" = "#46B8DAFF", "NoCG" = "gray30"),
                     isCG = c("NoCG" = "cornsilk", "CG" = "chartreuse3"))

color_lim <- max(abs(log2(pheat_matrix)))

Pheat_ClustAll <- pheatmap(
  log2(pheat_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  #cutree_rows = 4,
  annotation_row = Annotation |> select(Call, isCG),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 0,
  treeheight_col = 30,
  treeheight_row = 30,
  #breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
Pheat_ClustAll
# ggsave(filename = "Pheatmap_Set-0_WT_MATAC5clusters_wNoCG_OneMotifPerFam.pdf", plot = Pheat_ClustAll, width = 12, height = 16, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-0_WT_MATAC5clusters/")
dev.off()


## -- Export legend
Legend <- pheatmap(
  log2(pheat_matrix[1:10,]),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(20)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  #cutree_rows = 4,
  annotation_row = Annotation |> select(Call, isCG),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 0,
  treeheight_col = 30,
  treeheight_row = 30,
  #breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = TRUE
)
# ggsave(filename = "Pheatmap_Set-0_WT_MATAC5clusters_Legend.pdf", plot = Legend, width = 8, height = 6, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-0_WT_MATAC5clusters/")
dev.off()
