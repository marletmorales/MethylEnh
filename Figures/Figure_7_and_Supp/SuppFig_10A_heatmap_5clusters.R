## TF motif enrichment in 5 dynamic clusters
## Suplemental Figure 10A
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
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

## -- Keep only expressed TFs ----
## Import master table
TF_expr <- read.table(file = "../../Anais_source_files/MotifIDs_JASPAR2026_wRNAcounts.txt", header = TRUE, sep = '\t')

## We are only interested in WT at D0-D4, so drop the rest of the columns for simplicity
TF_expr <- TF_expr |> select(motif.ID, JASPAR_name, gene_symbol, WT_D0, WT_D2, WT_D4)
## Rename columns
TF_expr <- TF_expr |> rename(name = JASPAR_name)
## Merge tables
Cluster_all <- left_join(x = Cluster_all, y = TF_expr, by = c("motif.ID", "name"))
rm(TF_expr)

## Remove those motifs in which the count is 0 for all days
Cluster_all <- Cluster_all |> filter(!(WT_D0 <= 0 & WT_D2 <= 0 & WT_D4 <= 0))

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
Annotation <- Main_motif |>
  filter(motif.ID %in% rownames(pheat_matrix)) |>
  select(TF, motif.ID, cluster_number, name, Call, Category) |>
  distinct()
Annotation <- Annotation |> tibble::column_to_rownames("motif.ID")
Annotation <- Annotation |> mutate(isCG = case_when(
  Category == "NoCG" ~ "NoCG",
  TRUE ~ "CG"))
Annot_colors <- list(Call = c("MethylPlus" = "#CB2260",
                              "Dual" = "orange2",
                              "MethylMinus" = "#305FD2",
                              "Little effect" = "slategray2",
                              "inconclusive" = "gray60",
                              "NoData" = "white"),
                     Category = c("Primary" = "#D43F3AFF",
                                  "Secondary" = "#E18727FF",
                                  "Tertiary" = "#EEA236FF",
                                  "Quaternary" = "#46B8DAFF",
                                  "NoCG" = "gray30"),
                     isCG = c("NoCG" = "cornsilk",
                              "CG" = "chartreuse3"))

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
  treeheight_row = 30,
  #breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
# ggsave(filename = "Pheatmap_Set-0_WT_MATAC5clusters_wNoCG_OneMotifPerFam.pdf", plot = Pheat_ClustAll, width = 12, height = 16, units = "cm",
#        path = "Set-0_WT_MATAC5clusters/")
dev.off()
