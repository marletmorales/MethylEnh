## TF motif enrichment in Dnmt KO differentially-accessible regions
## Figure 7A,B and Supplemental Figure 10B
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TF_motif_enrichment/TFmotifView/Jaspar2026/Set-1_DA_TKO_MATAC5clusters/Results/")

## -- Part 1: Import data ----
TKO_motifs <- read.table(file = "DA_TKO_MATAC5clusters_allDays_wRNA-Cat-Call.txt", header = TRUE, sep = '\t')

## Factorize
TKO_motifs$Gen <- factor(TKO_motifs$Gen, levels = c("TKO", "DKO"))
TKO_motifs$Day <- factor(TKO_motifs$Day, levels = c("D0", "D2", "D4"))
TKO_motifs$direction <- factor(TKO_motifs$direction, levels = c("up", "down"))
TKO_motifs$Sig_subclass <- factor(TKO_motifs$Sig_subclass, levels = c("Sig_both", "Sig_up", "Sig_down", "Sig_none"))
TKO_motifs$Category <- factor(TKO_motifs$Category, levels = c("Primary", "Secondary", "Tertiary", "Quaternary", "NoCG"))
TKO_motifs$Call <- factor(TKO_motifs$Call, levels = c("MethylPlus", "Dual", "MethylMinus", "Little effect", "inconclusive", "NotAssessed"))

## Print stats
cat("Number of motifs:\t",length(unique(TKO_motifs$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(TKO_motifs$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(TKO_motifs$cluster_number)),"\n",sep = "")

## Remove motifs that are "Sig_both" in all three days (e.g. ATF2)
Motif_SigBoth <- TKO_motifs |>
  group_by(motif.ID) |>
  summarise("Num_SigBoth" = sum(Sig_subclass == "Sig_both")) |>
  filter(Num_SigBoth == 6) |>
  distinct() |>
  pull(motif.ID)
cat("Number of motifs significant in all six targets:\t",length(Motif_SigBoth),"\n",sep = "")
TKO_motifs <- TKO_motifs |> filter(!motif.ID %in% Motif_SigBoth)

cat("Remaining motifs:\t",length(unique(TKO_motifs$motif.ID)),"\n",sep = "")
cat("Remaining TFs:\t",length(unique(TKO_motifs$TF)),"\n",sep = "")
cat("Remaining JASPAR clusters:\t",length(unique(TKO_motifs$cluster_number)),"\n",sep = "")
rm(Motif_SigBoth)


## -- Part 2: Create common variables for plotting ----
Annotation <- TKO_motifs |> select(TF, motif.ID, cluster_number, JASPAR_name, Day, Sig_subclass, Category, Call) |> distinct()
Annotation <- Annotation |> pivot_wider(names_from = Day, values_from = Sig_subclass) |> distinct()
Annotation <- Annotation |> tibble::column_to_rownames("motif.ID")
Annotation <- Annotation |> mutate(isCG = case_when(
  Category == "NoCG" ~ "NoCG",
  TRUE ~ "CG"))
Annot_colors <- list(Call = c("MethylPlus" = "#CB2260", "Dual" = "orange2", "MethylMinus" = "#305FD2", "Little effect" = "slategray2", "inconclusive" = "gray60", "NotAssessed" = "white"),
                     Category = c("Primary" = "#D43F3AFF", "Secondary" = "#E18727FF", "Tertiary" = "#EEA236FF", "Quaternary" = "#46B8DAFF", "NoCG" = "gray30"),
                     D0 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     D2 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     D4 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     isCG = c("NoCG" = "cornsilk", "CG" = "chartreuse3"))


## -- Part 3: Plot all TFs (including NoG) ----
## Supplemental Figure 10B

## Keep only one motif.ID per cluster_number (the motif with the most variability in enrichment)
Main_motif <- TKO_motifs |> group_by(motif.ID, cluster_number) |> summarise("Enrichment_sd" = sd(enrichment_peaks_vs_control))
Main_motif <- Main_motif |> group_by(cluster_number) |> slice_max(Enrichment_sd) |> slice_head(n = 1)
Main_motif <- TKO_motifs |> filter(motif.ID %in% Main_motif$motif.ID)

## Make matrix for plotting
Main_matrix <- Main_motif |> select(motif.ID, enrichment_peaks_vs_control, Day, direction)
Main_matrix <- Main_matrix |> pivot_wider(names_from = c(Day, direction), values_from = enrichment_peaks_vs_control, names_sep = ".")
## Reorder columns
Main_matrix <- Main_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
Main_matrix <- as.matrix(Main_matrix |> tibble::column_to_rownames("motif.ID"))

## Plot
Pheat_Main <- pheatmap(
  log2(Main_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = Annotation |> select(Call, isCG, D4, D2, D0),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_row = 30,
  breaks = seq(-0.7, 1.2, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
# ggsave(filename = "Pheatmap_Set-1_TKO_MATAC5clusters_wNoCG_OneMotifPerFam.pdf", plot = Pheat_Main, width = 12, height = 18, units = "cm")
dev.off()
rm(Main_matrix, Main_motif, Pheat_Main)


## -- Part 4: Focus in CG with annotations ----
CG_motifs <- TKO_motifs |> filter(Call != "NotAssessed")

## Print stats
cat("Number of motifs:\t",length(unique(CG_motifs$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(CG_motifs$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(CG_motifs$cluster_number)),"\n",sep = "")

## Pivot wider for plotting
CG_motifs <- CG_motifs |> select(TF, motif.ID, cluster_number, enrichment_peaks_vs_control, p_value, Gen:Call)
CG_motifs <- CG_motifs |> pivot_wider(names_from = direction,
                                      values_from = c(enrichment_peaks_vs_control, p_value),
                                      names_glue = "{.value}_{direction}")

## -- Part 5: Scatterplots ----
## Figure 7A

## Enrichment up vs down, color by Significance
Scatter_sig <- ggplot(CG_motifs, aes(x = enrichment_peaks_vs_control_down, y = enrichment_peaks_vs_control_up, color = Sig_subclass)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Sig_both" = "gray40",
                                "Sig_up" = "#F68532",
                                "Sig_down" = "#266799",
                                "Sig_none" = "gray70")) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Down peaks enrichment", y = "Up peaks enrichment", color = "") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap(~Day)
Scatter_sig
# ggsave(filename = "Enrichment_Up_vs_Down_D0-D2-D4_color_Significance.pdf", plot = Scatter_sig, width = 20, height = 8, units = "cm")

## Enrichment up vs down, color by Call
Scatter_call <- ggplot(CG_motifs, aes(x = enrichment_peaks_vs_control_down, y = enrichment_peaks_vs_control_up, color = Call)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("MethylPlus" = "#CB2260",
                                "Dual" = "orange2",
                                "MethylMinus" = "#305FD2",
                                "Little effect" = "slategray2",
                                "inconclusive" = "gray60")) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Down peaks enrichment", y = "Up peaks enrichment", color = "") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap(~Day)
Scatter_call
# ggsave(filename = "Enrichment_Up_vs_Down_D0-D2-D4_color_Call.pdf", plot = Scatter_call, width = 20, height = 8, units = "cm")


## -- Part 6: Heatmap ----
## Figure 7B

CG_matrix <- CG_motifs |> select(motif.ID, Day, enrichment_peaks_vs_control_up, enrichment_peaks_vs_control_down)
colnames(CG_matrix) <- c("motif.ID", "Day", "up","down")
CG_matrix <- CG_matrix |> pivot_wider(names_from = Day,
                                      values_from = c(up, down),
                                      names_glue = "{Day}.{.value}")
## Reorder columns
CG_matrix <- CG_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
CG_matrix <- as.matrix(CG_matrix |> tibble::column_to_rownames("motif.ID"))

## Plot
Pheat_CG <- pheatmap(
  log2(CG_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = 5,
  annotation_row = Annotation |> select(Call, D4, D2, D0),
  border_color = NA,
  drop_levels = FALSE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_row = 30,
  #breaks = seq(-0.7, 1.2, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
# ggsave(filename = "Pheatmap_Set-1_TKO_MATAC5clusters_OnlywCall_All.pdf", plot = Pheat_CG, width = 12, height = 18, units = "cm")
dev.off()


## -- Part 7: Extract legend ----
Legend <- pheatmap(
  log2(CG_matrix[1:10,]),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(20)),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row = Annotation |> select(D4,Call),
  border_color = NA, drop_levels = TRUE, show_rownames = FALSE,
  #breaks = seq(-0.7, 1.2, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = TRUE,
  legend = FALSE, silent = TRUE)

## Extract grob
index_annot <- which(Legend$gtable$layout$name == "annotation_legend")
Legend <- Legend$gtable$grobs[index_annot]
Legend <- grid.arrange(grobs = Legend)

## Export
#ggsave(filename = "Pheatmap_Set-1_TKO_MATAC5clusters_Legend.pdf", plot = Legend, width = 3, height = 7, units = "cm")
dev.off()

