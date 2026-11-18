## 2026.06.30
## Set-1: Parsing TFmotifView results for TKO DA (including TF homology)
## Edit 2026.08.05: Change coloring for isCG_lax for NoCG (Supp 10B)
## Edit 2026.11.12: Change the motif annotation input file (simplified heterodimer expression). Change code for adding annotation to heatmap.
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TF_motif_enrichment/TFmotifView/Jaspar2026/Set-1_DA_TKO_MATAC5clusters/Results/")

## ---------
## >> Part 1: Data import -----------------

## -- Read and transform results
pval_cutoff <- (1/4)^6

Col_names <- c("motif",
               "cluster_number",
               "number_peaks_with_motif",
               "number_of_peaks",
               "number_control_peaks_with_motif",
               "number_of_control_peaks",
               "percent_peaks_with_motif",
               "percent_control_peaks_with_motif",
               "enrichment_peaks_vs_control",
               "p_value")

TKO_D0_up <- read.table(file = "WT_D0_vs_TKO_D0_DEpeaks_MATAC-5clusters_up_vs_Static_original_motif_table.txt", header = FALSE, sep = '\t')
TKO_D0_down <- read.table(file = "WT_D0_vs_TKO_D0_DEpeaks_MATAC-5clusters_down_vs_Static_original_motif_table.txt", header = FALSE, sep = '\t')
TKO_D2_up <- read.table(file = "WT_D2_vs_TKO_D2_DEpeaks_MATAC-5clusters_up_vs_Static_original_motif_table.txt", header = FALSE, sep = '\t')
TKO_D2_down <- read.table(file = "WT_D2_vs_TKO_D2_DEpeaks_MATAC-5clusters_down_vs_Static_original_motif_table.txt", header = FALSE, sep = '\t')
TKO_D4_up <- read.table(file = "WT_D4_vs_TKO_D4_DEpeaks_MATAC-5clusters_up_vs_Static_original_motif_table.txt", header = FALSE, sep = '\t')
TKO_D4_down <- read.table(file = "WT_D4_vs_TKO_D4_DEpeaks_MATAC-5clusters_down_vs_Static_original_motif_table.txt", header = FALSE, sep = '\t')

colnames(TKO_D0_up) <- Col_names
colnames(TKO_D0_down) <- Col_names
colnames(TKO_D2_up) <- Col_names
colnames(TKO_D2_down) <- Col_names
colnames(TKO_D4_up) <- Col_names
colnames(TKO_D4_down) <- Col_names

TKO_D0_up <- TKO_D0_up |> mutate(Gen = "TKO", Day = "D0", direction = "up")
TKO_D0_down <- TKO_D0_down |> mutate(Gen = "TKO", Day = "D0", direction = "down")
TKO_D2_up <- TKO_D2_up |> mutate(Gen = "TKO", Day = "D2", direction = "up")
TKO_D2_down <- TKO_D2_down |> mutate(Gen = "TKO", Day = "D2", direction = "down")
TKO_D4_up <- TKO_D4_up |> mutate(Gen = "TKO", Day = "D4", direction = "up")
TKO_D4_down <- TKO_D4_down |> mutate(Gen = "TKO", Day = "D4", direction = "down")

## Merge
TKO_all <- rbind(TKO_D0_up, TKO_D0_down, TKO_D2_up, TKO_D2_down, TKO_D4_up, TKO_D4_down)
TKO_all$Gen <- factor(TKO_all$Gen, levels = c("TKO", "DKO"))
TKO_all$Day <- factor(TKO_all$Day, levels = c("D0", "D2", "D4"))
TKO_all$direction <- factor(TKO_all$direction, levels = c("up", "down"))
rm(TKO_D0_up, TKO_D0_down, TKO_D2_up, TKO_D2_down, TKO_D4_up, TKO_D4_down, Col_names)

## Separate TF from ID
tmp <- as.data.frame(str_split_fixed(string = TKO_all$motif, pattern = "\\.", n = 2))
colnames(tmp) <- c("TF", "motif.ID")
TKO_all <- cbind(tmp, TKO_all[,colnames(TKO_all) != "motif"])
rm(tmp)

## -- Import motif annotation
Consensus_seqs <- read.table(file = "../../Anais_source_files/MotifIDs_JASPAR2026_annotated_wHomologs_wRNA_wCall_orig_simplified.txt", header = TRUE, sep = "\t", na.strings = "")
## Rename column
Consensus_seqs <- Consensus_seqs |> rename(name_orig = name)

## -- Merge Consensus_seqs with results
TKO_all <- left_join(x = TKO_all, y = Consensus_seqs, by = "motif.ID")
## Reorder columns
TKO_all <- TKO_all |> relocate(c(motif.ID, TF, name_orig))
rm(Consensus_seqs)

cat("Number of motifs:\t",length(unique(TKO_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(TKO_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(TKO_all$cluster_number)),"\n",sep = "")

## -- Export
#write.table(x = TKO_all, file = "TKO_all_table_v3_wHomologs_orig.txt", append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)

## ---------
## >> Part 2: Filtering -----------------

## Keep only motifs that are significant in at least _one_ DA comparison
Summ_sig <- TKO_all |>
  select(motif.ID, p_value, Day, direction) |>
  group_by(motif.ID) |>
  summarise("num_significant" = sum(p_value < pval_cutoff))
table(Summ_sig$num_significant)
Not_sig <- Summ_sig |> filter(num_significant == 0) |> pull(motif.ID)
TKO_all <- TKO_all |> filter(!motif.ID %in% Not_sig)

cat("Number of motifs:\t",length(unique(TKO_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(TKO_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(TKO_all$cluster_number)),"\n",sep = "")
rm(Not_sig, Summ_sig)

## -- How many motifs are significant in each target (Day x direction)?
## Significant motifs
Summ_sig_code <- TKO_all |> select(motif.ID, Day, direction, p_value) |>
  group_by(Day, direction) |>
  summarise("num_significant" = sum(p_value < pval_cutoff))
Summ_sig_code$Perc <- round(Summ_sig_code$num_significant*100/length(unique(TKO_all$motif.ID)), digits = 2)
Summ_sig_code <- Summ_sig_code |>
  pivot_wider(names_from = direction, values_from = c(num_significant, Perc)) |>
  relocate(Day, num_significant_down, Perc_down, num_significant_up, Perc_up)
Summ_sig_code
## D2_down has the least significant, followed by D4_down

## -- Annotate whether motif is significant only in up, only in down, both, or none
Summ_sig_code <- TKO_all |> select(motif.ID, Day, direction, p_value) |>
  group_by(motif.ID, Day, direction)
## Pivot wider
Summ_sig_code <- Summ_sig_code |> pivot_wider(names_from = direction, values_from = p_value)
Summ_sig_code$up <- Summ_sig_code$up < pval_cutoff
Summ_sig_code$down <- Summ_sig_code$down < pval_cutoff
## Make annotation column
Summ_sig_code <- Summ_sig_code |> mutate("Sig_subclass" = case_when(
  up == TRUE & down == TRUE ~ "Sig_both",
  up == TRUE & down == FALSE ~ "Sig_up",
  up == FALSE & down == TRUE ~ "Sig_down",
  up == FALSE & down == FALSE ~ "Sig_none"))
## Merge
TKO_all <- left_join(x = TKO_all, y = Summ_sig_code[,c("motif.ID", "Day", "Sig_subclass")], by = c("motif.ID", "Day"))
## Factorize
TKO_all$Sig_subclass <- factor(TKO_all$Sig_subclass, levels = c("Sig_both", "Sig_up", "Sig_down", "Sig_none"))
## Print stats
tmp <- TKO_all |> select(motif.ID, Day, Sig_subclass) |> distinct()
table(tmp[, c("Day", "Sig_subclass")])
rm(Summ_sig_code, tmp)

## -- Remove motifs that are "Sig_both" in all three days (e.g. ATF2)
Motif_SigBoth <- TKO_all |> select(motif.ID, Day, direction, Sig_subclass) |>
  group_by(motif.ID) |>
  summarise("Num_SigBoth" = sum(Sig_subclass == "Sig_both")) |>
  filter(Num_SigBoth == 6) |>
  pull(motif.ID)
cat("Number of motifs significant in all six targets:\t",length(Motif_SigBoth),"\n",sep = "")
TKO_all <- TKO_all |> filter(!motif.ID %in% Motif_SigBoth)

cat("Number of motifs:\t",length(unique(TKO_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(TKO_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(TKO_all$cluster_number)),"\n",sep = "")
rm(Motif_SigBoth)

## -- Remove non-expressed TFs
## TFs whose expression is < 1 in all days in WT and TKO
Non_expressed <- TKO_all |> select(motif.ID, mouse_homolog, WT_D0:WT_D4, TKO_D0:TKO_D4) |> distinct()
Non_expressed <- Non_expressed |> filter(if_all(WT_D0:TKO_D4, ~.x < 1 | is.na(.x)))
TKO_all <- TKO_all |> filter(!motif.ID %in% Non_expressed$motif.ID)

cat("Number of motifs:\t",length(unique(TKO_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(TKO_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(TKO_all$cluster_number)),"\n",sep = "")
rm(Non_expressed)

## ---------
## >> Part 3: Heatmap annotation -----------------

## -- Make annotation data.frame
Annotation <- TKO_all |> select(motif.ID, Day, Call:Sig_subclass) |> distinct()
## Pivot wider according to the day's significance
Annotation <- Annotation |> pivot_wider(names_from = Day, values_from = Sig_subclass)
## Name rows
Annotation <- Annotation |> tibble::column_to_rownames("motif.ID")
## Factorize columns
Annotation$Call <- factor(Annotation$Call, levels = c("MethylPlus", "MethylMinus", "Dual", "Little effect", "inconclusive","NoData"))
Annotation$CG_category <- factor(Annotation$CG_category, levels = c("Primary", "Secondary", "Tertiary", "Quaternary", "NoCG_sequence"))
Annotation$isCG_strict <- factor(Annotation$isCG_strict, levels = c("CpG", "Non_CpG"))
Annotation$isCG_lax <- factor(Annotation$isCG_lax, levels = c("CpG", "Non_CpG"))
Annotation$isCT <- factor(Annotation$isCT, levels = c("CpT", "Non_CpT"))
Annotation$D0 <- factor(Annotation$D0, levels = c("Sig_both", "Sig_up", "Sig_down", "Sig_none"))
Annotation$D2 <- factor(Annotation$D2, levels = c("Sig_both", "Sig_up", "Sig_down", "Sig_none"))
Annotation$D4 <- factor(Annotation$D4, levels = c("Sig_both", "Sig_up", "Sig_down", "Sig_none"))

## Specify annotation colors
Annot_colors <- list(Call = c("MethylPlus" = "#CB2260", "Dual" = "#c5d107", "MethylMinus" = "#305FD2", "Little effect" = "slategray2", "inconclusive" = "gray60", "NoData" = "white"),
                     CG_category = c("Primary" = "#574571", "Secondary" = "#90719F", "Tertiary" = "#90719F", "Quaternary" = "#DEE0E8", "NoCG_sequence" = "white"),
                     D0 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     D2 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     D4 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     isCG_strict = c("CpG" = "chartreuse3", "Non_CpG" = "cornsilk"),
                     isCG_lax = c("CpG" = "chartreuse3", "Non_CpG" = "cornsilk"))

## ---------
## >> Part 4: Plotting -----------------

## > Plot all TFs (including NoCG) ----
## Supplemental Figure 10B

## -- Keep only one motif.ID per cluster_number (the motif with the most variability in enrichment)
Sub_motifs <- TKO_all |>
  summarise("Enrichment_sd" = sd(enrichment_peaks_vs_control), .by = c(motif.ID, cluster_number)) |>
  slice_max(Enrichment_sd, by = cluster_number) |>
  slice_head(n = 1, by = cluster_number) |>
  pull(motif.ID)

## -- Make matrix for plotting
Sub_matrix <- TKO_all |>
  filter(motif.ID %in% Sub_motifs) |>
  select(motif.ID, enrichment_peaks_vs_control, Day, direction)
Sub_matrix <- Sub_matrix |> pivot_wider(names_from = c(Day, direction), values_from = enrichment_peaks_vs_control, names_sep = ".")
## Reorder columns
Sub_matrix <- Sub_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
Sub_matrix <- as.matrix(Sub_matrix |> tibble::column_to_rownames("motif.ID"))

## -- Define color limit
# plot(x = 1:length(Sub_matrix), y = sort(abs(log2(Sub_matrix))))
# quantile(abs(log2(Sub_matrix)), seq(0.8, 1, 0.01))

#color_lim <- max(abs(log2(Sub_matrix)))
#color_lim <- color_lim * 0.9
color_lim <- 1.1

## -- Plot
k <- 1
Pheat_Sub <- pheatmap(
  log2(Sub_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = k,
  annotation_row = Annotation |> select(Call, isCG_strict, D4, D2, D0),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_row = 25,
  breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
cat("Number of motifs plotted:\t",nrow(Sub_matrix),"\n",sep = "")

# ggsave(filename = "Pheatmap_Set-1_DA_TKO_wNoCG_OneMotifPerFam_wHomologs_v2.pdf", plot = Pheat_Sub, width = 12, height = 18, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")
dev.off()
rm(Sub_matrix, Sub_motifs, Pheat_Sub)

## > Plot TFs with strong CG ----
## Filter for primary CG

CG_all <- TKO_all |> filter(CG_category == "Primary")
## Print stats
cat("Number of motifs:\t",length(unique(CG_all$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(CG_all$TF)),"\n",sep = "")
cat("Number of JASPAR clusters:\t",length(unique(CG_all$cluster_number)),"\n",sep = "")

## -- Make matrix for plotting
CG_matrix <- CG_all |> select(motif.ID, enrichment_peaks_vs_control, Day, direction)
CG_matrix <- CG_matrix |> pivot_wider(names_from = c(Day, direction), values_from = enrichment_peaks_vs_control, names_sep = ".")
## Reorder columns
CG_matrix <- CG_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
CG_matrix <- as.matrix(CG_matrix |> tibble::column_to_rownames("motif.ID"))

## Define color limit
#plot(x = 1:length(CG_matrix), y = sort(abs(log2(CG_matrix))))
#quantile(abs(log2(CG_matrix)), seq(0.8, 1, 0.01))

#color_lim <- max(abs(log2(CG_matrix)))
#color_lim <- color_lim * 0.9
color_lim <- 1.1

## -- Plot
k <- 1
Pheat_CG <- pheatmap(
  log2(CG_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = k,
  annotation_row = Annotation |> select(Call, D4, D2, D0),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_row = 25,
  breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
cat("Number of motifs plotted:\t",nrow(CG_matrix),"\n",sep = "")

# ggsave(filename = "Pheatmap_Set-1_DA_TKO_CGprimary_wHomologs_v2.pdf", plot = Pheat_CG, width = 12, height = 18, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")
dev.off()
rm(Pheat_CG)

## > Plot TFs with strong CG and annotation ----
## Filter for Yin et al. annotation

CG_annot <- CG_all |> select(motif.ID, Call) |> filter(Call != "NoData") |> distinct() |> pull(motif.ID)
CG_matrix_annot <- CG_matrix[CG_annot,]

## -- Plot
k <- 6
Pheat_CG_annot <- pheatmap(
  log2(CG_matrix_annot),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = k,
  annotation_row = Annotation |> select(Call, D4, D2, D0),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_row = 25,
  breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE
)
cat("Number of motifs plotted:\t",nrow(CG_matrix_annot),"\n",sep = "")

# ggsave(filename = "Pheatmap_Set-1_DA_TKO_CGprimaryAnnot_wHomologs_v2.pdf", plot = Pheat_CG, width = 12, height = 18, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")
dev.off()
rm(Pheat_CG_annot)
