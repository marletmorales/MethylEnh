## 2026.04.15
## Analyze TF enrichment in up vs down regions of TKO DiffAcc
## Edit 2026.04.18: Make supplementary figure with NoCG (adapted code from 2026_04_16_TFmotifView_Set-0_MATAC5clusters_pheatmap_withNoCG.R)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TF_motif_enrichment/TFmotifView/Jaspar2026/Set-1_DA_TKO_MATAC5clusters/Results/")

## -- Part 1: Import data ----
TKO_sub <- read.table(file = "DA_TKO_MATAC5clusters_allDays_wRNA-Cat-Call.txt", header = TRUE, sep = '\t')
## Factorize
TKO_sub$Gen <- factor(TKO_sub$Gen, levels = c("TKO", "DKO"))
TKO_sub$Day <- factor(TKO_sub$Day, levels = c("D0", "D2", "D4"))
TKO_sub$direction <- factor(TKO_sub$direction, levels = c("up", "down"))
TKO_sub$Sig_subclass <- factor(TKO_sub$Sig_subclass, levels = c("Sig_both", "Sig_up", "Sig_down", "Sig_none"))
TKO_sub$Category <- factor(TKO_sub$Category, levels = c("Primary", "Secondary", "Tertiary", "Quaternary", "NoCG"))
TKO_sub$Call <- factor(TKO_sub$Call, levels = c("MethylPlus", "Dual", "MethylMinus", "Little effect", "inconclusive", "NotAssessed"))
## Filter for CG motifs
CG_all <- TKO_sub |> filter(Category != "NoCG")
CG_sub <- CG_all |> filter(Call != "NotAssessed")

## Print stats
tmp <- CG_sub |> select(TF, motif.ID, Category, Call) |> distinct()
cat("Number of motifs:\t",length(unique(tmp$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(tmp$TF)),"\n",sep = "")
cat("Distribution of CG category\n",sep = "")
table(tmp$Category, exclude = "NoCG")
cat("Distribution of MethylSELEX call\n",sep = "")
table(tmp$Call, exclude = "NotAssessed")
cat("Distribution of Significance call\n",sep = "")
table(CG_sub |> filter(direction == "up") |> select(Day, Sig_subclass))
rm(tmp)

## -- Remove motifs that are "Sig_both" in all three days (e.g. ATF2)
Motif_OnlySigBoth <- CG_sub |>
  group_by(motif.ID) |>
  summarise("Num_SigBoth" = sum(Sig_subclass == "Sig_both")) |>
  filter(Num_SigBoth == 6) |>
  distinct() |>
  pull(motif.ID)
cat("Number of motifs significant in all six targets:\t",length(Motif_OnlySigBoth),"\n",sep = "")
CG_sub <- CG_sub |> filter(!motif.ID %in% Motif_OnlySigBoth)

## Print stats
tmp <- CG_sub |> select(TF, motif.ID, Category, Call) |> distinct()
cat("Number of motifs:\t",length(unique(tmp$motif.ID)),"\n",sep = "")
cat("Number of TFs:\t",length(unique(tmp$TF)),"\n",sep = "")
cat("Distribution of CG category\n",sep = "")
table(tmp$Category, exclude = "NoCG")
cat("Distribution of MethylSELEX call\n",sep = "")
table(tmp$Call, exclude = "NotAssessed")
cat("Distribution of Significance call\n",sep = "")
table(CG_sub |> filter(direction == "up") |> select(Day, Sig_subclass))
rm(tmp, Motif_OnlySigBoth)

## -- Part 2: Transform for plotting ----

## Pivot wider, only for relevant columns
CG_sub_wider <- CG_sub |> select(TF, motif.ID, cluster_number, enrichment_peaks_vs_control, p_value, Gen:Call)
CG_sub_wider <- CG_sub_wider |> pivot_wider(names_from = direction,
                                            values_from = c(enrichment_peaks_vs_control, p_value),
                                            names_glue = "{.value}_{direction}")

## -- Stats
## Relationship between Significance and Call (single normalization)
Summ_Sig_Call <- CG_sub_wider |> group_by(Day, Sig_subclass, Call) |> summarise("N" = n())
## Total motifs per class
Tot_class <- CG_sub_wider |> select(motif.ID, Call) |> distinct() |> group_by(Call) |> summarise("Total" = n())
## Calculate percentage
Summ_Sig_Call <- left_join(x = Summ_Sig_Call, y = Tot_class)
Summ_Sig_Call <- Summ_Sig_Call |> mutate(Perc = round(N*100/Total, digits = 2))
## Pivot wider
Summ_Sig_Call <- Summ_Sig_Call |> pivot_wider(names_from = Day, values_from = c(N, Perc), values_fill = 0)
## Clean table
Summ_Sig_Call <- Summ_Sig_Call |> arrange(Sig_subclass, Call)
Summ_Sig_Call <- Summ_Sig_Call |> relocate(Total, .after = last_col())
View(Summ_Sig_Call |> arrange(Call, Sig_subclass))
rm(Summ_Sig_Call)

## Relationship between Significance and Call (double normalization)
Summ_Sig_Call_v2 <- CG_sub_wider |> group_by(Day, Sig_subclass, Call) |> summarise("N" = n())
## Total motifs per class
Tot_class <- CG_sub_wider |> select(motif.ID, Call) |> distinct() |> group_by(Call) |> summarise("Total_call" = n())
Tot_sig <- CG_sub_wider |> select(motif.ID, Day, Sig_subclass) |> distinct() |> group_by(Day, Sig_subclass)  |> summarise("Total_sig" = n())
Summ_Sig_Call_v2 <- left_join(x = Summ_Sig_Call_v2, y = Tot_class, by = "Call")
Summ_Sig_Call_v2 <- left_join(x = Summ_Sig_Call_v2, y = Tot_sig, by = c("Day", "Sig_subclass"))
## Calculate percentage
Summ_Sig_Call_v2 <- Summ_Sig_Call_v2 |> mutate(Perc_call = round(N*100/Total_call, digits = 2))
Summ_Sig_Call_v2 <- Summ_Sig_Call_v2 |> mutate(Perc_sig = round(N*100/Total_sig, digits = 2))
## Clean table
Summ_Sig_Call_v2 <- Summ_Sig_Call_v2 |> arrange(Day, Sig_subclass, Call)
Summ_Sig_Call_v2 <- Summ_Sig_Call_v2 |> relocate(Day, Sig_subclass, Call, N, Perc_call, Total_call, Perc_sig, Total_sig)
View(Summ_Sig_Call_v2 |> arrange(Call, Sig_subclass))

rm(Tot_class, Tot_sig)

## -- Part 3: Scatter-Plotting ----

## Enrichment up vs down, color by Significance
Scatter_sig <- ggplot(CG_sub_wider, aes(x = enrichment_peaks_vs_control_down, y = enrichment_peaks_vs_control_up, color = Sig_subclass)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Sig_both" = "gray40",
                                "Sig_up" = "#F68532",
                                "Sig_down" = "#266799",
                                "Sig_none" = "gray70")) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Down: Enrichment vs background", y = "Up: Enrichment vs background", color = "") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap(~Day)
Scatter_sig
# ggsave(filename = "Enrichment_Up_vs_Down_D0-D2-D4_color_Significance.pdf", plot = Scatter_sig, width = 20, height = 8, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")

## Enrichment up vs down, color by Call
Scatter_call <- ggplot(CG_sub_wider, aes(x = enrichment_peaks_vs_control_down, y = enrichment_peaks_vs_control_up, color = Call)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("MethylPlus" = "#CB2260",
                                "Dual" = "orange2",
                                "MethylMinus" = "#305FD2",
                                "Little effect" = "slategray2",
                                "inconclusive" = "gray60")) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Down: Enrichment vs background", y = "Up: Enrichment vs background", color = "") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap(~Day)
Scatter_call
# ggsave(filename = "Enrichment_Up_vs_Down_D0-D2-D4_color_Call.pdf", plot = Scatter_call, width = 20, height = 8, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")

## -- Part 4: Heatmapping ----

CG_sub_df <- CG_sub_wider |> select(JASPAR_name, motif.ID, cluster_number, Day, Sig_subclass, Category, Call, enrichment_peaks_vs_control_up, enrichment_peaks_vs_control_down, p_value_up, p_value_down)
CG_sub_matrix <- CG_sub_df |> select(motif.ID, Day, enrichment_peaks_vs_control_up, enrichment_peaks_vs_control_down)
colnames(CG_sub_matrix) <- c("motif.ID", "Day", "up","down")
CG_sub_matrix <- CG_sub_matrix |> pivot_wider(names_from = Day,
                                              values_from = c(up, down),
                                              names_glue = "{Day}.{.value}")
## Reorder columns
CG_sub_matrix <- CG_sub_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
CG_sub_matrix <- as.matrix(CG_sub_matrix |> tibble::column_to_rownames("motif.ID"))

## Generate annotation df
Annotation <- CG_sub_df |> select(JASPAR_name:Call) |> pivot_wider(names_from = Day, values_from = Sig_subclass) |> distinct()
Annotation <- Annotation |> tibble::column_to_rownames("motif.ID")
Annot_colors <- list(Call = c("MethylPlus" = "#CB2260", "Dual" = "orange2", "MethylMinus" = "#305FD2", "Little effect" = "slategray2", "inconclusive" = "gray60"),
                     Category = c("Primary" = "#D43F3AFF", "Secondary" = "#E18727FF", "Tertiary" = "#EEA236FF", "Quaternary" = "#46B8DAFF"),
                     D0 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     D2 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"),
                     D4 = c("Sig_both" = "gray40", "Sig_up" = "#F68532", "Sig_down" = "#266799", "Sig_none" = "gray70"))

color_lim <- max(abs(log2(CG_sub_matrix)))

## -- Default plot
Pheat_CGsub <- pheatmap(
  log2(CG_sub_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = 5,
  annotation_row = Annotation |> select(Call, D4, D2, D0),
  border_color = NA,
  drop_levels = FALSE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_col = 30,
  treeheight_row = 30,
  #breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE,
)
# ggsave(filename = "Pheatmap_Set-1_TKO_MATAC5clusters_OnlywCall_All.pdf", plot = Pheat_CGsub, width = 12, height = 18, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")
dev.off()

## -- Part 4.5: Extract legend
Legend <- pheatmap(
  log2(CG_sub_matrix[1:10,]),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(20)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  #cutree_rows = 4,
  annotation_row = Annotation |> select(D4,Call),
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
# ggsave(filename = "Pheatmap_Set-1_TKO_MATAC5clusters_Legend.pdf", plot = Legend, width = 8, height = 6, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")
#dev.off()

## -- Part 5: Keep only one motif per cluster ----

cat("Total number of motifs:\t",length(unique(CG_sub_df$motif.ID)),"\n",sep = "")
cat("Total number of clusters:\t",length(unique(CG_sub_df$cluster_number)),"\n",sep = "")

## -- Do all motifs of the same cluster have the same call?
Repeated_clust <- Annotation |> select(cluster_number, Call) |> distinct() |> pull(cluster_number)
Repeated_clust <- Repeated_clust[duplicated(Repeated_clust)]
Repeated_clust <- Annotation |> filter(cluster_number %in% Repeated_clust) |> select(cluster_number, Call) |> distinct()
cat("Number of clusters with more than one Call:\t",length(unique(Repeated_clust$cluster_number)),"\n",sep = "")
View(Repeated_clust |> arrange(cluster_number, Call))

## 264 potentially problematic (IRFs)
## 355 potentially problematic (SNAI1, TCF4)
## All other clusters are either MethylPlus or MethylMinus combined with another lesser category (Dual, Little effect, inconclusive)

## -- Retain one motif per cluster using sequential prioritization rules
## Code adapted from 2026_04_02_TFmotifView_Jaspar2026_parsing_v2.R
## Written by Claude Opus, adapted by Marlet

# Define prioritization functions
call_priority <- function(call) { case_when(
  call %in% c("MethylPlus", "MethylMinus", "Dual") ~ 1,
  call == "Little effect" ~ 2,
  call == "inconclusive" ~ 3,
  TRUE ~ 4)}
sig_priority <- function(sig) { case_when(
  sig %in% c("Sig_up", "Sig_down") ~ 2,
  sig == "Sig_both" ~ 1,
  sig == "Sig_none" ~ 0,
  TRUE ~ 0)}
cat_priority <- function(cat) { case_when(
  cat == "Primary" ~ 2,
  cat %in% c("Secondary", "Tertiary") ~ 1,
  cat == "Quaternary" ~ 0,
  TRUE ~ 0)}

# Calculate scores for each motif.ID within each cluster_number
motif_scores <- CG_sub_df |> group_by(cluster_number, motif.ID) |> summarise(
  # Rule 1: Best (lowest) Call priority
  best_call_priority = min(call_priority(Call)),
  # Rule 2: Count direction of significance
  sig_score = sum(sig_priority(Sig_subclass)),
  # Rule 3: Best (highest) CG Category
  best_cat_priority = max(cat_priority(Category)),
  # Rule 4: Biggest difference between enr_up and enr_down
  enr_dif = max(abs(enrichment_peaks_vs_control_up - enrichment_peaks_vs_control_down)),
  .groups = "drop")

# Apply sequential prioritization rules
selected_motifs <- motif_scores |> group_by(cluster_number) |>
  # Rule 1: Keep only those with the best (lowest) call priority
  filter(best_call_priority == min(best_call_priority)) |>
  # Rule 2: Among remaining, keep those with the best pattern of significance
  filter(sig_score == max(sig_score)) |>
  # Rule 3: Among remaining, keep those with the best CG category
  filter(best_cat_priority == max(best_cat_priority)) |>
  # Rule 4: Among remaining, keep the one with the highest enrichment difference
  filter(enr_dif == max(enr_dif)) |>
  ungroup()

# Check for ties after all filters
ties_check <- selected_motifs |> group_by(cluster_number) |> summarise(n_motifs = n(), .groups = "drop") |> filter(n_motifs > 1)
if (nrow(ties_check) > 0) {
  cat("Clusters with more than 1 observation:\n")
  selected_motifs |> filter(cluster_number %in% ties_check$cluster_number)
}

## -- Re-do heatmap with selected motifs
CG_sub_matrix_select <- CG_sub_matrix[selected_motifs$motif.ID,]
color_lim_select <- max(abs(log2(CG_sub_matrix_select)))

## Default plot
pheatmap(
  log2(CG_sub_matrix_select),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  #cutree_rows = 5,
  annotation_row = Annotation |> select(Call, D4, D2, D0),
  border_color = NA,
  drop_levels = FALSE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_col = 30,
  treeheight_row = 30,
  #breaks = seq(-color_lim, color_lim, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE,
)


## --------------------------------------

## ---------------------------------------------------------------------- ##
## -- Extra: Defining top motifs per target, independent of CG content -- ##
## ---------------------------------------------------------------------- ##

## -- Top 50 by p-value
pval_cutoff <- (1/4)^6
Top_pval <- TKO_sub |> group_by(Day, direction) |> filter(p_value < pval_cutoff) |> slice_min(p_value, n = 50)
D0_up_pval <- sort(unique(Top_pval |> filter(Day == "D0" & direction == "up") |> pull(motif.ID)))
D0_down_pval <- sort(unique(Top_pval |> filter(Day == "D0" & direction == "down") |> pull(motif.ID)))
D2_up_pval <- sort(unique(Top_pval |> filter(Day == "D2" & direction == "up") |> pull(motif.ID)))
D2_down_pval <- sort(unique(Top_pval |> filter(Day == "D2" & direction == "down") |> pull(motif.ID)))
D4_up_pval <- sort(unique(Top_pval |> filter(Day == "D4" & direction == "up") |> pull(motif.ID)))
D4_down_pval <- sort(unique(Top_pval |> filter(Day == "D4" & direction == "down") |> pull(motif.ID)))

## -- Top 50 by enrichment
Top_enr <- TKO_sub |> group_by(Day, direction) |> filter(p_value < pval_cutoff) |> slice_max(enrichment_peaks_vs_control, n = 50)
D0_up_enr <- sort(unique(Top_enr |> filter(Day == "D0" & direction == "up") |> pull(motif.ID)))
D0_down_enr <- sort(unique(Top_enr |> filter(Day == "D0" & direction == "down") |> pull(motif.ID)))
D2_up_enr <- sort(unique(Top_enr |> filter(Day == "D2" & direction == "up") |> pull(motif.ID)))
D2_down_enr <- sort(unique(Top_enr |> filter(Day == "D2" & direction == "down") |> pull(motif.ID)))
D4_up_enr <- sort(unique(Top_enr |> filter(Day == "D4" & direction == "up") |> pull(motif.ID)))
D4_down_enr <- sort(unique(Top_enr |> filter(Day == "D4" & direction == "down") |> pull(motif.ID)))

## -- Assessing intersects ----
## -- Intersect of up vs down; same day
## P-value
length(intersect(D0_up_pval, D0_down_pval))
length(intersect(D2_up_pval, D2_down_pval))
length(intersect(D4_up_pval, D4_down_pval))
## Enrichment
length(intersect(D0_up_enr, D0_down_enr))
length(intersect(D2_up_enr, D2_down_enr))
length(intersect(D4_up_enr, D4_down_enr))

## -- Intersect of up; diff day
## P-value
length(intersect(D0_up_pval, D2_up_pval))
length(intersect(D0_up_pval, D4_up_pval))
length(intersect(D2_up_pval, D4_up_pval))
## Enrichment
length(intersect(D0_up_enr, D2_up_enr))
length(intersect(D0_up_enr, D4_up_enr))
length(intersect(D2_up_enr, D4_up_enr))

## -- Intersect of down; diff day
## P-value
length(intersect(D0_down_pval, D2_down_pval))
length(intersect(D0_down_pval, D4_down_pval))
length(intersect(D2_down_pval, D4_down_pval))
## Enrichment
length(intersect(D0_down_enr, D2_down_enr))
length(intersect(D0_down_enr, D4_down_enr))
length(intersect(D2_down_enr, D4_down_enr))

## Clean workspace
rm(D0_down_enr, D0_down_pval, D0_up_enr, D0_up_pval,
   D2_down_enr, D2_down_pval, D2_up_enr, D2_up_pval,
   D4_down_enr, D4_down_pval, D4_up_enr, D4_up_pval,
   Top_enr, Top_pval)

## ------------------------------------

## ---------------------- ##
## -- Extra: plot NoCG -- ##
## ---------------------- ##

## -- Plot all motif.IDs; remove those that are "Sig_both" in all three days
Motif_OnlySigBoth <- TKO_sub |>
  group_by(motif.ID) |>
  summarise("Num_SigBoth" = sum(Sig_subclass == "Sig_both")) |>
  filter(Num_SigBoth == 6) |>
  distinct() |>
  pull(motif.ID)
cat("Number of motifs significant in all six targets:\t",length(Motif_OnlySigBoth),"\n",sep = "")
TKO_sub <- TKO_sub |> filter(!motif.ID %in% Motif_OnlySigBoth)
cat("Remaining motifs:\t",length(unique(TKO_sub$motif.ID)),"\n",sep = "")

## -- Generate matrix
pheat_matrix <- TKO_sub |> select(motif.ID, enrichment_peaks_vs_control, Day, direction)
pheat_matrix <- pheat_matrix |> pivot_wider(names_from = c(Day, direction), values_from = enrichment_peaks_vs_control, names_sep = ".")
## Reorder columns
pheat_matrix <- pheat_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
pheat_matrix <- as.matrix(pheat_matrix |> tibble::column_to_rownames("motif.ID"))

## Generate annotation df
Annotation <- TKO_sub |>
  filter(motif.ID %in% rownames(pheat_matrix) & direction == "up") |>
  select(TF, motif.ID, cluster_number, JASPAR_name, Day, Sig_subclass, Category, Call) |> distinct()
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

## -- Default plot
pheatmap(
  log2(pheat_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = Annotation |> select(Call, isCG, D4, D2, D0),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_col = 30,
  treeheight_row = 30,
  breaks = seq(-0.7, 1.2, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE,
)

## ---------------------------------
## ---------------------------------
## Keep only one motif.ID per cluster_number (the motif with the most variability in enrichment) (n = 219)
Motif_var <- TKO_sub |> group_by(motif.ID, cluster_number) |> summarise("Enrichment_sd" = sd(enrichment_peaks_vs_control))
Main_motif <- Motif_var |> group_by(cluster_number) |> slice_max(Enrichment_sd) |> slice_head(n = 1)
Main_motif <- TKO_sub |> filter(motif.ID %in% Main_motif$motif.ID)
## Keep only relevant columns
Main_motif <- Main_motif |> select(TF, motif.ID, cluster_number, enrichment_peaks_vs_control, p_value, Day, direction, Call, Category)

## Make matrix for plotting
pheat_matrix <- Main_motif |> select(motif.ID, enrichment_peaks_vs_control, Day, direction)
pheat_matrix <- pheat_matrix |> pivot_wider(names_from = c(Day, direction), values_from = enrichment_peaks_vs_control, names_sep = ".")
## Reorder columns
pheat_matrix <- pheat_matrix[,c("motif.ID", "D0.down", "D2.down", "D4.down", "D0.up", "D2.up", "D4.up")]
pheat_matrix <- as.matrix(pheat_matrix |> tibble::column_to_rownames("motif.ID"))

## -- Default plot
Pheat_TKOsub <- pheatmap(
  log2(pheat_matrix),
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = Annotation |> select(Call, isCG, D4, D2, D0),
  border_color = NA,
  drop_levels = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  treeheight_col = 30,
  treeheight_row = 30,
  breaks = seq(-0.7, 1.2, length.out = 101),
  annotation_colors = Annot_colors,
  annotation_legend = FALSE,
)
# ggsave(filename = "Pheatmap_Set-1_TKO_MATAC5clusters_wNoCG_OneMotifPerFam.pdf", plot = Pheat_TKOsub, width = 12, height = 18, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-1_DA_TKO_MATAC5clusters/")
dev.off()
