## 2026.04.27
## Analyze TF enrichment in open regions with high and low 5hmC
## Adapted from "2026_04_23_TFmotifView_Set-3_5hmC_rich_v2.Rmd"
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TF_motif_enrichment/TFmotifView/Jaspar2026/Set-3_5hmC_rich/Results/")

## -- Part 1: Import 5hmC data ----
hmC_data <- read.table(file = "WT_D4_hMATAC_Top10perc_vs_WT_D4_hMATAC_Bottom10perc_motif_table_-u.txt", header = FALSE, sep = '\t')
colnames(hmC_data) <- c("motif",
                        "number_peaks_with_motif",
                        "number_of_peaks",
                        "number_control_peaks_with_motif",
                        "number_of_control_peaks",
                        "percent_peaks_with_motif",
                        "percent_control_peaks_with_motif",
                        "enrichment_peaks_vs_control",
                        "pval_OR",
                        "pval_UR")
## Add cluster number
clust_num <- read.table(file = "WT_D4_hMATAC_Top10perc_vs_WT_D4_hMATAC_Bottom10perc_motif_table.txt", header = FALSE, sep = '\t')[,1:2]
colnames(clust_num) <- c("motif", "cluster_number")
hmC_data <- left_join(hmC_data, clust_num) |> relocate(cluster_number, .after = "motif")
rm(clust_num)

## Split TF and motif.ID
tmp <- as.data.frame(str_split_fixed(string = hmC_data$motif, pattern = "\\.", n = 2))
colnames(tmp) <- c("TF", "motif.ID")
hmC_data <- cbind(tmp, hmC_data[,colnames(hmC_data) != "motif"])
rm(tmp)

## 1,019 starting motifs, 77 not found in control peaks
## Remove absent motifs (no enrichment calculated)
hmC_data <- hmC_data |> filter(number_control_peaks_with_motif != 0)


## -- Part 2: Import the consensus sequences ----
Consensus <- read.table(file = "../../Anais_source_files/MotifIDs_JASPAR2026_annotated.txt", header = TRUE, sep = "\t", na.strings = "")
colnames(Consensus) <- c("motif.ID", "name", "degenerate_consensus","error")
## Remove rows where motif was not retrieved
Consensus <- Consensus |> filter(!is.na(name))

## Merge Consensus with original table
hmC_data <- left_join(x = hmC_data, y = Consensus[,c("motif.ID", "name", "degenerate_consensus")], by = "motif.ID")

## Convert TF name to all caps to agree with Taipale
hmC_data$TF <- str_to_upper(hmC_data$TF)
hmC_data$name <- str_to_upper(hmC_data$name)

rm(Consensus)


## -- Part 3: Keep only expressed TFs ----
## Import master table from Anais_source_files
TF_expr <- read.table(file = "../../Anais_source_files/MotifIDs_JASPAR2026_wRNAcounts.txt", header = TRUE, sep = '\t')
## In my particular case (Set-3) we are only interested in D4, so I can drop the rest of the columns for simplicity
TF_expr <- TF_expr |> select(motif.ID, JASPAR_name, gene_symbol, WT_D4)
## Rename columns
colnames(TF_expr) <- c("motif.ID", "name", "gene_symbol", "WT_D4")

## Merge tables
hmC_data <- left_join(x = hmC_data, y = TF_expr, by = c("motif.ID", "name"))
rm(TF_expr)

## Remove those motifs with no quantification (-1)
hmC_data <- hmC_data |> filter(WT_D4 != -1)


## -- Part 4: Assess significance ----

## Mark significant motifs
pval_cutoff <- (1/4)^6
hmC_data <- hmC_data |> mutate(
  Significant_OR = case_when(
    pval_OR < pval_cutoff ~ "Yes",
    .default = "No"),
  Significant_UR = case_when(
    pval_UR < pval_cutoff ~ "Yes",
    .default = "No"),
  Significance = factor(case_when(
    Significant_OR == "Yes" & Significant_UR == "No"  ~ "Over",
    Significant_OR == "No"  & Significant_UR == "Yes" ~ "Under",
    Significant_OR == "Yes" & Significant_UR == "Yes" ~ "Both",
    Significant_OR == "No"  & Significant_UR == "No"  ~ "None"),
    levels = c("Over", "Under", "Both", "None")))
cat("Number of over-represented motifs:\t",  sum(hmC_data$Significant_OR == "Yes"), "\n", sep = "")
cat("Number of under-represented motifs:\t", sum(hmC_data$Significant_UR == "Yes"), "\n", sep = "")
cat("Number of non-significant motifs:\t", sum(hmC_data$Significance == "None"), "\n", sep = "")

## Keep only significant motifs
hmC_data <- hmC_data |> filter(Significance != "None")


## -- Part 5: Identify CG presence ----
## -- Primary
## Motifs that contain *specifically* a CpG dinucleotide
CG_primary <- hmC_data[str_detect(string = hmC_data$degenerate_consensus, pattern = "CG"), c("TF", "motif.ID")] |> distinct()

## -- Secondary
## Motifs that _may_ contain a CpG (non-fixed nucleotide)
## R (A|G), S (C|G), K (G|T), V (!T), D (!C), B (!A), N
CG_secondary <- hmC_data[str_detect(string = hmC_data$degenerate_consensus, pattern = "C[RSKVDBN]"), c("TF", "motif.ID")] |> distinct()
## Remove from secondary those motifs that were already considered in primary (CG-containing)
CG_secondary <- dplyr::setdiff(x = CG_secondary, y = CG_primary)

## -- Tertiary
## Motifs that _may_ contain a CpG (non-fixed nucleotide, but from the C side)
## Some TFs analyzed by Taipale are characterized as "MethylPlus" (e.g. HOXA13; seed "NYTYGTAAAAN") despite the C not being fixed
## M (A|C), S (C|G), Y (C|T), V (!T), H (!G), B (!A), N
CG_tertiary <- hmC_data[str_detect(string = hmC_data$degenerate_consensus, pattern = "[MSYVHBN]G"), c("TF", "motif.ID")] |> distinct()
## Remove from tertiary those motifs that were already considered in primary or secondary
CG_tertiary <- dplyr::setdiff(x = CG_tertiary, y = rbind(CG_primary, CG_secondary))

## Add Category to hmC_data
CG_primary$Category <- "Primary"
CG_secondary$Category <- "Secondary"
CG_tertiary$Category <- "Tertiary"
hmC_data <- left_join(x = hmC_data, y = rbind(CG_primary, CG_secondary, CG_tertiary))
rm(CG_primary, CG_secondary, CG_tertiary)


## -- Part 6: Annotating 5mC responsiveness ----
## -- Importing methylation calls from Taipale
MethSelex <- read.table(file = "../../../../Taipale_2017/Taipale_meth_call_clean.txt", header = TRUE, sep = '\t')
colnames(MethSelex) <- c("name", "Call")

## Annotate TFmotifView results
hmC_data <- left_join(x = hmC_data, y = MethSelex, by = "name")
rm(MethSelex)
## Fill in missing Calls
hmC_data$Call[is.na(hmC_data$Call)] <- "NotAssessed"

## For motifs without a Category yet, fill with:
##  a) Quaternary, if assessed by Taipale
##  b) NoCG, if not assesed by Taipale and not Primary/Secondary/Tertiary
hmC_data <- hmC_data |> mutate(Category = case_when(
  !is.na(Category) ~ Category,
  is.na(Category) & Call == "NotAssessed" ~ "NoCG",
  is.na(Category) & Call != "NotAssessed" ~ "Quaternary"))

## Factorize Category and Call
hmC_data$Category <- factor(hmC_data$Category, levels = c("Primary", "Secondary", "Tertiary", "Quaternary", "NoCG"))
hmC_data$Call <- factor(hmC_data$Call, levels = c("MethylPlus", "Dual", "MethylMinus", "Little effect", "inconclusive", "NotAssessed"))


## -- Part 7: Plot all TFs, regardless of CG ----
## Sub: Expressed, significant

## -- Prepare matrix
hmC_data_matrix <- hmC_data |> select(motif.ID, enrichment_peaks_vs_control)
colnames(hmC_data_matrix) <- c("motif.ID", "enrich")
hmC_data_matrix <- as.matrix(hmC_data_matrix |> tibble::column_to_rownames("motif.ID"))

## -- Prepare annotation df
Annotation <- hmC_data |> select(TF, motif.ID, cluster_number, name, Call, Category, Significance)
## Rename columns and factors for simplicity
colnames(Annotation) <- c("TF", "motif.ID", "clust_num", "name", "Call", "Category", "Sig")
Annotation$Sig <- forcats::fct_recode(Annotation$Sig,
                                      "High 5hmC" = "Over",
                                      "Low 5hmC" = "Under")
Annotation$Sig <- forcats::fct_drop(Annotation$Sig)
Annotation <- Annotation |> tibble::column_to_rownames("motif.ID")
Annotation <- Annotation |> mutate(isCG = case_when(
  Category == "NoCG" ~ "NoCG",
  TRUE ~ "CG"))

Annot_colors <- list(Call = c("MethylPlus" = "#CB2260",
                              "Dual" = "orange2",
                              "MethylMinus" = "#305FD2",
                              "Little effect" = "slategray2",
                              "inconclusive" = "gray60",
                              "NotAssessed" = "white"),
                     Category = c("Primary" = "#D43F3AFF",
                                  "Secondary" = "#E18727FF",
                                  "Tertiary" = "#EEA236FF",
                                  "Quaternary" = "#46B8DAFF",
                                  "NoCG" = "gray30"),
                     isCG = c("NoCG" = "cornsilk",
                              "CG" = "chartreuse3"),
                     Sig = c("High 5hmC" = "#AF0663",
                             "Low 5hmC" = "gray90"))

## -- Plot
Pheat_wNoCG <- pheatmap(
  hmC_data_matrix,
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows    = TRUE,
  cluster_cols    = FALSE,
  cutree_rows     = 2,
  cellwidth       = 70,
  border_color    = NA,
  drop_levels     = TRUE,
  show_rownames   = FALSE,
  angle_col       = 45,
  treeheight_row  = 40,
  breaks          = seq(0.3, 1.83, length.out = 101),
  annotation_row  = Annotation |> select(Call, isCG, Sig),
  annotation_colors = Annot_colors,
  annotation_legend = TRUE
)
# ggsave(filename = "Supp_Pheatmap_Set-3_5hmC_rich_wNoCG.pdf", plot = Pheat_wNoCG, width = 12, height = 16, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-3_5hmC_rich/")
dev.off()


## -- Part 8: Plot CG only ----
## Sub: Expressed, significant

## -- Filter for CG content, prepare matrix
CG_data <- hmC_data |> filter(Category != "NoCG")
CG_data_matrix <- CG_data |> select(motif.ID, enrichment_peaks_vs_control)
colnames(CG_data_matrix) <- c("motif.ID", "enrich")
CG_data_matrix <- as.matrix(CG_data_matrix |> tibble::column_to_rownames("motif.ID"))

## -- Plot
Pheat_CG <- pheatmap(
  CG_data_matrix,
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows    = TRUE,
  cluster_cols    = FALSE,
  cutree_rows     = 2,
  cellwidth       = 70,
  border_color    = NA,
  drop_levels     = TRUE,
  show_rownames   = FALSE,
  angle_col       = 45,
  treeheight_row  = 40,
  breaks          = seq(0.3, 1.83, length.out = 101),
  annotation_row  = Annotation |> select(Call, Sig),
  annotation_colors = Annot_colors,
  annotation_legend = TRUE
)
# ggsave(filename = "Supp_Pheatmap_Set-3_5hmC_rich_CGonly.pdf", plot = Pheat_CG, width = 12, height = 16, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-3_5hmC_rich/")
dev.off()


## -- Part 9: Plot annotated CG only ----

Annot_data <- CG_data |> filter(Call != "NotAssessed")
Annot_data_matrix <- Annot_data |> select(motif.ID, enrichment_peaks_vs_control)
colnames(Annot_data_matrix) <- c("motif.ID", "enrich")
Annot_data_matrix <- as.matrix(Annot_data_matrix |> tibble::column_to_rownames("motif.ID"))

## -- Plot
Pheat_AnnotCG <- pheatmap(
  Annot_data_matrix,
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "BrBG", n = 9))(101)),
  cluster_rows    = TRUE,
  cluster_cols    = FALSE,
  cutree_rows     = 2,
  cellwidth       = 70,
  border_color    = NA,
  drop_levels     = TRUE,
  show_rownames   = FALSE,
  angle_col       = 45,
  treeheight_row  = 40,
  breaks          = seq(0.3, 1.83, length.out = 101),
  annotation_row  = Annotation |> select(Call, Sig),
  annotation_colors = Annot_colors,
  annotation_legend = TRUE
)
# ggsave(filename = "Pheatmap_Set-3_5hmC_rich_CG-annot.pdf", plot = Pheat_AnnotCG, width = 12, height = 16, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/TF_motif_enrichment/TFmotifView/Set-3_5hmC_rich/")
dev.off()



