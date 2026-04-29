## 2026.04.08
## Link gene expression to levels of 5mC and 5hmC
## Based on 2026_04_08_5mC-5hmC_vs_RNAseq.R; adapted for paper figure
library(dplyr)
library(ggplot2)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/MEP_hMATAC_RNAseq/")

## -- Read data
meC_data <- read.table(file = "../Measure_DNAme/DNAme_ROIs/Enhancers_intersect_v3/meanDNAme_WT_D0-D2-D4_K4me1-MATAC-K27ac-WGBS__Coordinates_WT_Master_Enhancer_Peakome_v3.txt", header = TRUE, sep = '\t')
hmC_data <- read.table(file = "../Measure_DNAme/DNAme_ROIs/WT_Master_Enhancer_Peakome_v6/WT_D0-D2-D4_hMATAC__Coordinates_WT_Master_Enhancer_Peakome_v6.txt", header = TRUE, sep = '\t')
RNA_data <- read.table(file = "../Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_annotated_TargetGene_wExpression_wActvM.txt", header = TRUE, sep = '\t')

## -- Clean data
## 5mC
meC_data <- meC_data |> select(chr:end, D0_class:D4_class, matches("num.cpgs_WT_D[024]_MATAC"), matches("meth.wmean_WT_D[024]_MATAC")) |> arrange(chr, start, end)
## 5hmC
hmC_data <- hmC_data |> select(!coord)
colnames(hmC_data) <- sub(pattern = "num.cpgs_D", replacement = "num.cpgs_WT_D",x = colnames(hmC_data))
colnames(hmC_data) <- sub(pattern = "meth.wmean_D", replacement = "meth.wmean_WT_D",x = colnames(hmC_data))
## dplyr version:
## hmC_data <- hmC_data |> rename_with(~ sub("num.cpgs_D", "num.cpgs_WT_D", .x))

## -- Merge data
## Keep in mind that meC_data has additional rows as it comes from Enhancers_intersect_v3
RNA_data <- left_join(x = RNA_data, y = meC_data)
RNA_data <- left_join(x = RNA_data, y = hmC_data)
RNA_data <- RNA_data |> relocate(contains("num.cpgs_"), .after = last_col())
rm(meC_data, hmC_data)

## -- Pivot to longer for plotting
RNA_data_longer <- RNA_data |> tidyr::pivot_longer(cols = matches("^(meth\\.wmean|num\\.cpgs)_WT_"),
                                                   names_to = c(".value", "Day", "Mark"),
                                                   names_pattern = "^(meth\\.wmean|num\\.cpgs)_WT_(.+)_(.+)$")
## Drop missing 5mC/5hmC
RNA_data_longer <- RNA_data_longer |> filter(!is.na(meth.wmean))
## Factorize mark
RNA_data_longer$Mark <- factor(RNA_data_longer$Mark, levels = c("hMATAC", "MATAC"))

## ------------------- ##
## ---- Plot time ---- ##
## ------------------- ##

Output_folder <- "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_6/Supp_5mC_5hmC_vs_RNA"

## -------
## -- D0

RNA_D0_longer <- RNA_data_longer[RNA_data_longer$Day == "D0" & RNA_data_longer$D0_class %in% c("actv", "prim"),]
RNA_D0_longer$D0_class <- factor(RNA_D0_longer$D0_class, levels = c("actv", "prim"))

## Facet by enhancer class and technique
p <- ggplot(RNA_D0_longer, aes(x = meth.wmean, y = log2(WT_D0_avg + 1), color = D0_class)) +
  geom_point(na.rm = TRUE, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = c(actv = "#E5614E", prim = "#F0AC4B")) +
  geom_smooth(method = "lm", na.rm = TRUE, se = TRUE, alpha = 0.4, color = "gray30", linewidth = 0.3) +
  labs(x = "WT D0 5(h)mC%", y = "RNA log2(WT D0)") +
  coord_cartesian(xlim = c(0,100)) +
  theme_bw() + theme(legend.position = "none", panel.grid.minor = element_blank())+
  facet_grid(rows = vars(Mark), cols = vars(D0_class))
p
Name_file <- "WT_D0_hMATAC_MATAC_vs_RNAseq_byEnhClassExt_ap.png"
#ggsave(filename = Name_file, plot = p, device = "png", path = Output_folder, width = 12, height = 12, units = "cm", create.dir = TRUE)

## -------
## -- D2

RNA_D2_longer <- RNA_data_longer[RNA_data_longer$Day == "D2" & RNA_data_longer$D2_class_ext %in% c("actv", "actvM", "prim"),]
RNA_D2_longer$D2_class_ext <- factor(RNA_D2_longer$D2_class_ext, levels = c("actv", "actvM", "prim"))

## Facet by enhancer class and technique
p <- ggplot(RNA_D2_longer, aes(x = meth.wmean, y = log2(WT_D2_avg + 1), color = D2_class_ext)) +
  geom_point(na.rm = TRUE, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = c(actv = "#E5614E", actvM = "#97C124", prim = "#F0AC4B")) +
  geom_smooth(method = "lm", na.rm = TRUE, se = TRUE, alpha = 0.4, color = "gray30", linewidth = 0.3) +
  labs(x = "WT D2 5(h)mC%", y = "RNA log2(WT D2)") +
  coord_cartesian(xlim = c(0,100)) +
  theme_bw() + theme(legend.position = "none", panel.grid.minor = element_blank())+
  facet_grid(rows = vars(Mark), cols = vars(D2_class_ext))
p
Name_file <- "WT_D2_hMATAC_MATAC_vs_RNAseq_byEnhClassExt_aap.png"
#ggsave(filename = Name_file, plot = p, device = "png", path = Output_folder, width = 17, height = 12, units = "cm", create.dir = TRUE)

## -------
## -- D4

RNA_D4_longer <- RNA_data_longer[RNA_data_longer$Day == "D4" & RNA_data_longer$D4_class_ext %in% c("actv", "actvM", "prim"),]
RNA_D4_longer$D4_class_ext <- factor(RNA_D4_longer$D4_class_ext, levels = c("actv", "actvM", "prim"))

## Facet by enhancer class and technique
p <- ggplot(RNA_D4_longer, aes(x = meth.wmean, y = log2(WT_D4_avg + 1), color = D4_class_ext)) +
  geom_point(na.rm = TRUE, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = c(actv = "#E5614E", actvM = "#97C124", prim = "#F0AC4B")) +
  geom_smooth(method = "lm", na.rm = TRUE, se = TRUE, alpha = 0.4, color = "gray30", linewidth = 0.3) +
  labs(x = "WT D4 5(h)mC%", y = "RNA log2(WT D4)") +
  coord_cartesian(xlim = c(0,100)) +
  theme_bw() + theme(legend.position = "none", panel.grid.minor = element_blank())+
  facet_grid(rows = vars(Mark), cols = vars(D4_class_ext))
p
Name_file <- "WT_D4_hMATAC_MATAC_vs_RNAseq_byEnhClassExt_aap.png"
#ggsave(filename = Name_file, plot = p, device = "png", path = Output_folder, width = 17, height = 12, units = "cm", create.dir = TRUE)


## ----------------------------------------------
## ----------------------------------------------

## ---- Additional plots and analysis ---- ##

## -- Relationship 5hmC vs 5mC
## D0
ggplot(RNA_data, aes(x = meth.wmean_WT_D0_MATAC, y = meth.wmean_WT_D0_hMATAC)) +
  geom_point(na.rm = TRUE, alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", na.rm = TRUE, se = TRUE, alpha = 0.4, color = "blue", linewidth = 0.3) +
  theme_bw() +
  facet_wrap(~D0_class)

## D2
ggplot(RNA_data, aes(x = meth.wmean_WT_D2_MATAC, y = meth.wmean_WT_D2_hMATAC)) +
  geom_point(na.rm = TRUE, alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", na.rm = TRUE, se = TRUE, alpha = 0.4, color = "blue", linewidth = 0.3) +
  theme_bw() +
  facet_wrap(~D2_class_ext)

## D4
ggplot(RNA_data, aes(x = meth.wmean_WT_D4_MATAC, y = meth.wmean_WT_D4_hMATAC)) +
  geom_point(na.rm = TRUE, alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", na.rm = TRUE, se = TRUE, alpha = 0.4, color = "blue", linewidth = 0.3) +
  theme_bw() +
  facet_wrap(~D4_class_ext)


## -- Correlation coefficient 5(h)mC and RNA expression

tmp_D0 <- RNA_D0_longer |> select(chr, start, end, D0_class, gene_symbol, WT_D0_avg, Day, Mark, meth.wmean) |> filter(!is.na(WT_D0_avg))
tmp_D2 <- RNA_D2_longer |> select(chr, start, end, D2_class_ext, gene_symbol, WT_D2_avg, Day, Mark, meth.wmean) |> filter(!is.na(WT_D2_avg))
tmp_D4 <- RNA_D4_longer |> select(chr, start, end, D4_class_ext, gene_symbol, WT_D4_avg, Day, Mark, meth.wmean) |> filter(!is.na(WT_D4_avg))
tmp_D0 <- tmp_D0 |> rename(Day_class = D0_class, WT_Day_avg = WT_D0_avg)
tmp_D2 <- tmp_D2 |> rename(Day_class = D2_class_ext, WT_Day_avg = WT_D2_avg)
tmp_D4 <- tmp_D4 |> rename(Day_class = D4_class_ext, WT_Day_avg = WT_D4_avg)
tmp_allDays <- rbind(tmp_D0, tmp_D2, tmp_D4)
rm(tmp_D0, tmp_D2, tmp_D4)

## Factorize
tmp_allDays$Day_class <- factor(tmp_allDays$Day_class, levels = c("actv", "actvM", "prim"))
tmp_allDays$Day <- factor(tmp_allDays$Day, levels = c("D0", "D2", "D4"))
tmp_allDays$Mark <- factor(tmp_allDays$Mark, levels = c("hMATAC", "MATAC"))
tmp_allDays <- tmp_allDays |> group_by(Day_class, Day, Mark)

## Extract Pearson correlations
Summ_cor <- tmp_allDays |> summarise("tStat" = cor.test(meth.wmean, log2(WT_Day_avg+1))$statistic,
                                     "DFreedom" = cor.test(meth.wmean, log2(WT_Day_avg+1))$parameter,
                                     "Correlation" = cor.test(meth.wmean, log2(WT_Day_avg+1))$estimate,
                                     "CI.left" = cor.test(meth.wmean, log2(WT_Day_avg+1))$conf.int[1],
                                     "CI.right" = cor.test(meth.wmean, log2(WT_Day_avg+1))$conf.int[2],
                                     "Pval" = cor.test(meth.wmean, log2(WT_Day_avg+1))$p.value)
Summ_cor$Paj <- p.adjust(Summ_cor$Pval, method = "bonferroni")
Summ_cor$Sig <- Summ_cor$Paj < 0.05
Summ_cor <- Summ_cor |> arrange(Day)

## Make simplified version for exporting
Summ_cor_simple <- Summ_cor |> select(Day_class, Day, Mark, Correlation, CI.left, CI.right, Pval, Paj, Sig)
Summ_cor_simple$Correlation <- round(Summ_cor_simple$Correlation, digits = 4)
Summ_cor_simple$CI.left <- round(Summ_cor_simple$CI.left, digits = 4)
Summ_cor_simple$CI.right <- round(Summ_cor_simple$CI.right, digits = 4)

## Round to four significant decimals
Summ_cor_simple$Pval <- signif(Summ_cor_simple$Pval, digits = 5)
Summ_cor_simple$Paj <- signif(Summ_cor_simple$Paj, digits = 5)

## Export
#write.table(x = Summ_cor_simple, file = paste0(Output_folder,"/Info.txt"), append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)

