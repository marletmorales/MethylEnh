## 22.12.2025
## Plotting 5hmC hMATAC_D2 vs hMATAC_D4
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_hMATAC_peakome_RmProm/")

## -- Read 5hmC data
WT_D0_hMATAC <- read.table(file = "meanDNAme_WT_D0_hMATAC_mergedreps__Coordinates_WT_D0-D2-D4_hMATAC_peakome_RmProm__2025_12_22.txt", header = TRUE, sep = '\t')
WT_D2_hMATAC <- read.table(file = "meanDNAme_WT_D2_hMATAC_mergedreps__Coordinates_WT_D0-D2-D4_hMATAC_peakome_RmProm__2025_12_22.txt", header = TRUE, sep = '\t')
WT_D4_hMATAC <- read.table(file = "meanDNAme_WT_D4_hMATAC_mergedreps__Coordinates_WT_D0-D2-D4_hMATAC_peakome_RmProm__2025_12_22.txt", header = TRUE, sep = '\t')
WT_D0_hMATAC$Dataset <- "WT_D0_hMATAC"
WT_D2_hMATAC$Dataset <- "WT_D2_hMATAC"
WT_D4_hMATAC$Dataset <- "WT_D4_hMATAC"

## -- Merge all data
All_data <- rbind(WT_D0_hMATAC, WT_D2_hMATAC, WT_D4_hMATAC)
All_data <- All_data[, c("chr", "start", "end", "num.cpgs", "meth.wmean", "Dataset")]
All_data <- pivot_wider(All_data, names_from = "Dataset", values_from = c("num.cpgs", "meth.wmean"))
All_data$coord <- str_c(str_c(All_data$chr, All_data$start, sep = ":"), All_data$end, sep = "-")
length(unique(All_data$coord))
rm(list = c("WT_D0_hMATAC", "WT_D2_hMATAC", "WT_D4_hMATAC"))

## ---- Plot 5hmC: D2_hMATAC vs D4_hMATAC (NonProm peaks)
#Folder_path <- "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_5/Panel-D_scatterplot_hMATAC_D2_vs_D4/"

plot <- ggplot(All_data, aes(x = meth.wmean_WT_D2_hMATAC, y = meth.wmean_WT_D4_hMATAC)) +
  stat_bin2d(na.rm = TRUE, bins = 50) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,2000), oob = scales::squish) +
  geom_abline(intercept = 0, slope = 1, color = 'gray40') +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "5hmC hMATAC - D2", y = "5hmC hMATAC - D4") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
plot
# ## png
# ggsave(filename = "Adj_WT_D2_hMATAC-5hmC_v_WT_D4_hMATAC-5hmC__Coordinates_WT_hMATAC_NonPromPeakome.png",
#        plot = plot, device = "png", path = Folder_path,
#        width = 10, height = 10, units = "cm")
# ## pdf
# ggsave(filename = "Adj_WT_D2_hMATAC-5hmC_v_WT_D4_hMATAC-5hmC__Coordinates_WT_hMATAC_NonPromPeakome.pdf",
#        plot = plot, device = "pdf", path = Folder_path,
#        width = 10, height = 10, units = "cm")

## Number of peaks and pearson correlation
nrow(All_data |> filter(!is.na(meth.wmean_WT_D2_hMATAC) & !is.na(meth.wmean_WT_D4_hMATAC)))
cor(x = All_data$meth.wmean_WT_D2_hMATAC, y = All_data$meth.wmean_WT_D4_hMATAC, use = "pairwise.complete.obs", method = "pearson")
cor(x = All_data$meth.wmean_WT_D2_hMATAC, y = All_data$meth.wmean_WT_D4_hMATAC, use = "pairwise.complete.obs", method = "spearman")
# 63,514 peaks in common
# pearson = 0.2860364
# spearman = 0.3446889

## Export color scale
p <- ggplot(All_data, aes(x = meth.wmean_WT_D2_hMATAC, y = meth.wmean_WT_D4_hMATAC)) +
  stat_bin2d(na.rm = TRUE, bins = 50) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,2000), oob = scales::squish) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right")
p <- cowplot::ggdraw(cowplot::get_legend(p))
p
# ggsave(filename = "Color_scale_hMATAC-5hmC_D2_vs_D4.pdf",
#        plot = p, device = "pdf", path = Folder_path, width = 3, height = 5, units = "cm")

