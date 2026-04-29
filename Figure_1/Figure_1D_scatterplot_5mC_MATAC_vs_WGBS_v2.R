## 26.08.2025
## Plotting DNAme MATAC vs WGBS
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_MATAC_NonPromPeakome/")

## -- Read DNAme data
WT_D0_MATAC <- read.table(file = "meanDNAme_WT_D0_MATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", header = TRUE, sep = '\t')
WT_D2_MATAC <- read.table(file = "meanDNAme_WT_D2_MATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", header = TRUE, sep = '\t')
WT_D4_MATAC <- read.table(file = "meanDNAme_WT_D4_MATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", header = TRUE, sep = '\t')

WT_D0_WGBS <- read.table(file = "meanDNAme_WT_D0_WGBS__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", header = TRUE, sep = '\t')
WT_D2_WGBS <- read.table(file = "meanDNAme_WT_D2_WGBS__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", header = TRUE, sep = '\t')
WT_D4_WGBS <- read.table(file = "meanDNAme_WT_D4_WGBS__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", header = TRUE, sep = '\t')

WT_D0_MATAC$Dataset <- "WT_D0_MATAC"
WT_D2_MATAC$Dataset <- "WT_D2_MATAC"
WT_D4_MATAC$Dataset <- "WT_D4_MATAC"
WT_D0_WGBS$Dataset <- "WT_D0_WGBS"
WT_D2_WGBS$Dataset <- "WT_D2_WGBS"
WT_D4_WGBS$Dataset <- "WT_D4_WGBS"

## -- Merge all data
All_data <- rbind(WT_D0_MATAC, WT_D2_MATAC, WT_D4_MATAC, WT_D0_WGBS, WT_D2_WGBS, WT_D4_WGBS)
All_data <- All_data[, c("chr", "start", "end", "num.cpgs", "meth.wmean", "Dataset")]
All_data <- pivot_wider(All_data, names_from = "Dataset", values_from = c("num.cpgs", "meth.wmean"))
All_data$coord <- str_c(str_c(All_data$chr, All_data$start, sep = ":"), All_data$end, sep = "-")
length(unique(All_data$coord))
rm(list = c("WT_D0_MATAC", "WT_D2_MATAC", "WT_D4_MATAC", "WT_D0_WGBS", "WT_D2_WGBS", "WT_D4_WGBS"))

## -- Read Differentially-accessible peaks (44K)
Peaks_DA <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT_MATAC_v1/Dynamic/K5/bed/All_clusters_labeled.bed", header = FALSE, sep = '\t')
colnames(Peaks_DA) <- c("chr","start","end","cluster")

## -- Add peak and cluster information to All_data
All_data <- merge(x = All_data, y = Peaks_DA, all.x = TRUE)
All_data$cluster[is.na(All_data$cluster)] <- "none"

## Reorder columns
All_data <- All_data[, c("chr", "start", "end", "coord", "cluster",
                         "num.cpgs_WT_D0_MATAC", "meth.wmean_WT_D0_MATAC",
                         "num.cpgs_WT_D2_MATAC", "meth.wmean_WT_D2_MATAC",
                         "num.cpgs_WT_D4_MATAC", "meth.wmean_WT_D4_MATAC",
                         "num.cpgs_WT_D0_WGBS", "meth.wmean_WT_D0_WGBS",
                         "num.cpgs_WT_D2_WGBS", "meth.wmean_WT_D2_WGBS",
                         "num.cpgs_WT_D4_WGBS", "meth.wmean_WT_D4_WGBS")]


## --------------------------------------------------- ##
## -- Plot DNAme: MATAC vs WGBS (all NonProm peaks) -- ##
## --------------------------------------------------- ##

#Folder_path <- "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Scatterplot_5mC_MATACvWGBS/"

## ------ Day 0 ------
## -- Adjusted scale
D0_adj <- ggplot(All_data, aes(x = meth.wmean_WT_D0_MATAC, y = meth.wmean_WT_D0_WGBS)) +
  stat_bin2d(na.rm = TRUE, bins = 50) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(10,400), oob = scales::squish) +
  geom_abline(intercept = 0, slope = 1, color = 'gray40') +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "5mC MATAC - D0", y = "5mC WGBS - D0") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                     axis.title = element_blank(), axis.text = element_blank())
D0_adj

## png
ggsave(filename = "Adj_WT_D0_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.png",
       plot = D0_adj, device = "png",
       path = Folder_path,
       width = 10, height = 10, units = "cm")
## pdf
ggsave(filename = "Adj_WT_D0_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.pdf",
       plot = D0_adj, device = "pdf",
       path = Folder_path,
       width = 10, height = 10, units = "cm")

## ------ Day 2 ------
## -- Adjusted scale
D2_adj <- ggplot(All_data, aes(x = meth.wmean_WT_D2_MATAC, y = meth.wmean_WT_D2_WGBS)) +
  stat_bin2d(na.rm = TRUE, bins = 50) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(10,400), oob = scales::squish) +
  geom_abline(intercept = 0, slope = 1, color = 'gray40') +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "5mC MATAC - D2", y = "5mC WGBS - D2") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                     axis.title = element_blank(), axis.text = element_blank())
D2_adj

## Number of peaks and pearson correlation
nrow(All_data |> filter(!is.na(meth.wmean_WT_D2_MATAC) & !is.na(meth.wmean_WT_D2_WGBS)))
cor(x = All_data$meth.wmean_WT_D2_MATAC, y = All_data$meth.wmean_WT_D2_WGBS, use = "pairwise.complete.obs", method = "pearson")

## png
ggsave(filename = "Adj_WT_D2_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.png",
       plot = D2_adj, device = "png",
       path = Folder_path,
       width = 10, height = 10, units = "cm")
## pdf
ggsave(filename = "Adj_WT_D2_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.pdf",
       plot = D2_adj, device = "pdf",
       path = Folder_path,
       width = 10, height = 10, units = "cm")


## ------ Day 4 ------
## -- Adjusted scale
D4_adj <- ggplot(All_data, aes(x = meth.wmean_WT_D4_MATAC, y = meth.wmean_WT_D4_WGBS)) +
  stat_bin2d(na.rm = TRUE, bins = 50) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(10,400), oob = scales::squish) +
  geom_abline(intercept = 0, slope = 1, color = 'gray40') +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "5mC MATAC - D4", y = "5mC WGBS - D4") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                     axis.title = element_blank(), axis.text = element_blank())
D4_adj

## Number of peaks and pearson correlation
nrow(All_data |> filter(!is.na(meth.wmean_WT_D4_MATAC) & !is.na(meth.wmean_WT_D4_WGBS)))
cor(x = All_data$meth.wmean_WT_D4_MATAC, y = All_data$meth.wmean_WT_D4_WGBS, use = "pairwise.complete.obs", method = "pearson")

## png
ggsave(filename = "Adj_WT_D4_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.png",
       plot = D4_adj, device = "png",
       path = Folder_path,
       width = 10, height = 10, units = "cm")
## pdf
ggsave(filename = "Adj_WT_D4_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.pdf",
       plot = D4_adj, device = "pdf",
       path = Folder_path,
       width = 10, height = 10, units = "cm")


## Export color scale
p <- ggplot(All_data, aes(x = meth.wmean_WT_D0_MATAC, y = meth.wmean_WT_D0_WGBS)) +
  stat_bin2d(na.rm = TRUE, bins = 50) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(10,400), oob = scales::squish) +
  labs(x = "5mC MATAC - D0", y = "5mC WGBS - D0") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right")
p <- cowplot::ggdraw(cowplot::get_legend(p))
p
ggsave(filename = "Color_scale_MATAC-5mC_v_WGBS_5mC__Coordinates_WT_MATAC_NonPromPeakome.pdf",
       plot = p, device = "pdf", path = Folder_path, width = 3, height = 5, units = "cm")


## ----------------------------------------------------------- ##
## -- Plot DNAme: MATAC vs WGBS (all NonProm peaks, by CpG) -- ##
## ----------------------------------------------------------- ##

## Nope, it doesn't look better. Ignore.

# ## --------------------- ##
# ## -- DNAme - per CpG -- ##
# ## --------------------- ##
#
# ## -- Read DNAme data
# WT_D0_MATAC <- read.table(file = "WT_D0_MATAC_mergedreps_trimV6_mm10_cov5_i_WT_MATAC_macs2_broad_D0toD4_nonprom.txt", header = FALSE, sep = '\t')[,1:4]
# WT_D2_MATAC <- read.table(file = "WT_D2_MATAC_mergedreps_trimV6_mm10_cov5_i_WT_MATAC_macs2_broad_D0toD4_nonprom.txt", header = FALSE, sep = '\t')[,1:4]
# WT_D4_MATAC <- read.table(file = "WT_D4_MATAC_mergedreps_trimV6_mm10_cov5_i_WT_MATAC_macs2_broad_D0toD4_nonprom.txt", header = FALSE, sep = '\t')[,1:4]
#
# WT_D0_WGBS <- read.table(file = "WT_D0_WGBS_cov5_i_WT_MATAC_macs2_broad_D0toD4_nonprom.txt", header = FALSE, sep = '\t')[,1:4]
# WT_D2_WGBS <- read.table(file = "WT_D2_WGBS_cov5_i_WT_MATAC_macs2_broad_D0toD4_nonprom.txt", header = FALSE, sep = '\t')[,1:4]
# WT_D4_WGBS <- read.table(file = "WT_D4_WGBS_cov5_i_WT_MATAC_macs2_broad_D0toD4_nonprom.txt", header = FALSE, sep = '\t')[,1:4]
#
# colnames(WT_D0_MATAC) <- c("chr", "start", "end", "meth.wmean")
# colnames(WT_D2_MATAC) <- c("chr", "start", "end", "meth.wmean")
# colnames(WT_D4_MATAC) <- c("chr", "start", "end", "meth.wmean")
# colnames(WT_D0_WGBS) <- c("chr", "start", "end", "meth.wmean")
# colnames(WT_D2_WGBS) <- c("chr", "start", "end", "meth.wmean")
# colnames(WT_D4_WGBS) <- c("chr", "start", "end", "meth.wmean")
#
# WT_D0_MATAC$Dataset <- "WT_D0_MATAC"
# WT_D2_MATAC$Dataset <- "WT_D2_MATAC"
# WT_D4_MATAC$Dataset <- "WT_D4_MATAC"
# WT_D0_WGBS$Dataset <- "WT_D0_WGBS"
# WT_D2_WGBS$Dataset <- "WT_D2_WGBS"
# WT_D4_WGBS$Dataset <- "WT_D4_WGBS"
#
# ## -- Merge all data
# All_CpG <- rbind(WT_D0_MATAC, WT_D2_MATAC, WT_D4_MATAC, WT_D0_WGBS, WT_D2_WGBS, WT_D4_WGBS)
# All_CpG <- pivot_wider(All_CpG, names_from = "Dataset", values_from = c("meth.wmean"))
#
# All_CpG$coord <- str_c(str_c(All_CpG$chr, All_CpG$start, sep = ":"), All_CpG$end, sep = "-")
# length(unique(All_CpG$coord))
# rm(list = c("WT_D0_MATAC", "WT_D2_MATAC", "WT_D4_MATAC", "WT_D0_WGBS", "WT_D2_WGBS", "WT_D4_WGBS"))
#
#
# ## ------ Day 0 ------
# ggplot(All_CpG, aes(x = WT_D0_MATAC, y = WT_D0_WGBS)) +
#   stat_bin2d(aes(fill = after_stat(density)), na.rm = TRUE, bins = 50) +
#   scale_fill_viridis_c(option = "rocket", direction = -1) +
#   geom_abline(intercept = 0, slope = 1, color = 'gray40') +
#   coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
#   labs(x = "5mC MATAC - D0", y = "5mC WGBS - D0") +
#   theme_bw() + theme(legend.position = "none")
#
# # ggplot(All_CpG, aes(x = WT_D0_MATAC, y = WT_D0_WGBS)) +
# #   stat_density2d(na.rm = TRUE, bins = 100, contour = TRUE) +
# #   geom_abline(intercept = 0, slope = 1, color = 'gray40') +
# #   coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
# #   labs(x = "5mC MATAC - D0", y = "5mC WGBS - D0") +
# #   theme_bw() + theme(legend.position = "none")
#
#
# ## ------ Day 2 ------
# ggplot(All_CpG, aes(x = WT_D2_MATAC, y = WT_D2_WGBS)) +
#   stat_bin2d(aes(fill = after_stat(density)), na.rm = TRUE, bins = 50) +
#   scale_fill_viridis_c(option = "rocket", direction = -1) +
#   geom_abline(intercept = 0, slope = 1, color = 'gray40') +
#   coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
#   labs(x = "5mC MATAC - D2", y = "5mC WGBS - D2") +
#   theme_bw() + theme(legend.position = "none")
#
# ggplot(All_CpG, aes(x = WT_D2_MATAC, y = WT_D2_WGBS)) +
#   stat_density2d(na.rm = TRUE, bins = 40, contour = TRUE) +
#   geom_abline(intercept = 0, slope = 1, color = 'gray40') +
#   coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
#   labs(x = "5mC MATAC - D2", y = "5mC WGBS - D2") +
#   theme_bw() + theme(legend.position = "none")
#
#
# ## ------ Day 4 ------
# ggplot(All_CpG, aes(x = WT_D4_MATAC, y = WT_D4_WGBS)) +
#   stat_bin2d(aes(fill = after_stat(density)), na.rm = TRUE, bins = 50) +
#   scale_fill_viridis_c(option = "rocket", direction = -1) +
#   geom_abline(intercept = 0, slope = 1, color = 'gray40') +
#   coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
#   labs(x = "5mC MATAC - D4", y = "5mC WGBS - D4") +
#   theme_bw() + theme(legend.position = "none")
#
# ggplot(All_CpG, aes(x = WT_D4_MATAC, y = WT_D4_WGBS)) +
#   stat_density2d(na.rm = TRUE, bins = 40, contour = TRUE) +
#   geom_abline(intercept = 0, slope = 1, color = 'gray40') +
#   coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
#   labs(x = "5mC MATAC - D4", y = "5mC WGBS - D4") +
#   theme_bw() + theme(legend.position = "none")
