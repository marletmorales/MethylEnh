## 25.11.2025
## Compare WT DNAme at K4me1-defined peakome (59,233 enhancers)
## Based on 06.08.2025
library(ggplot2)
library(tidyr)
library(stringr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_Master_Enhancer_Peakome_v5/")

## -- Read DNAme results
WT_D0_K4me1 <- read.table(file = "meanDNAme_WT_D0_K4me1_MCnT_mergedreps__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D2_K4me1 <- read.table(file = "meanDNAme_WT_D2_K4me1_MCnT_mergedreps__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D4_K4me1 <- read.table(file = "meanDNAme_WT_D4_K4me1_MCnT_mergedreps__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')

WT_D0_MATAC <- read.table(file = "meanDNAme_WT_D0_MATAC_mergedreps__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D2_MATAC <- read.table(file = "meanDNAme_WT_D2_MATAC_mergedreps__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D4_MATAC <- read.table(file = "meanDNAme_WT_D4_MATAC_mergedreps__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')

WT_D0_K27ac <- read.table(file = "meanDNAme_WT_D0_MCnT-EMseq_K27ac_mergedreps_0mCH__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D2_K27ac <- read.table(file = "meanDNAme_WT_D2_MCnT-EMseq_K27ac_mergedreps_0mCH__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D4_K27ac <- read.table(file = "meanDNAme_WT_D4_MCnT-EMseq_K27ac_mergedreps_0mCH__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')

WT_D0_WGBS <- read.table(file = "meanDNAme_WT_D0_WGBS__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D2_WGBS <- read.table(file = "meanDNAme_WT_D2_WGBS__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')
WT_D4_WGBS <- read.table(file = "meanDNAme_WT_D4_WGBS__Coordinates_WT_Master_Enhancer_Peakome_v5__2025_11_25.txt", header = TRUE, sep = '\t')

WT_D0_K4me1$Dataset <- "WT_D0_K4me1"
WT_D2_K4me1$Dataset <- "WT_D2_K4me1"
WT_D4_K4me1$Dataset <- "WT_D4_K4me1"

WT_D0_MATAC$Dataset <- "WT_D0_MATAC"
WT_D2_MATAC$Dataset <- "WT_D2_MATAC"
WT_D4_MATAC$Dataset <- "WT_D4_MATAC"

WT_D0_K27ac$Dataset <- "WT_D0_K27ac"
WT_D2_K27ac$Dataset <- "WT_D2_K27ac"
WT_D4_K27ac$Dataset <- "WT_D4_K27ac"

WT_D0_WGBS$Dataset <- "WT_D0_WGBS"
WT_D2_WGBS$Dataset <- "WT_D2_WGBS"
WT_D4_WGBS$Dataset <- "WT_D4_WGBS"

Peaks_All <- read.table(file = "WT_Master_Enhancer_Peakome_v5_single-class.bed", header = FALSE, sep = '\t')
colnames(Peaks_All) <- c("chr", "start", "end", "D0_class", "D2_class", "D4_class")

All_data <- rbind(WT_D0_K4me1, WT_D2_K4me1, WT_D4_K4me1,
                  WT_D0_MATAC, WT_D2_MATAC, WT_D4_MATAC,
                  WT_D0_K27ac, WT_D2_K27ac, WT_D4_K27ac,
                  WT_D0_WGBS, WT_D2_WGBS, WT_D4_WGBS)

## Convert dataset to factor
All_data$Dataset <- factor(All_data$Dataset, levels = c("WT_D0_K4me1", "WT_D2_K4me1", "WT_D4_K4me1",
                                                        "WT_D0_MATAC", "WT_D2_MATAC", "WT_D4_MATAC",
                                                        "WT_D0_K27ac", "WT_D2_K27ac", "WT_D4_K27ac",
                                                        "WT_D0_WGBS", "WT_D2_WGBS", "WT_D4_WGBS"))
## Add annotation from peaks
All_data <- merge(x = All_data, y = Peaks_All, all.x = TRUE, all.y = FALSE)
## Convert annotation to factors
All_data$D0_class <- factor(All_data$D0_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D2_class <- factor(All_data$D2_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D4_class <- factor(All_data$D4_class, levels = c("actv", "prim", "inac", "othr"))
## Add coordinate
All_data$coord <- str_c(str_c(All_data$chr, All_data$start, sep = ":"), All_data$end, sep = "-")

## Pivot
All_data <- pivot_wider(All_data, names_from = "Dataset", values_from = c("num.cpgs", "meth.mean", "c.meth", "c.unmeth", "meth.wmean"))
length(unique(All_data$coord))

## Reorder columns
All_data <- All_data[, c("chr", "start", "end", "coord", "D0_class", "D2_class", "D4_class",
                         "num.cpgs_WT_D0_K4me1", "meth.mean_WT_D0_K4me1", "c.meth_WT_D0_K4me1", "c.unmeth_WT_D0_K4me1", "meth.wmean_WT_D0_K4me1",
                         "num.cpgs_WT_D2_K4me1", "meth.mean_WT_D2_K4me1", "c.meth_WT_D2_K4me1", "c.unmeth_WT_D2_K4me1", "meth.wmean_WT_D2_K4me1",
                         "num.cpgs_WT_D4_K4me1", "meth.mean_WT_D4_K4me1", "c.meth_WT_D4_K4me1", "c.unmeth_WT_D4_K4me1", "meth.wmean_WT_D4_K4me1",
                         "num.cpgs_WT_D0_MATAC", "meth.mean_WT_D0_MATAC", "c.meth_WT_D0_MATAC", "c.unmeth_WT_D0_MATAC", "meth.wmean_WT_D0_MATAC",
                         "num.cpgs_WT_D2_MATAC", "meth.mean_WT_D2_MATAC", "c.meth_WT_D2_MATAC", "c.unmeth_WT_D2_MATAC", "meth.wmean_WT_D2_MATAC",
                         "num.cpgs_WT_D4_MATAC", "meth.mean_WT_D4_MATAC", "c.meth_WT_D4_MATAC", "c.unmeth_WT_D4_MATAC", "meth.wmean_WT_D4_MATAC",
                         "num.cpgs_WT_D0_K27ac", "meth.mean_WT_D0_K27ac", "c.meth_WT_D0_K27ac", "c.unmeth_WT_D0_K27ac", "meth.wmean_WT_D0_K27ac",
                         "num.cpgs_WT_D2_K27ac", "meth.mean_WT_D2_K27ac", "c.meth_WT_D2_K27ac", "c.unmeth_WT_D2_K27ac", "meth.wmean_WT_D2_K27ac",
                         "num.cpgs_WT_D4_K27ac", "meth.mean_WT_D4_K27ac", "c.meth_WT_D4_K27ac", "c.unmeth_WT_D4_K27ac", "meth.wmean_WT_D4_K27ac",
                         "num.cpgs_WT_D0_WGBS", "meth.mean_WT_D0_WGBS", "c.meth_WT_D0_WGBS", "c.unmeth_WT_D0_WGBS", "meth.wmean_WT_D0_WGBS",
                         "num.cpgs_WT_D2_WGBS", "meth.mean_WT_D2_WGBS", "c.meth_WT_D2_WGBS", "c.unmeth_WT_D2_WGBS", "meth.wmean_WT_D2_WGBS",
                         "num.cpgs_WT_D4_WGBS", "meth.mean_WT_D4_WGBS", "c.meth_WT_D4_WGBS", "c.unmeth_WT_D4_WGBS", "meth.wmean_WT_D4_WGBS")]

rm(list = c("WT_D0_K4me1", "WT_D2_K4me1", "WT_D4_K4me1",
            "WT_D0_MATAC", "WT_D2_MATAC", "WT_D4_MATAC",
            "WT_D0_K27ac", "WT_D2_K27ac", "WT_D4_K27ac",
            "WT_D0_WGBS", "WT_D2_WGBS", "WT_D4_WGBS"))

# ## Export table
# write.table(x = All_data, file = "meanDNAme_WT_D0-D2-D4_K4me1-MATAC-K27ac-WGBS__Coordinates_WT_Master_Enhancer_Peakome_v5.txt",
#             append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE)

Folder_path <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/Measure_DNAme_ROIs/WT_Master_Enhancer_Peakome_v5/"

## ----------------------- ##
## -------- WT D0 -------- ##

Sub_data <- All_data[All_data$D0_class != "othr",]
## K4me1
p1 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_WT_D0_K4me1, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 K4me1 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## MATAC
p2 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv","prim"),],
             aes(x = meth.wmean_WT_D0_MATAC, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 MATAC 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## K27ac
p3 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv"),],
             aes(x = meth.wmean_WT_D0_K27ac, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 K27ac 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## WGBS
p4 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_WT_D0_WGBS, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 WGBS 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)

D0_alltechs <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)
ggsave(filename = "WT_D0_DNAme_per_enhancer_type.pdf", plot = D0_alltechs, device = "pdf", path = Folder_path, width = 26, height = 7, units = "cm")

## Save individual plots
ggsave(filename = "WT_D0_DNAme_per_enhancer_type_K4me1.pdf", plot = p1, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D0_DNAme_per_enhancer_type_MATAC.pdf", plot = p2, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D0_DNAme_per_enhancer_type_K27ac.pdf", plot = p3, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D0_DNAme_per_enhancer_type_WGBS.pdf", plot = p4, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")

table(Sub_data$D0_class)
colSums(!is.na(Sub_data[Sub_data$D0_class=="actv",c("meth.wmean_WT_D0_K4me1", "meth.wmean_WT_D0_MATAC", "meth.wmean_WT_D0_K27ac", "meth.wmean_WT_D0_WGBS")]))
colSums(!is.na(Sub_data[Sub_data$D0_class=="prim",c("meth.wmean_WT_D0_K4me1", "meth.wmean_WT_D0_MATAC", "meth.wmean_WT_D0_WGBS")]))
colSums(!is.na(Sub_data[Sub_data$D0_class=="inac",c("meth.wmean_WT_D0_K4me1", "meth.wmean_WT_D0_WGBS")]))


## ----------------------- ##
## -------- WT D2 -------- ##

Sub_data <- All_data[All_data$D2_class != "othr",]
## K4me1
p1 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_WT_D2_K4me1, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 K4me1 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## MATAC
p2 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv","prim"),],
             aes(x = meth.wmean_WT_D2_MATAC, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 MATAC 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## K27ac
p3 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv"),],
             aes(x = meth.wmean_WT_D2_K27ac, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 K27ac 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## WGBS
p4 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_WT_D2_WGBS, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 WGBS 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)

D2_alltechs <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)
ggsave(filename = "WT_D2_DNAme_per_enhancer_type.pdf", plot = D2_alltechs, device = "pdf", path = Folder_path, width = 26, height = 7, units = "cm")

## Save individual plots
ggsave(filename = "WT_D2_DNAme_per_enhancer_type_K4me1.pdf", plot = p1, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D2_DNAme_per_enhancer_type_MATAC.pdf", plot = p2, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D2_DNAme_per_enhancer_type_K27ac.pdf", plot = p3, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D2_DNAme_per_enhancer_type_WGBS.pdf", plot = p4, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")

table(Sub_data$D2_class)
colSums(!is.na(Sub_data[Sub_data$D2_class=="actv",c("meth.wmean_WT_D2_K4me1", "meth.wmean_WT_D2_MATAC", "meth.wmean_WT_D2_K27ac", "meth.wmean_WT_D2_WGBS")]))
colSums(!is.na(Sub_data[Sub_data$D2_class=="prim",c("meth.wmean_WT_D2_K4me1", "meth.wmean_WT_D2_MATAC", "meth.wmean_WT_D2_WGBS")]))
colSums(!is.na(Sub_data[Sub_data$D2_class=="inac",c("meth.wmean_WT_D2_K4me1", "meth.wmean_WT_D2_WGBS")]))


## ----------------------- ##
## -------- WT D4 -------- ##

Sub_data <- All_data[All_data$D4_class != "othr",]
## K4me1
p1 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_WT_D4_K4me1, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 K4me1 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## MATAC
p2 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv","prim"),],
             aes(x = meth.wmean_WT_D4_MATAC, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 MATAC 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## K27ac
p3 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv"),],
             aes(x = meth.wmean_WT_D4_K27ac, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 K27ac 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## WGBS
p4 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_WT_D4_WGBS, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 WGBS 5mC") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)

D4_alltechs <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)
ggsave(filename = "WT_D4_DNAme_per_enhancer_type.pdf", plot = D4_alltechs, device = "pdf", path = Folder_path, width = 26, height = 7, units = "cm")

## Save individual plots
ggsave(filename = "WT_D4_DNAme_per_enhancer_type_K4me1.pdf", plot = p1, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D4_DNAme_per_enhancer_type_MATAC.pdf", plot = p2, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D4_DNAme_per_enhancer_type_K27ac.pdf", plot = p3, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")
ggsave(filename = "WT_D4_DNAme_per_enhancer_type_WGBS.pdf", plot = p4, device = "pdf", path = Folder_path, width = 12, height = 12, units = "cm")

table(Sub_data$D4_class)
colSums(!is.na(Sub_data[Sub_data$D4_class=="actv",c("meth.wmean_WT_D4_K4me1", "meth.wmean_WT_D4_MATAC", "meth.wmean_WT_D4_K27ac", "meth.wmean_WT_D4_WGBS")]))
colSums(!is.na(Sub_data[Sub_data$D4_class=="prim",c("meth.wmean_WT_D4_K4me1", "meth.wmean_WT_D4_MATAC", "meth.wmean_WT_D4_WGBS")]))
colSums(!is.na(Sub_data[Sub_data$D4_class=="inac",c("meth.wmean_WT_D4_K4me1", "meth.wmean_WT_D4_WGBS")]))
