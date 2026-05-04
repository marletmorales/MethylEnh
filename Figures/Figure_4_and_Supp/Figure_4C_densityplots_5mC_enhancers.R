## Density plots of 5mC at different enhancer types
## Figure 4C and Supplemental Figure 5C
library(ggplot2)
library(tidyr)
library(stringr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/")

## -- Read DNAme results
All_data <- read.table(file = "WT_DNAme_allSamples__Coordinates_WT_Master_Enhancer_Peakome_v6.txt", header = TRUE, sep = "\t")
## Convert annotation to factors
All_data$D0_class <- factor(All_data$D0_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D2_class <- factor(All_data$D2_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D4_class <- factor(All_data$D4_class, levels = c("actv", "prim", "inac", "othr"))

#Folder_path <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/Measure_DNAme_ROIs/WT_Master_Enhancer_Peakome_v6/"

## ----------------------- ##
## -------- WT D0 -------- ##

Sub_data <- All_data[All_data$D0_class != "othr",]
## K4me1
p1 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_D0_K4me1, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 K4me1 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## MATAC
p2 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv","prim"),],
             aes(x = meth.wmean_D0_MATAC, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 MATAC 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## K27ac
p3 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv"),],
             aes(x = meth.wmean_D0_K27ac, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 K27ac 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## WGBS
p4 <- ggplot(Sub_data[Sub_data$D0_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_D0_WGBS, color = D0_class, fill = D0_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D0 WGBS 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

D0_alltechs <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)
D0_alltechs
#ggsave(filename = "WT_D0_DNAme_per_enhancer_type.pdf", plot = D0_alltechs, device = "pdf", path = Folder_path, width = 26, height = 7, units = "cm")

## Save individual plots
#ggsave(filename = "WT_D0_DNAme_per_enhancer_type_K4me1.pdf", plot = p1, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D0_DNAme_per_enhancer_type_MATAC.pdf", plot = p2, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D0_DNAme_per_enhancer_type_K27ac.pdf", plot = p3, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D0_DNAme_per_enhancer_type_WGBS.pdf", plot = p4, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")

## ----------------------- ##
## -------- WT D2 -------- ##

Sub_data <- All_data[All_data$D2_class != "othr",]
## K4me1
p1 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_D2_K4me1, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 K4me1 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## MATAC
p2 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv","prim"),],
             aes(x = meth.wmean_D2_MATAC, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 MATAC 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## K27ac
p3 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv"),],
             aes(x = meth.wmean_D2_K27ac, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 K27ac 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## WGBS
p4 <- ggplot(Sub_data[Sub_data$D2_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_D2_WGBS, color = D2_class, fill = D2_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D2 WGBS 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

D2_alltechs <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)
D2_alltechs
#ggsave(filename = "WT_D2_DNAme_per_enhancer_type.pdf", plot = D2_alltechs, device = "pdf", path = Folder_path, width = 26, height = 7, units = "cm")

## Save individual plots
#ggsave(filename = "WT_D2_DNAme_per_enhancer_type_K4me1.pdf", plot = p1, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D2_DNAme_per_enhancer_type_MATAC.pdf", plot = p2, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D2_DNAme_per_enhancer_type_K27ac.pdf", plot = p3, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D2_DNAme_per_enhancer_type_WGBS.pdf", plot = p4, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")

## ----------------------- ##
## -------- WT D4 -------- ##

Sub_data <- All_data[All_data$D4_class != "othr",]
## K4me1
p1 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_D4_K4me1, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 K4me1 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## MATAC
p2 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv","prim"),],
             aes(x = meth.wmean_D4_MATAC, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 MATAC 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## K27ac
p3 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv"),],
             aes(x = meth.wmean_D4_K27ac, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 K27ac 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
## WGBS
p4 <- ggplot(Sub_data[Sub_data$D4_class %in% c("actv","prim","inac"),],
             aes(x = meth.wmean_D4_WGBS, color = D4_class, fill = D4_class)) +
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  scale_color_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
  coord_cartesian(xlim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(x = "WT D4 WGBS 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

D4_alltechs <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 1)
D4_alltechs
#ggsave(filename = "WT_D4_DNAme_per_enhancer_type.pdf", plot = D4_alltechs, device = "pdf", path = Folder_path, width = 26, height = 7, units = "cm")

## Save individual plots
#ggsave(filename = "WT_D4_DNAme_per_enhancer_type_K4me1.pdf", plot = p1, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D4_DNAme_per_enhancer_type_MATAC.pdf", plot = p2, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D4_DNAme_per_enhancer_type_K27ac.pdf", plot = p3, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")
#ggsave(filename = "WT_D4_DNAme_per_enhancer_type_WGBS.pdf", plot = p4, device = "pdf", path = Folder_path, width = 10, height = 10, units = "cm")

