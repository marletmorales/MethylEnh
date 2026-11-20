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

## -- Remove inac/othr enhancers
New_Enhs <- read.table(file = "WT_Master_Enhancer_Peakome_v6_single-class_wActvM_v2.txt", header = TRUE, sep = "\t")
## Make coordinates
New_Enhs$coord <- str_c(str_c(New_Enhs$chr, New_Enhs$start, sep = ":"), New_Enhs$end, sep = "-")
New_Enhs <- New_Enhs |> dplyr::relocate(coord, .after = end)
## Filter
All_data <- All_data |> dplyr::filter(coord %in% New_Enhs$coord)
rm(New_Enhs)

## -- Plotting variables
Folder_path <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Paper_additional_files/Figures_Revision_1/Figure_4/Panel-C_Density_plots_5mC_by_Enhancer_Class/"
Color_scale <- c("actv" = "#b1c848",
                 "actvM" = "#2cb11b",
                 "prim" = "#1bb6af",
                 "inac" = "#485d7b",
                 "othr" = "#d7d7d7")

gg_layers <- list(
  geom_density(na.rm = TRUE, outline.type = "full", linewidth = 0.75, alpha = 0.1, bounds = c(0,100)),
  scale_fill_manual(values = Color_scale, drop = FALSE),
  scale_color_manual(values = Color_scale, drop = FALSE),
  coord_cartesian(xlim = c(0,100)),
  scale_x_continuous(expand = c(0, 0)),
  scale_y_continuous(expand = expansion(mult = c(0, .1))),
  theme_bw(),
  theme(legend.position = "none", panel.grid = element_blank()))

## ----------------------- ##
## -------- WT D0 -------- ##

## K4me1
p1 <- ggplot(All_data[All_data$D0_class %in% c("actv","prim","inac"),], aes(x = meth.wmean_D0_K4me1, color = D0_class, fill = D0_class)) +
  labs(x = "WT D0 K4me1 5mC %") +
  gg_layers
## MATAC
p2 <- ggplot(All_data[All_data$D0_class %in% c("actv","prim"),], aes(x = meth.wmean_D0_MATAC, color = D0_class, fill = D0_class)) +
  labs(x = "WT D0 MATAC 5mC %") +
  gg_layers
## K27ac
p3 <- ggplot(All_data[All_data$D0_class %in% c("actv"),], aes(x = meth.wmean_D0_K27ac, color = D0_class, fill = D0_class)) +
  labs(x = "WT D0 K27ac 5mC %") +
  gg_layers
## WGBS
p4 <- ggplot(All_data[All_data$D0_class %in% c("actv","prim","inac"),], aes(x = meth.wmean_D0_WGBS, color = D0_class, fill = D0_class)) +
  labs(x = "WT D0 WGBS 5mC %") +
  gg_layers

## Save individual plots
# ggsave(filename = "WT_D0_DNAme_per_enhancer_type_K4me1_v2.pdf", plot = p1, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D0_DNAme_per_enhancer_type_MATAC_v2.pdf", plot = p2, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D0_DNAme_per_enhancer_type_K27ac_v2.pdf", plot = p3, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D0_DNAme_per_enhancer_type_WGBS_v2.pdf", plot = p4, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")

## ----------------------- ##
## -------- WT D2 -------- ##

## K4me1
p1 <- ggplot(All_data[All_data$D2_class %in% c("actv","prim","inac"),], aes(x = meth.wmean_D2_K4me1, color = D2_class, fill = D2_class)) +
  labs(x = "WT D2 K4me1 5mC %") +
  gg_layers
## MATAC
p2 <- ggplot(All_data[All_data$D2_class %in% c("actv","prim"),], aes(x = meth.wmean_D2_MATAC, color = D2_class, fill = D2_class)) +
  labs(x = "WT D2 MATAC 5mC %") +
  gg_layers
## K27ac
p3 <- ggplot(All_data[All_data$D2_class %in% c("actv"),], aes(x = meth.wmean_D2_K27ac, color = D2_class, fill = D2_class)) +
  labs(x = "WT D2 K27ac 5mC %") +
  gg_layers
## WGBS
p4 <- ggplot(All_data[All_data$D2_class %in% c("actv","prim","inac"),], aes(x = meth.wmean_D2_WGBS, color = D2_class, fill = D2_class)) +
  labs(x = "WT D2 WGBS 5mC %") +
  gg_layers

## Save individual plots
# ggsave(filename = "WT_D2_DNAme_per_enhancer_type_K4me1_v2.pdf", plot = p1, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D2_DNAme_per_enhancer_type_MATAC_v2.pdf", plot = p2, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D2_DNAme_per_enhancer_type_K27ac_v2.pdf", plot = p3, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D2_DNAme_per_enhancer_type_WGBS_v2.pdf", plot = p4, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")

## ----------------------- ##
## -------- WT D4 -------- ##

## K4me1
p1 <- ggplot(All_data[All_data$D4_class %in% c("actv","prim","inac"),], aes(x = meth.wmean_D4_K4me1, color = D4_class, fill = D4_class)) +
  labs(x = "WT D4 K4me1 5mC %") +
  gg_layers
## MATAC
p2 <- ggplot(All_data[All_data$D4_class %in% c("actv","prim"),], aes(x = meth.wmean_D4_MATAC, color = D4_class, fill = D4_class)) +
  labs(x = "WT D4 MATAC 5mC %") +
  gg_layers
## K27ac
p3 <- ggplot(All_data[All_data$D4_class %in% c("actv"),], aes(x = meth.wmean_D4_K27ac, color = D4_class, fill = D4_class)) +
  labs(x = "WT D4 K27ac 5mC %") +
  gg_layers
## WGBS
p4 <- ggplot(All_data[All_data$D4_class %in% c("actv","prim","inac"),], aes(x = meth.wmean_D4_WGBS, color = D4_class, fill = D4_class)) +
  labs(x = "WT D4 WGBS 5mC %") +
  gg_layers

## Save individual plots
# ggsave(filename = "WT_D4_DNAme_per_enhancer_type_K4me1_v2.pdf", plot = p1, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D4_DNAme_per_enhancer_type_MATAC_v2.pdf", plot = p2, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D4_DNAme_per_enhancer_type_K27ac_v2.pdf", plot = p3, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
# ggsave(filename = "WT_D4_DNAme_per_enhancer_type_WGBS_v2.pdf", plot = p4, device = "pdf", path = Folder_path, width = 8, height = 8, units = "cm")
