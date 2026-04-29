## Plotting DNAme of MATAC and WGBS, separated by cluster
## Figure 1E and Supplemental Figure 1E
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

## -- Read dynamic peaks
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

## ----------------------------------
## ---- Barplots clustered peaks ----

Clust_data <- filter(All_data, cluster != "none")
## Remove unnecessary columns
Clust_data <- Clust_data[,c("chr", "start", "end", "coord", "cluster",
                            "meth.wmean_WT_D0_MATAC", "meth.wmean_WT_D2_MATAC", "meth.wmean_WT_D4_MATAC",
                            "meth.wmean_WT_D0_WGBS", "meth.wmean_WT_D2_WGBS", "meth.wmean_WT_D4_WGBS")]
## Pivot from wide to long
Clust_data <- pivot_longer(data = Clust_data, cols = c("meth.wmean_WT_D0_MATAC", "meth.wmean_WT_D2_MATAC", "meth.wmean_WT_D4_MATAC", "meth.wmean_WT_D0_WGBS", "meth.wmean_WT_D2_WGBS", "meth.wmean_WT_D4_WGBS"),
                           names_to = "Dataset", values_to = "meth.wmean", values_drop_na = TRUE)

## Split $Dataset into two: day and technique
Day_tech <- gsub(pattern = "meth.wmean_WT_", replacement = "", x = Clust_data$Dataset)
Day_tech <- str_split_fixed(string = Day_tech, pattern = "_", n = 2)
Clust_data$Day <- Day_tech[,1]
Clust_data$Tech <- Day_tech[,2]

## -- Plot per cluster

p1 <- ggplot(Clust_data[Clust_data$cluster == "c1",], aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(y = "% DNAme", x = "c1 - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p1
#ggsave(filename = "Cluster_1_boxplot.pdf", plot = p1, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 8, height = 8, units = "cm")

p2 <- ggplot(Clust_data[Clust_data$cluster == "c2",], aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(y = "% DNAme", x = "c2 - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p2
#ggsave(filename = "Cluster_2_boxplot.pdf", plot = p2, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 8, height = 8, units = "cm")

p3 <- ggplot(Clust_data[Clust_data$cluster == "c3",], aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(y = "% DNAme", x = "c3 - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p3
#ggsave(filename = "Cluster_3_boxplot.pdf", plot = p3, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 8, height = 8, units = "cm")

p4 <- ggplot(Clust_data[Clust_data$cluster == "c4",], aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(y = "% DNAme", x = "c4 - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p4
#ggsave(filename = "Cluster_4_boxplot.pdf", plot = p4, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 8, height = 8, units = "cm")

p5 <- ggplot(Clust_data[Clust_data$cluster == "c5",], aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(y = "% DNAme", x = "c5 - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p5
#ggsave(filename = "Cluster_5_boxplot.pdf", plot = p5, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 8, height = 8, units = "cm")


## Get color legend
p <- ggplot(Clust_data[Clust_data$cluster == "c1",], aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  labs(fill = "Technique") +
  theme_bw() + theme(legend.position = "right", panel.grid = element_blank())
p <- cowplot::ggdraw(cowplot::get_legend(p))
p
#ggsave(filename = "Boxplot_legend.pdf", plot = p, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 4, height = 4, units = "cm")


## ---------------------------------
## ---- Non-clustered peaks (static)

Static_data <- filter(All_data, cluster == "none")
## Remove unnecessary columns
Static_data <- Static_data[,c("chr", "start", "end", "coord", "cluster",
                            "meth.wmean_WT_D0_MATAC", "meth.wmean_WT_D2_MATAC", "meth.wmean_WT_D4_MATAC",
                            "meth.wmean_WT_D0_WGBS", "meth.wmean_WT_D2_WGBS", "meth.wmean_WT_D4_WGBS")]
## Pivot from wide to long
Static_data <- pivot_longer(data = Static_data, cols = c("meth.wmean_WT_D0_MATAC", "meth.wmean_WT_D2_MATAC", "meth.wmean_WT_D4_MATAC", "meth.wmean_WT_D0_WGBS", "meth.wmean_WT_D2_WGBS", "meth.wmean_WT_D4_WGBS"),
                           names_to = "Dataset", values_to = "meth.wmean", values_drop_na = TRUE)

## Split $Dataset into two: day and technique
Day_tech <- gsub(pattern = "meth.wmean_WT_", replacement = "", x = Static_data$Dataset)
Day_tech <- str_split_fixed(string = Day_tech, pattern = "_", n = 2)
Static_data$Day <- Day_tech[,1]
Static_data$Tech <- Day_tech[,2]

p6 <- ggplot(Static_data, aes(x = Day, y = meth.wmean, fill = Tech)) +
  geom_boxplot(outliers = FALSE, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(y = "% DNAme", x = "Static - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p6
#ggsave(filename = "Static_boxplot.pdf", plot = p6, device = "pdf", path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_1/Panel-D_Spaghetti/", width = 8, height = 8, units = "cm")
