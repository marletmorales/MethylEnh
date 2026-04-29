## 22.12.2025
## Figure 5 panel E: boxplots of 5hmC trends (related to Fig 1E)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_MATAC_NonPromPeakome/")

## -- Read DNAme data
WT_D0_hMATAC <- read.table(file = "meanDNAme_WT_D0_hMATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__2025_12_22.txt", header = TRUE, sep = '\t')
WT_D2_hMATAC <- read.table(file = "meanDNAme_WT_D2_hMATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__2025_12_22.txt", header = TRUE, sep = '\t')
WT_D4_hMATAC <- read.table(file = "meanDNAme_WT_D4_hMATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__2025_12_22.txt", header = TRUE, sep = '\t')

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

## -- Read Dynamic peaks (44K)
Peaks_DA <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT_MATAC_v1/Dynamic/K5/bed/All_clusters_labeled.bed", header = FALSE, sep = '\t')
colnames(Peaks_DA) <- c("chr","start","end","cluster")

## -- Add peak and cluster information to All_data
All_data <- merge(x = All_data, y = Peaks_DA, all.x = TRUE)
All_data$cluster[is.na(All_data$cluster)] <- "static"
## Factorize
All_data$cluster <- factor(All_data$cluster, levels = c("c1", "c2", "c3", "c4", "c5", "static"))


## Reorder columns
All_data <- All_data[, c("chr", "start", "end", "coord", "cluster",
                         "num.cpgs_WT_D0_hMATAC", "meth.wmean_WT_D0_hMATAC",
                         "num.cpgs_WT_D2_hMATAC", "meth.wmean_WT_D2_hMATAC",
                         "num.cpgs_WT_D4_hMATAC", "meth.wmean_WT_D4_hMATAC")]

## Organize data for plotting
#Clust_data <- filter(All_data, cluster != "none")
Clust_data <- All_data
## Remove unnecessary columns
Clust_data <- Clust_data[,c("chr", "start", "end", "coord", "cluster",
                            "meth.wmean_WT_D0_hMATAC", "meth.wmean_WT_D2_hMATAC", "meth.wmean_WT_D4_hMATAC")]
## Pivot from wide to long
Clust_data <- pivot_longer(data = Clust_data, cols = c("meth.wmean_WT_D0_hMATAC", "meth.wmean_WT_D2_hMATAC", "meth.wmean_WT_D4_hMATAC"),
                           names_to = "Dataset", values_to = "meth.wmean", values_drop_na = TRUE)

## Split $Dataset into two: day and technique
Day_tech <- gsub(pattern = "meth.wmean_WT_", replacement = "", x = Clust_data$Dataset)
Day_tech <- str_split_fixed(string = Day_tech, pattern = "_", n = 2)
Clust_data$Day <- Day_tech[,1]
Clust_data$Tech <- Day_tech[,2]


## ---- Plot
Folder_path <- "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_5/Panel-C_Boxplots_5hmC_MATAC5clusters/"

## All clusters, by Day
p <- ggplot(Clust_data, aes(x = Day, y = meth.wmean, fill = cluster)) +
  geom_boxplot(outliers = FALSE, alpha = 0.8) +
  coord_cartesian(ylim = c(0,40)) +
  scale_fill_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67", "#95B6D6")) +
  labs(y = "% 5hmC", x = "clusters & static - Day") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p
#ggsave(filename = "All_Clusters_boxplot_5hmC_v2.pdf", plot = p, device = "pdf", path = Folder_path, width = 12, height = 8, units = "cm")


# ## ---------------------------------
# ## ---- Non-clustered peaks (static)
#
# Static_data <- filter(All_data, cluster == "none")
# ## Remove unnecessary columns
# Static_data <- Static_data[,c("chr", "start", "end", "coord", "cluster",
#                               "meth.wmean_WT_D0_hMATAC", "meth.wmean_WT_D2_hMATAC", "meth.wmean_WT_D4_hMATAC")]
#
# ## Pivot from wide to long
# Static_data <- pivot_longer(data = Static_data, cols = c("meth.wmean_WT_D0_hMATAC", "meth.wmean_WT_D2_hMATAC", "meth.wmean_WT_D4_hMATAC"),
#                             names_to = "Dataset", values_to = "meth.wmean", values_drop_na = TRUE)
#
# ## Split $Dataset into two: day and technique
# Day_tech <- gsub(pattern = "meth.wmean_WT_", replacement = "", x = Static_data$Dataset)
# Day_tech <- str_split_fixed(string = Day_tech, pattern = "_", n = 2)
# Static_data$Day <- Day_tech[,1]
# Static_data$Tech <- Day_tech[,2]
#
# ## Static
# p6 <- ggplot(Static_data, aes(x = Day, y = meth.wmean)) +
#   geom_boxplot(outliers = FALSE, alpha = 8, fill = "#95B6D6") +
#   coord_cartesian(ylim = c(0,40)) +
#   labs(y = "% 5hmC", x = "Static - Day") +
#   theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
# p6
# #ggsave(filename = "Static_boxplot.pdf", plot = p6, device = "pdf", path = Folder_path, width = 4, height = 8, units = "cm")
