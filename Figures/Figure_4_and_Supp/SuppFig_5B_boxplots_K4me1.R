## K4me1 signal in WT enhancers
## Supplemental Figure 5B
## Edit 2026.11.23: Changed input file (v2 of enhancer classification). Changed color coding. Make code cleaner (actually I ended up changing the entire code).
library(ggplot2)
library(dplyr)
library(stringr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/")

## --------------------- ##
## -- Data processing -- ##
## --------------------- ##

## -- Read K4me1 signal
## Import counts
D0_counts <- read.table(file = "../../../htseq_count/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_D0_K4me1_MCnT_mergedreps_K4me1_counts.txt", header = FALSE, sep = '\t')
D2_counts <- read.table(file = "../../../htseq_count/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_D2_K4me1_MCnT_mergedreps_K4me1_counts.txt", header = FALSE, sep = '\t')
D4_counts <- read.table(file = "../../../htseq_count/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_D4_K4me1_MCnT_mergedreps_K4me1_counts.txt", header = FALSE, sep = '\t')

## Rename columns
colnames(D0_counts) <- c("ID", "coord", "Counts")
colnames(D2_counts) <- c("ID", "coord", "Counts")
colnames(D4_counts) <- c("ID", "coord", "Counts")

## Calculate and normalize by library size
D0_counts$Counts <- D0_counts$Counts/sum(D0_counts$Counts) * 1000000
D2_counts$Counts <- D2_counts$Counts/sum(D2_counts$Counts) * 1000000
D4_counts$Counts <- D4_counts$Counts/sum(D4_counts$Counts) * 1000000

## Combine all tables into one
D0_counts$Day <- "D0"
D2_counts$Day <- "D2"
D4_counts$Day <- "D4"
All_counts <- rbind(D0_counts, D2_counts, D4_counts)
## Remove additional htseq rows and clean table
All_counts <- All_counts |> filter(grepl(pattern = "chr", x = coord)) |> select(-ID)
rm(D0_counts, D2_counts, D4_counts)

## -- Read enhancers classification
Enhancers <- read.table(file = "WT_Master_Enhancer_Peakome_v6_single-class_wActvM_v2.txt", header = TRUE, sep = '\t')
Enhancers$coord <- str_c(str_c(Enhancers$chr, Enhancers$start, sep = ":"), Enhancers$end, sep = "-")
Enhancers <- Enhancers |> select(chr, start, end, coord, D0_class, D2_class, D4_class)
## Pivot longer
Enhancers <- Enhancers |> tidyr::pivot_longer(cols = c(D0_class, D2_class, D4_class), names_to = "Day", values_to = "class")
Enhancers$Day <- sub(pattern = "_class", replacement = "", x = Enhancers$Day)

## -- Append enhancer information to counts
Enhancers <- inner_join(x = Enhancers, y = All_counts, by = c("coord", "Day"))
## Normalized counts by enhancer width
Enhancers$Counts_bp <- Enhancers$Counts/(Enhancers$end - Enhancers$start + 1)
## Factorize
Enhancers$Day <- factor(Enhancers$Day, levels = c("D0", "D2", "D4"))
Enhancers$class <- factor(Enhancers$class, levels = c("actv", "prim", "inac", "othr"))

rm(All_counts)

## -------------- ##
## -- Plotting -- ##
## -------------- ##

## Normalized K4me1 counts (CPM/bp)
Plot_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Paper_additional_files/Figures_Revision_1/Figure_4/Supp_5B/"

Color_scale <- c("actv" = "#b1c848",
                 "actvM" = "#2cb11b",
                 "prim" = "#1bb6af",
                 "inac" = "#485d7b",
                 "othr" = "#d7d7d7")

ylim <- Enhancers |> group_by(Day, class) |>
  summarise("IQR" = 1.5*IQR(Counts_bp),
            "Q75" = quantile(Counts_bp, 0.75),
            "TopWhisker" = IQR + Q75)
ylim <- c(0, max(ylim$TopWhisker))

## -- D0
p1 <- ggplot(Enhancers[Enhancers$Day == "D0",], aes(x = class, y = Counts_bp, fill = class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 1, outliers = FALSE) +
  scale_fill_manual(values = Color_scale) +
  labs(x = "Enhancer class", y = "Counts (CPM/bp)") +
  coord_cartesian(ylim = ylim) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank(), axis.title.x = element_blank())
p1
#ggsave(filename = "WT_D0_K4me1_htseqcounts_CPMbp_per_enhancer_class.pdf", plot = p1, width = 9, height = 9, units = "cm", path = Plot_folder)

## -- D2
p2 <- ggplot(Enhancers[Enhancers$Day == "D2",], aes(x = class, y = Counts_bp, fill = class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 1, outliers = FALSE) +
  scale_fill_manual(values = Color_scale) +
  labs(x = "Enhancer class", y = "Counts (CPM/bp)") +
  coord_cartesian(ylim = ylim) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank(), axis.title.x = element_blank())
p2
#ggsave(filename = "WT_D2_K4me1_htseqcounts_CPMbp_per_enhancer_class.pdf", plot = p2, width = 9, height = 9, units = "cm", path = Plot_folder)

## -- D4
p3 <- ggplot(Enhancers[Enhancers$Day == "D4",], aes(x = class, y = Counts_bp, fill = class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 1, outliers = FALSE) +
  scale_fill_manual(values = Color_scale) +
  labs(x = "Enhancer class", y = "Counts (CPM/bp)") +
  coord_cartesian(ylim = ylim) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank(), axis.title.x = element_blank())
p3
#ggsave(filename = "WT_D4_K4me1_htseqcounts_CPMbp_per_enhancer_class.pdf", plot = p3, width = 9, height = 9, units = "cm", path = Plot_folder)


## ----------- ##
## -- Stats -- ##
## ----------- ##

## -- Wilcox testing comparing 'othr' to the rest
Summary_stats <- Enhancers |> group_by(Day, class) |> summarize("Num_enh" = n(),
                                                                "Mean_norm_counts" = round(mean(Counts_bp), digits = 6),
                                                                "Median_norm_counts" = round(median(Counts_bp), digits = 6))
## D0
C_actv <- Enhancers |> filter(Day == "D0" & class == "actv") |> pull(Counts_bp)
C_prim <- Enhancers |> filter(Day == "D0" & class == "prim") |> pull(Counts_bp)
C_inac <- Enhancers |> filter(Day == "D0" & class == "inac") |> pull(Counts_bp)
C_othr <- Enhancers |> filter(Day == "D0" & class == "othr") |> pull(Counts_bp)
t_actv <- wilcox.test(x = C_actv, y = C_othr)
t_prim <- wilcox.test(x = C_prim, y = C_othr)
t_inac <- wilcox.test(x = C_inac, y = C_othr)
T_results <- data.frame("Day" = "D0",
                        "class" = c("actv", "prim", "inac", "othr"),
                        "paj" = p.adjust(p = c(t_actv$p.value, t_prim$p.value, t_inac$p.value, NA), method = "bonferroni"))

## D2
C_actv <- Enhancers |> filter(Day == "D2" & class == "actv") |> pull(Counts_bp)
C_prim <- Enhancers |> filter(Day == "D2" & class == "prim") |> pull(Counts_bp)
C_inac <- Enhancers |> filter(Day == "D2" & class == "inac") |> pull(Counts_bp)
C_othr <- Enhancers |> filter(Day == "D2" & class == "othr") |> pull(Counts_bp)
t_actv <- wilcox.test(x = C_actv, y = C_othr)
t_prim <- wilcox.test(x = C_prim, y = C_othr)
t_inac <- wilcox.test(x = C_inac, y = C_othr)
T_results <- rbind(T_results, data.frame("Day" = "D2",
                                         "class" = c("actv", "prim", "inac", "othr"),
                                         "paj" = p.adjust(p = c(t_actv$p.value, t_prim$p.value, t_inac$p.value, NA), method = "bonferroni")))

## D4
C_actv <- Enhancers |> filter(Day == "D4" & class == "actv") |> pull(Counts_bp)
C_prim <- Enhancers |> filter(Day == "D4" & class == "prim") |> pull(Counts_bp)
C_inac <- Enhancers |> filter(Day == "D4" & class == "inac") |> pull(Counts_bp)
C_othr <- Enhancers |> filter(Day == "D4" & class == "othr") |> pull(Counts_bp)
t_actv <- wilcox.test(x = C_actv, y = C_othr)
t_prim <- wilcox.test(x = C_prim, y = C_othr)
t_inac <- wilcox.test(x = C_inac, y = C_othr)
T_results <- rbind(T_results, data.frame("Day" = "D4",
                                         "class" = c("actv", "prim", "inac", "othr"),
                                         "paj" = p.adjust(p = c(t_actv$p.value, t_prim$p.value, t_inac$p.value, NA), method = "bonferroni")))

## Append results to summary stats
Summary_stats <- merge(Summary_stats, T_results)
## Sort
Summary_stats <- Summary_stats |> arrange(Day, class)
## Export
# write.table(x = Summary_stats, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE,
#             file = paste0(Plot_folder, "/Info.txt"))

