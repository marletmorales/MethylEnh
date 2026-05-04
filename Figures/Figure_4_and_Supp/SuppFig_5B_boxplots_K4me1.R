## K4me1 signal in WT enhancers
## Supplemental Figure 5B
library(ggplot2)
library(dplyr)
library(stringr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/")

## ------------------- ##
## -- Data analysis -- ##
## ------------------- ##

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
D0_libsize <- sum(D0_counts$Counts)
D2_libsize <- sum(D2_counts$Counts)
D4_libsize <- sum(D4_counts$Counts)
D0_counts$Counts <- D0_counts$Counts/D0_libsize * 1000000
D2_counts$Counts <- D2_counts$Counts/D2_libsize * 1000000
D4_counts$Counts <- D4_counts$Counts/D4_libsize * 1000000

## Remove additional htseq rows and clean table
D0_counts <- D0_counts |> filter(grepl(pattern = "peak_", x = ID)) |> select(coord, Counts)
D2_counts <- D2_counts |> filter(grepl(pattern = "peak_", x = ID)) |> select(coord, Counts)
D4_counts <- D4_counts |> filter(grepl(pattern = "peak_", x = ID)) |> select(coord, Counts)

## -- Read enhancers classification
Enhancers <- read.table(file = "WT_Master_Enhancer_Peakome_v6_single-class.bed", header = FALSE, sep = '\t')
colnames(Enhancers) <- c("chr", "start", "end", "D0_class", "D2_class", "D4_class")
Enhancers$coord <- str_c(str_c(Enhancers$chr, Enhancers$start, sep = ":"), Enhancers$end, sep = "-")

## -- Append enhancer information to counts
## Merge tables
D0_counts <- merge(x = D0_counts, y = Enhancers[,c("chr", "start", "end", "coord", "D0_class")], all.x = TRUE, all.y = FALSE)
D2_counts <- merge(x = D2_counts, y = Enhancers[,c("chr", "start", "end", "coord", "D2_class")], all.x = TRUE, all.y = FALSE)
D4_counts <- merge(x = D4_counts, y = Enhancers[,c("chr", "start", "end", "coord", "D4_class")], all.x = TRUE, all.y = FALSE)
rm(Enhancers)

## Re-name columns
D0_counts <- D0_counts |> rename(class = D0_class)
D2_counts <- D2_counts |> rename(class = D2_class)
D4_counts <- D4_counts |> rename(class = D4_class)

## Add Day
D0_counts$Day <- "D0"
D2_counts$Day <- "D2"
D4_counts$Day <- "D4"

## -- Merge all tables into a long format
Counts_l <- rbind(D0_counts, D2_counts, D4_counts)
## Clean table
Counts_l <- Counts_l[,c("chr", "start", "end", "Day", "class", "Counts", "coord")]
Counts_l <- Counts_l |> arrange(chr, start, end)
## Factorize
Counts_l$Day <- factor(Counts_l$Day, levels = c("D0", "D2", "D4"))
Counts_l$class <- factor(Counts_l$class, levels = c("actv", "prim", "inac", "othr"))
## Normalized counts by enhancer width (divide by num_bp)
Counts_l$Counts_bp <- Counts_l$Counts/(Counts_l$end - Counts_l$start + 1)

rm(D0_counts, D2_counts, D4_counts)

## -------------- ##
## -- Plotting -- ##
## -------------- ##

## Plot normalized K4me1 counts (CPM/bp)
p <- ggplot(Counts_l, aes(x = class, y = Counts_bp, fill = class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 1, outliers = FALSE) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbfe6")) +
  labs(x = "Enhancer class", y = "Counts (CPM/bp)") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~Day)
p
# ggsave(filename = "WT_D0-D2-D4_K4me1_htseqcounts_CPMbp_per_enhancer_class.pdf", plot = p, device = "pdf", width = 14, height = 6, units = "cm")

## ----------- ##
## -- Stats -- ##
## ----------- ##

## -- Wilcox testing comparing 'othr' to the rest
Summary_stats <- Counts_l |> group_by(Day, class) |> summarize("Num_enh" = n(),
                                                               "Mean_norm_counts" = round(mean(Counts_bp), digits = 6),
                                                               "Median_norm_counts" = round(median(Counts_bp), digits = 6))
## D0
C_actv <- Counts_l |> filter(Day == "D0" & class == "actv") |> pull(Counts_bp)
C_prim <- Counts_l |> filter(Day == "D0" & class == "prim") |> pull(Counts_bp)
C_inac <- Counts_l |> filter(Day == "D0" & class == "inac") |> pull(Counts_bp)
C_othr <- Counts_l |> filter(Day == "D0" & class == "othr") |> pull(Counts_bp)
t_actv <- wilcox.test(x = C_actv, y = C_othr)
t_prim <- wilcox.test(x = C_prim, y = C_othr)
t_inac <- wilcox.test(x = C_inac, y = C_othr)
T_results <- data.frame("Day" = "D0",
                        "class" = c("actv", "prim", "inac", "othr"),
                        "paj" = p.adjust(p = c(t_actv$p.value, t_prim$p.value, t_inac$p.value, NA), method = "bonferroni"))

## D2
C_actv <- Counts_l |> filter(Day == "D2" & class == "actv") |> pull(Counts_bp)
C_prim <- Counts_l |> filter(Day == "D2" & class == "prim") |> pull(Counts_bp)
C_inac <- Counts_l |> filter(Day == "D2" & class == "inac") |> pull(Counts_bp)
C_othr <- Counts_l |> filter(Day == "D2" & class == "othr") |> pull(Counts_bp)
t_actv <- wilcox.test(x = C_actv, y = C_othr)
t_prim <- wilcox.test(x = C_prim, y = C_othr)
t_inac <- wilcox.test(x = C_inac, y = C_othr)
T_results <- rbind(T_results, data.frame("Day" = "D2",
                                         "class" = c("actv", "prim", "inac", "othr"),
                                         "paj" = p.adjust(p = c(t_actv$p.value, t_prim$p.value, t_inac$p.value, NA), method = "bonferroni")))

## D4
C_actv <- Counts_l |> filter(Day == "D4" & class == "actv") |> pull(Counts_bp)
C_prim <- Counts_l |> filter(Day == "D4" & class == "prim") |> pull(Counts_bp)
C_inac <- Counts_l |> filter(Day == "D4" & class == "inac") |> pull(Counts_bp)
C_othr <- Counts_l |> filter(Day == "D4" & class == "othr") |> pull(Counts_bp)
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
#             file = "Supp-B_Boxplot_K4me1_read_counts_per enhancer_type/Info.txt")

