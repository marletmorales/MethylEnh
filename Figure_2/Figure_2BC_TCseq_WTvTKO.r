## 25.09.2025
## Compare logFC accessibility and WT DNAme
## Analyze only WT vs TKO; export plots
## 13.11.2025: Change logFC to 1 to be congruent with the other analysis (I originally set it to 0).
## 13.11.2025: Improved writing of log file.
## 13.11.2025 (v2): Added quantification of WT_WGBS for the stats.
## 13.11.2025 (v2): Remove unnecessary plots (they are still in v1).
library(TCseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT-DKO-TKO_MATAC/")

## Specify Day of comparison
day <- "D4"
group_1 <- paste("WT_", day, sep = "")
group_2 <- paste("TKO_", day, sep = "")

Name_log <- paste("log_files/",group_1,"_vs_",group_2,"_MATAC-5clusters_stats.txt", sep = "")
cat("Comparing ",group_1," vs ",group_2,"\n",file = Name_log, sep = "", append = FALSE)

## Load environment
load(paste("environments/TCseq_",group_1,"_vs_",group_2,"_2025_09_10_MATAC-5clusters.RData", sep = ""))
## Remove unnecessary variables carried by the environment
rm(arguments, DE_peaks, experiment_design, header_eval, peaks_ref, Sub_experiment_design, Name_env, Name_file, Name_peaks, reads_path, TimeCourse)
## Specify plots folder
Plot_folder <- paste("../../../Plots/TCseq/WT-DKO-TKO_MATAC/MATAC-5clusters/WT_v_TKO/",day,"/", sep = "")

## Read results and remove unused columns
res <- as.data.frame(DBresult(tca_res, group1 = group_1, group2 = group_2,
                              p.adjust = "fdr", top.sig = FALSE, pvalue = "paj",
                              pvalue.threshold = 0.05, abs.fold = 1, direction = "both"))
res <- res[, c("seqnames", "start", "end", "width", "logFC", "PValue", "paj")]
colnames(res) <- c("chr", "start", "end", "width", "logFC", "PValue", "paj")

## Label peaks by cluster
peaks_clust <- read.table(file = "../WT_MATAC_v1/Dynamic/K5/bed/All_clusters_labeled.bed", header = FALSE, sep = '\t')
colnames(peaks_clust) <- c("chr", "start", "end", "cluster")
res <- merge(x = res, y = peaks_clust, by = c("chr", "start", "end"), all.x = TRUE, all.y = FALSE)
rm(peaks_clust)
## Sort rows and columns
res <- arrange(res, chr, start, end)
## Export number of observations
cat("n = ",nrow(res),"\n\n",file = Name_log, sep = "", append = TRUE)
cat("Peaks per cluster\n",file = Name_log, sep = "", append = TRUE)
write.table(x = data.frame(table(res$cluster)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)

## Identify differentially-accessible
res$direction <- "none"
res$direction[res$paj < 0.05 & res$logFC > 1] <- "up"
res$direction[res$paj < 0.05 & res$logFC < -1] <- "down"
res$direction <- factor(res$direction, levels = c("down", "up", "none"))

## Write number of DA peaks
cat("DA peaks\n",file = Name_log, sep = "", append = TRUE)
write.table(x = data.frame(table(res$direction)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)
cat("DA peaks per cluster\n",file = Name_log, sep = "", append = TRUE)
t <- as.data.frame(table(res$direction, res$cluster, exclude = "none"))
t <- pivot_wider(t, names_from = Var1, values_from = Freq)
cat("\tdown\tup\n",file = Name_log, sep = "\t", append = TRUE)
write.table(x = t, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)
rm(t)

## -- Before importing DNAme
## Volcano plot, labelling DA dots by cluster
res$class <- res$cluster
res$class[res$direction == "none"] <- "Not DA"
res$class <- factor(res$class, levels = c("c1", "c2", "c3", "c4", "c5", "Not DA"))
lim_fc_res <- max(abs(res$logFC))

plot <- ggplot(res, aes(x = logFC, y = -log10(paj), color = class)) +
  geom_point(size = 1.2, alpha = 0.4) +
  scale_color_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67", "gray20"), drop = FALSE) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
  geom_vline(xintercept = c(-1,1), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
  labs(x = "logFC(TKO/WT)", color = "Cluster") +
  coord_cartesian(xlim = c(-lim_fc_res, lim_fc_res)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
plot
Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DA_volcano_color_cluster.pdf", sep = "")
ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

Name_png <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DA_volcano_color_cluster.png", sep = "")
ggsave(filename = Name_png, device = "png", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

# ## Get cluster color legend
# p <- ggplot(res, aes(x = logFC, y = -log10(paj), color = class)) +
#   geom_point(size = 3, alpha = 1, shape = 15) + labs(color = "Cluster") + theme_bw() +
#   scale_color_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67", "gray20"), drop = FALSE)
# p <- cowplot::ggdraw(cowplot::get_legend(p))
# Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_legend_cluster.pdf", sep = "")
# ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = p, width = 4, height = 6, units = "cm")
# rm(p)


## ------------------ ##
## -- Import DNAme -- ##
## ------------------ ##

Name_WT <- paste("../../Measure_DNAme/DNAme_ROIs/WT_MATAC_NonPromPeakome/meanDNAme_WT_",day,"_MATAC_mergedreps__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", sep = "")
Name_KO <- paste("../../Measure_DNAme/DNAme_ROIs/WT_MATAC_NonPromPeakome/meanDNAme_TKO_",day,"_MATAC_mergedreps__Coordinates_MATAC-5clusters__v1_2025_09_10.txt", sep = "")
Name_WT_WGBS <- paste("../../Measure_DNAme/DNAme_ROIs/WT_MATAC_NonPromPeakome/meanDNAme_WT_",day,"_WGBS__Coordinates_WT_MATAC_NonPromPeakome__v1_2025_07_28.txt", sep = "")
DNAme_WT <- read.table(file = Name_WT, header = TRUE, sep = '\t')
DNAme_KO <- read.table(file = Name_KO, header = TRUE, sep = '\t')
DNAme_WT_WGBS <- read.table(file = Name_WT_WGBS, header = TRUE, sep = '\t')

## Remove unused columns
DNAme_WT <- DNAme_WT[, c("chr", "start", "end", "num.cpgs", "meth.wmean")]
DNAme_KO <- DNAme_KO[, c("chr", "start", "end", "num.cpgs", "meth.wmean")]
DNAme_WT_WGBS <- DNAme_WT_WGBS[, c("chr", "start", "end", "num.cpgs", "meth.wmean")]

## Rename columns by genotype and merge with results
colnames(DNAme_WT) <- c("chr", "start", "end", "num.cpgs.WT", "meth.wmean.WT")
colnames(DNAme_KO) <- c("chr", "start", "end", "num.cpgs.KO", "meth.wmean.KO")
colnames(DNAme_WT_WGBS) <- c("chr", "start", "end", "num.cpgs.WT_WGBS", "meth.wmean.WT_WGBS")
res <- merge(x = res, y = DNAme_WT, by = c("chr", "start", "end"), all.x = TRUE)
res <- merge(x = res, y = DNAme_KO, by = c("chr", "start", "end"), all.x = TRUE)
res <- merge(x = res, y = DNAme_WT_WGBS, by = c("chr", "start", "end"), all.x = TRUE)
rm(DNAme_WT, DNAme_KO, DNAme_WT_WGBS)


## -------------------------------------------- ##
## -- Plot logFC vs DNAme; scatter + density -- ##
## -------------------------------------------- ##
## Density plots don't work with coord_cartesian (misaligned)

## Grouping by direction (up, down, none)
plot <- ggplot(res, aes(x = logFC, y = meth.wmean.WT, color = direction)) +
  geom_vline(xintercept = 0, color = 'gray30', linewidth = 0.3, linetype = "solid") +
  geom_point(size = 0.9, alpha = 0.2, na.rm = TRUE) +
  scale_color_manual(values = c("#018571","#A6611A","gray30")) +
  #coord_cartesian(xlim = c(-lim_fc_res, lim_fc_res), ylim = c(0,100)) +
  labs(x = "logFC(TKO/WT)", y = "DNAme WT") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
plot <- ggMarginal(plot, type = 'density', margins = 'both', size = 5, groupColour = TRUE, groupFill = TRUE, alpha = 0.4)
plot
# Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_logFC_vs_DNAmeWT_wDens_all_narm_byDA.pdf", sep = "")
# ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

## Grouping by cluster (c1:c5 + no DA)
plot <- ggplot(res, aes(x = logFC, y = meth.wmean.WT, color = class)) +
  geom_vline(xintercept = 0, color = 'gray30', linewidth = 0.3, linetype = "solid") +
  geom_point(size = 1.1, alpha = 0.2, na.rm = TRUE) +
  scale_color_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67", "gray30"), drop = FALSE) +
  #scale_color_manual(values = c("#D2D365", "#85C661", "#52DED7", "#46A4DB", "#1C62AB", "gray30"), drop = FALSE) +
  #coord_cartesian(xlim = c(-lim_fc_res, lim_fc_res), ylim = c(0,100)) +
  labs(x = "logFC(TKO/WT)", y = "DNAme WT") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
plot <- ggMarginal(plot, type = 'density', margins = 'both', size = 5, groupColour = TRUE, groupFill = TRUE, alpha = 0.15)
plot
# Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_logFC_vs_DNAmeWT_wDens_all_narm_byclust.pdf", sep = "")
# ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")


## -------------------------------------- ##
## -- Plot DNAme distribution; boxplot -- ##
## -------------------------------------- ##

res_longer <- pivot_longer(data = res, cols = c(meth.wmean.WT, meth.wmean.KO, meth.wmean.WT_WGBS), names_to = "genotype", values_to = "wmean.DNAme")
res_longer$genotype <- gsub(pattern = "meth.wmean.", replacement = "", x = res_longer$genotype)
res_longer$genotype <- gsub(pattern = "KO", replacement = "TKO", x = res_longer$genotype)
res_longer$genotype <- factor(res_longer$genotype, levels = c("WT", "WT_WGBS", "TKO"))
## Drop observations in which DNAme was not measured
res_longer <- filter(res_longer, !is.na(wmean.DNAme))
table(res_longer$genotype)

## Write number of peaks with measured DNAme
cat("\nPeaks with 5mC measured\n",file = Name_log, sep = "", append = TRUE)
write.table(x = data.frame(table(res_longer$genotype)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)
cat("Peaks with 5mC measured per cluster\n",file = Name_log, sep = "", append = TRUE)
t <- data.frame(table(res_longer$cluster, res_longer$genotype))
t <- pivot_wider(t, names_from = Var2, values_from = Freq)
cat("\tWT\tWT_WGBS\tTKO\n",file = Name_log, sep = "\t", append = TRUE)
write.table(x = t, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)
rm(t)

## Boxplots
plot <- ggplot(res_longer, aes(x = cluster, y = wmean.DNAme, fill = genotype)) +
  geom_boxplot(outlier.shape = 21, outliers = FALSE, na.rm = TRUE) +
  scale_fill_manual(values = c("#266799", "#3B90AF", "#F86532")) +
  labs(y = "5mC %", x = "Cluster") +
  coord_cartesian(ylim = c(0,100)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
plot
Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DNAme_boxplots.pdf", sep = "")
ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 9, units = "cm")

# ## Get color legend
# p <- ggplot(res_longer, aes(x = cluster, y = wmean.DNAme, fill = genotype)) +
#   geom_boxplot(na.rm = TRUE) +
#   scale_fill_manual(values = c("#266799", "#3B90AF", "#F86532")) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "right")
# p <- cowplot::ggdraw(cowplot::get_legend(p))
# Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DNAme_boxplots_legend.pdf", sep = "")
# ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = p, width = 4, height = 6, units = "cm")
# rm(p)

## -- Calculate DNAme per cluster; perform wilcox test between WT and TKO
output <- data.frame("Cluster" = character(),
                     "WT_n" = numeric(),
                     "TKO_n" = numeric(),
                     "WT_5mC" = numeric(),
                     "TKO_5mC" = numeric(),
                     "pval" = numeric())
for(i in 1:5){
  Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]
  WT <- res_longer |> filter(cluster == Current_cluster & genotype == "WT") |> pull(wmean.DNAme)
  KO <- res_longer |> filter(cluster == Current_cluster & genotype == "TKO") |> pull(wmean.DNAme)
  test <- wilcox.test(x = WT, y = KO, paired = FALSE, na.rm = TRUE)
  tmp <- data.frame("Cluster" = Current_cluster,
                    "WT_n" = length(WT),
                    "TKO_n" = length(KO),
                    "WT_5mC" = mean(WT, na.rm = TRUE),
                    "TKO_5mC" = mean(KO, na.rm = TRUE),
                    "pval" = test$p.value)
  output <- rbind(output, tmp)
}
output$paj <- p.adjust(output$pval, method = "bonferroni")
## Format before printing
output$WT_5mC <- round(output$WT_5mC, digits = 2)
output$TKO_5mC <- round(output$TKO_5mC, digits = 2)
output$pval <- signif(output$pval, digits = 5)
output$paj <- signif(output$paj, digits = 5)

cat("\n5mC comparison per cluster (MATAC; all peaks)\n",file = Name_log, sep = "", append = TRUE)
suppressWarnings(write.table(x = output, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
cat("\n",file = Name_log, sep = "", append = TRUE)
rm(output, test, tmp, WT, KO)


## -------------------------------- ##
## -- WT 5mC by direction: MATAC -- ##
## -------------------------------- ##

## -- Comparing WT 5mC in down vs none vs up. Make a copy to reorder factors

res_WT <- res
## Remove NA (simplifies things downstream)
res_WT <- res_WT[!is.na(res_WT$meth.wmean.WT),]
cat("\nDA Peaks with 5mC in WT MATAC\n", file = Name_log, sep = "", append = TRUE)
write.table(x = data.frame(table(res_WT$direction)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)
## Reorder factors
res_WT$direction <- factor(res_WT$direction, levels = c("down", "none", "up"))

## -- Plot WT DNAme distributions by direction
plot <- ggplot(res_WT, aes(x = direction, y = meth.wmean.WT, fill = direction)) +
  geom_boxplot(outlier.shape = 21, outliers = TRUE) +
  scale_fill_manual(values = c("#266799","#266799","#266799")) +
  labs(y = "WT 5mC %") +
  scale_x_discrete(labels = c("Down", "No DA", "Up")) +
  coord_cartesian(ylim = c(0,100)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_blank())
plot
Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DNAmeWT_MATAC_all_byDA.pdf", sep = "")
ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 6, height = 10, units = "cm")

## -- Same, but now plot by cluster
plot <- plot + facet_wrap(~cluster, nrow = 1) + theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot
Name_plot <- sub(pattern = ".pdf", replacement = "_byclust.pdf", x = Name_plot)
ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 28, height = 5.5, units = "cm")


## -- Calculate DNAme per direction; perform wilcox test
## - All peaks
mC_down <- res_WT$meth.wmean.WT[res_WT$direction == "down"]
mC_none <- res_WT$meth.wmean.WT[res_WT$direction == "none"]
mC_up <- res_WT$meth.wmean.WT[res_WT$direction == "up"]
t_dn <- wilcox.test(x = mC_down, y = mC_none)
t_un <- wilcox.test(x = mC_up, y = mC_none)
t_du <- wilcox.test(x = mC_down, y = mC_up)

mC_output <- data.frame("A" = c("down", "up", "down"),
                        "B" = c("none", "none", "up"),
                        "A_n" = c(length(mC_down), length(mC_up), length(mC_down)),
                        "B_n" = c(length(mC_none), length(mC_none), length(mC_up)),
                        "A_5mC" = round(c(mean(mC_down), mean(mC_up), mean(mC_down)), digits = 2),
                        "B_5mC" = round(c(mean(mC_none), mean(mC_none), mean(mC_up)), digits = 2),
                        "pval" = c(t_dn$p.value, t_un$p.value, t_du$p.value))
mC_output$paj <- p.adjust(mC_output$pval, method = "bonferroni")
mC_output$pval <- signif(mC_output$pval, digits = 5)
mC_output$paj <- signif(mC_output$paj, digits = 5)

cat("DA Peaks with 5mC in WT MATAC: Comparisons\n",file = Name_log, sep = "", append = TRUE)
suppressWarnings(write.table(x = mC_output, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
cat("\n",file = Name_log, sep = "\t", append = TRUE)
rm(mC_down, mC_none, mC_up, t_dn, t_un, t_du, mC_output)

## - Peaks by cluster
mC_output_c <- data.frame("Cluster" = character(),
                          "A" = character(),
                          "B" = character(),
                          "A_n" = numeric(),
                          "B_n" = numeric(),
                          "A_5mC" = numeric(),
                          "B_5mC" = numeric(),
                          "pval" = numeric())

for(i in 1:5){
  Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]

  mC_down_c <- res_WT$meth.wmean.WT[res_WT$direction == "down" & res_WT$cluster == Current_cluster]
  mC_none_c <- res_WT$meth.wmean.WT[res_WT$direction == "none" & res_WT$cluster == Current_cluster]
  mC_up_c <- res_WT$meth.wmean.WT[res_WT$direction == "up" & res_WT$cluster == Current_cluster]

  t_dn_c <- wilcox.test(x = mC_down_c, y = mC_none_c)
  t_un_c <- wilcox.test(x = mC_up_c, y = mC_none_c)
  t_du_c <- wilcox.test(x = mC_down_c, y = mC_up_c)

  tmp <- data.frame("Cluster" = rep(Current_cluster, 3),
                    "A" = c("down", "up", "down"),
                    "B" = c("none", "none", "up"),
                    "A_n" = c(length(mC_down_c), length(mC_up_c), length(mC_down_c)),
                    "B_n" = c(length(mC_none_c), length(mC_none_c), length(mC_up_c)),
                    "A_5mC" = round(c(mean(mC_down_c), mean(mC_up_c), mean(mC_down_c)), digits = 2),
                    "B_5mC" = round(c(mean(mC_none_c), mean(mC_none_c), mean(mC_up_c)), digits = 2),
                    "pval" = c(t_dn_c$p.value, t_un_c$p.value, t_du_c$p.value))
  mC_output_c <- rbind(mC_output_c, tmp)
}
mC_output_c$paj <- p.adjust(mC_output_c$pval, method = "bonferroni")
mC_output_c$pval <- signif(mC_output_c$pval, digits = 5)
mC_output_c$paj <- signif(mC_output_c$paj, digits = 5)

cat("DA Peaks with 5mC in WT MATAC: Comparisons by cluster\n",file = Name_log, sep = "", append = TRUE)
suppressWarnings(write.table(x = mC_output_c[mC_output_c$Cluster == "c1",], file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
cat("\n",file = Name_log, sep = "\t", append = TRUE)
for(i in 2:5){
  Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]
  write.table(x = mC_output_c[mC_output_c$Cluster == Current_cluster,], file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
  cat("\n",file = Name_log, sep = "\t", append = TRUE)
}
rm(mC_down_c, mC_none_c, mC_up_c, t_dn_c, t_un_c, t_du_c, tmp, mC_output_c)


## -------------------------------- ##
## -- WT 5mC by direction: WGBS -- ##
## -------------------------------- ##

## -- Comparing WT 5mC in down vs none vs up. Make a copy to reorder factors

res_WT_WGBS <- res
## Remove NA (simplifies things downstream)
res_WT_WGBS <- res_WT_WGBS[!is.na(res_WT_WGBS$meth.wmean.WT_WGBS),]
cat("\nDA Peaks with 5mC in WT WGBS\n", file = Name_log, sep = "", append = TRUE)
write.table(x = data.frame(table(res_WT_WGBS$direction)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)
## Reorder factors
res_WT_WGBS$direction <- factor(res_WT_WGBS$direction, levels = c("down", "none", "up"))

## -- Plot WT_WGBS DNAme distributions by direction
plot <- ggplot(res_WT_WGBS, aes(x = direction, y = meth.wmean.WT_WGBS, fill = direction)) +
  geom_boxplot(outlier.shape = 21, outliers = TRUE) +
  scale_fill_manual(values = c("#3B90AF","#3B90AF","#3B90AF")) +
  labs(y = "WT_WGBS 5mC %") +
  scale_x_discrete(labels = c("Down", "No DA", "Up")) +
  coord_cartesian(ylim = c(0,100)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_blank())
plot
Name_plot <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DNAmeWT_WGBS_all_byDA.pdf", sep = "")
ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 6, height = 10, units = "cm")

## -- Same, but now plot by cluster
plot <- plot + facet_wrap(~cluster, nrow = 1) + theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot
Name_plot <- sub(pattern = ".pdf", replacement = "_byclust.pdf", x = Name_plot)
ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 28, height = 5.5, units = "cm")


## -- Calculate DNAme per direction; perform wilcox test
## - All peaks
mC_down <- res_WT_WGBS$meth.wmean.WT_WGBS[res_WT_WGBS$direction == "down"]
mC_none <- res_WT_WGBS$meth.wmean.WT_WGBS[res_WT_WGBS$direction == "none"]
mC_up <- res_WT_WGBS$meth.wmean.WT_WGBS[res_WT_WGBS$direction == "up"]
t_dn <- wilcox.test(x = mC_down, y = mC_none)
t_un <- wilcox.test(x = mC_up, y = mC_none)
t_du <- wilcox.test(x = mC_down, y = mC_up)

mC_output <- data.frame("A" = c("down", "up", "down"),
                        "B" = c("none", "none", "up"),
                        "A_n" = c(length(mC_down), length(mC_up), length(mC_down)),
                        "B_n" = c(length(mC_none), length(mC_none), length(mC_up)),
                        "A_5mC" = round(c(mean(mC_down), mean(mC_up), mean(mC_down)), digits = 2),
                        "B_5mC" = round(c(mean(mC_none), mean(mC_none), mean(mC_up)), digits = 2),
                        "pval" = c(t_dn$p.value, t_un$p.value, t_du$p.value))
mC_output$paj <- p.adjust(mC_output$pval, method = "bonferroni")
mC_output$pval <- signif(mC_output$pval, digits = 5)
mC_output$paj <- signif(mC_output$paj, digits = 5)

cat("DA Peaks with 5mC in WT WGBS: Comparisons\n",file = Name_log, sep = "", append = TRUE)
suppressWarnings(write.table(x = mC_output, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
cat("\n",file = Name_log, sep = "\t", append = TRUE)
rm(mC_down, mC_none, mC_up, t_dn, t_un, t_du, mC_output)

## - Peaks by cluster
mC_output_c <- data.frame("Cluster" = character(),
                          "A" = character(),
                          "B" = character(),
                          "A_n" = numeric(),
                          "B_n" = numeric(),
                          "A_5mC" = numeric(),
                          "B_5mC" = numeric(),
                          "pval" = numeric())

for(i in 1:5){
  Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]

  mC_down_c <- res_WT_WGBS$meth.wmean.WT_WGBS[res_WT_WGBS$direction == "down" & res_WT_WGBS$cluster == Current_cluster]
  mC_none_c <- res_WT_WGBS$meth.wmean.WT_WGBS[res_WT_WGBS$direction == "none" & res_WT_WGBS$cluster == Current_cluster]
  mC_up_c <- res_WT_WGBS$meth.wmean.WT_WGBS[res_WT_WGBS$direction == "up" & res_WT_WGBS$cluster == Current_cluster]

  t_dn_c <- wilcox.test(x = mC_down_c, y = mC_none_c)
  t_un_c <- wilcox.test(x = mC_up_c, y = mC_none_c)
  t_du_c <- wilcox.test(x = mC_down_c, y = mC_up_c)

  tmp <- data.frame("Cluster" = rep(Current_cluster, 3),
                    "A" = c("down", "up", "down"),
                    "B" = c("none", "none", "up"),
                    "A_n" = c(length(mC_down_c), length(mC_up_c), length(mC_down_c)),
                    "B_n" = c(length(mC_none_c), length(mC_none_c), length(mC_up_c)),
                    "A_5mC" = round(c(mean(mC_down_c), mean(mC_up_c), mean(mC_down_c)), digits = 2),
                    "B_5mC" = round(c(mean(mC_none_c), mean(mC_none_c), mean(mC_up_c)), digits = 2),
                    "pval" = c(t_dn_c$p.value, t_un_c$p.value, t_du_c$p.value))
  mC_output_c <- rbind(mC_output_c, tmp)
}
mC_output_c$paj <- p.adjust(mC_output_c$pval, method = "bonferroni")
mC_output_c$pval <- signif(mC_output_c$pval, digits = 5)
mC_output_c$paj <- signif(mC_output_c$paj, digits = 5)

cat("DA Peaks with 5mC in WT WGBS: Comparisons by cluster\n",file = Name_log, sep = "", append = TRUE)
suppressWarnings(write.table(x = mC_output_c[mC_output_c$Cluster == "c1",], file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
cat("\n",file = Name_log, sep = "\t", append = TRUE)
for(i in 2:5){
  Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]
  write.table(x = mC_output_c[mC_output_c$Cluster == Current_cluster,], file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
  cat("\n",file = Name_log, sep = "\t", append = TRUE)
}
rm(mC_down_c, mC_none_c, mC_up_c, t_dn_c, t_un_c, t_du_c, tmp, mC_output_c)
