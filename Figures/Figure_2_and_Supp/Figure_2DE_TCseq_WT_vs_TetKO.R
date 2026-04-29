## Plot DA regions between WT and Tet KO
## Figure 2D,E and Supplemental Figure 2F
library(TCseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT-DKO-TKO_MATAC/")

## Specify Day of comparison
day <- "D4"
group_1 <- paste("WT_", day, sep = "")
group_2 <- paste("DKO_", day, sep = "")

# Name_log <- paste("log_files/",group_1,"_vs_",group_2,"_MATAC-5clusters_stats.txt", sep = "")
# cat("Comparing ",group_1," vs ",group_2,"\n",file = Name_log, sep = "", append = FALSE)

## Load environment
if(day == "D2") {
  load("environments/TCseq_WT_D2_vs_DKO_D2_2025_09_10_MATAC-5clusters_TEST-FAIL.RData")
} else {
  load(paste("environments/TCseq_",group_1,"_vs_",group_2,"_2025_09_10_MATAC-5clusters.RData", sep = ""))
}
## Remove unnecessary variables carried by the environment
rm(arguments, experiment_design, header_eval, peaks_ref, Sub_experiment_design, Name_env, Name_peaks, reads_path, TimeCourse)
if(day != "D2") {rm(DE_peaks, Name_file)}

## Specify plots folder
Plot_folder <- paste("../../../Plots/TCseq/WT-DKO-TKO_MATAC/MATAC-5clusters/WT_v_DKO/",day,"/", sep = "")

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

if(sum(res$direction == "none") != nrow(res)) {
  cat("DA peaks per cluster\n",file = Name_log, sep = "", append = TRUE)
  t <- as.data.frame(table(res$direction, res$cluster, exclude = "none"))
  t <- pivot_wider(t, names_from = Var1, values_from = Freq)
  cat("\tdown\tup\n",file = Name_log, sep = "\t", append = TRUE)
  write.table(x = t, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
  cat("\n",file = Name_log, sep = "\t", append = TRUE)
  rm(t)
}

## -- Before importing DNAme
## Volcano plot, labelling DA dots by cluster
res$class <- res$cluster
res$class[res$direction == "none"] <- "Not DA"
res$class <- factor(res$class, levels = c("c1", "c2", "c3", "c4", "c5", "Not DA"))
res$cluster <- factor(res$cluster, levels = c("c1", "c2", "c3", "c4", "c5"))
lim_fc_res <- max(abs(res$logFC))

plot <- ggplot(res, aes(x = logFC, y = -log10(paj), color = class)) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67", "gray20"), drop = FALSE) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
  geom_vline(xintercept = c(-1,1), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
  labs(x = "logFC(DKO/WT)", color = "Cluster") +
  coord_cartesian(xlim = c(-lim_fc_res, lim_fc_res)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
plot
Name_plot <- paste("WT_",day,"_vs_DKO_",day,"_MATAC-5clusters_DA_volcano_color_cluster.pdf", sep = "")
#ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

Name_png <- paste("WT_",day,"_vs_DKO_",day,"_MATAC-5clusters_DA_volcano_color_cluster.png", sep = "")
#ggsave(filename = Name_png, device = "png", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

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
Name_KO <- paste("../../Measure_DNAme/DNAme_ROIs/WT_MATAC_NonPromPeakome/meanDNAme_DKO_",day,"_MATAC_mergedreps__Coordinates_MATAC-5clusters__v1_2025_09_10.txt", sep = "")
DNAme_WT <- read.table(file = Name_WT, header = TRUE, sep = '\t')
DNAme_KO <- read.table(file = Name_KO, header = TRUE, sep = '\t')

## Remove unused columns
DNAme_WT <- DNAme_WT[, c("chr", "start", "end", "num.cpgs", "meth.wmean")]
DNAme_KO <- DNAme_KO[, c("chr", "start", "end", "num.cpgs", "meth.wmean")]

## Rename columns by genotype and merge with results
colnames(DNAme_WT) <- c("chr", "start", "end", "num.cpgs.WT", "meth.wmean.WT")
colnames(DNAme_KO) <- c("chr", "start", "end", "num.cpgs.KO", "meth.wmean.KO")
res <- merge(x = res, y = DNAme_WT, by = c("chr", "start", "end"), all.x = TRUE)
res <- merge(x = res, y = DNAme_KO, by = c("chr", "start", "end"), all.x = TRUE)
rm(DNAme_WT, DNAme_KO)

res_both <- filter(res, !is.na(num.cpgs.WT) & !is.na(num.cpgs.KO))

## -------------------------------------------- ##
## -- Plot logFC vs DNAme; scatter + density -- ##
## -------------------------------------------- ##
## Density plots don't work with coord_cartesian (misaligned)

## Grouping by cluster (c1:c5)
plot <- ggplot(res_both, aes(x = logFC, y = meth.wmean.KO - meth.wmean.WT, color = cluster)) +
  geom_vline(xintercept = 0, color = 'gray30', linewidth = 0.3, linetype = "solid") +
  geom_hline(yintercept = 0, color = 'gray30', linewidth = 0.3, linetype = "solid") +
  geom_point(size = 1.1, alpha = 0.2, na.rm = TRUE) +
  scale_color_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67"), drop = FALSE) +
  #coord_cartesian(xlim = c(-lim_fc_res, lim_fc_res), ylim = c(0,100)) +
  labs(x = "logFC(DKO/WT)", y = "d5mC (DKO - WT)") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
plot
plot <- ggMarginal(plot, type = 'density', margins = 'both', size = 5, groupColour = TRUE, groupFill = TRUE, alpha = 0.15)
plot
Name_png <- paste("WT_",day,"_vs_DKO_",day,"_MATAC-5clusters_logFC_vs_deltaDNAme_wDens_both_byclust.png", sep = "")
#ggsave(filename = Name_png, device = "png", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

## -------------------------------------- ##
## -- Plot DNAme distribution; boxplot -- ##
## -------------------------------------- ##

# ## -------------------------
# ## -- All peaks measured
#
# res_longer <- pivot_longer(data = res, cols = c(meth.wmean.WT, meth.wmean.KO), names_to = "genotype", values_to = "wmean.DNAme")
# res_longer$genotype <- gsub(pattern = "meth.wmean.", replacement = "", x = res_longer$genotype)
# res_longer$genotype <- gsub(pattern = "KO", replacement = "DKO", x = res_longer$genotype)
# res_longer$genotype <- factor(res_longer$genotype, levels = c("WT", "DKO"))
# ## Drop observations in which DNAme was not measured
# res_longer <- filter(res_longer, !is.na(wmean.DNAme))
# table(res_longer$genotype)
#
# ## Write number of peaks with measured DNAme
# cat("\nPeaks with 5mC measured\n",file = Name_log, sep = "", append = TRUE)
# write.table(x = data.frame(table(res_longer$genotype)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
# cat("\n",file = Name_log, sep = "\t", append = TRUE)
# cat("Peaks with 5mC measured per cluster\n",file = Name_log, sep = "", append = TRUE)
# t <- data.frame(table(res_longer$cluster, res_longer$genotype))
# t <- pivot_wider(t, names_from = Var2, values_from = Freq)
# cat("\tWT\tDKO\n",file = Name_log, sep = "\t", append = TRUE)
# write.table(x = t, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
# cat("\n",file = Name_log, sep = "\t", append = TRUE)
# rm(t)
#
# ## Boxplots
# plot <- ggplot(res_longer, aes(x = cluster, y = wmean.DNAme, fill = genotype)) +
#   geom_boxplot(outlier.shape = 21, outliers = FALSE, na.rm = TRUE) +
#   scale_fill_manual(values = c("#266799", "#4c4c4c")) +
#   labs(y = "5mC %", x = "Cluster") +
#   coord_cartesian(ylim = c(0,100)) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
# plot
#
# ## -- Calculate DNAme per cluster; perform wilcox test between WT and DKO
# output <- data.frame("Cluster" = character(),
#                      "WT_n" = numeric(),
#                      "DKO_n" = numeric(),
#                      "WT_5mC" = numeric(),
#                      "DKO_5mC" = numeric(),
#                      "pval" = numeric())
# for(i in 1:5){
#   Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]
#   WT <- res_longer |> filter(cluster == Current_cluster & genotype == "WT") |> pull(wmean.DNAme)
#   KO <- res_longer |> filter(cluster == Current_cluster & genotype == "DKO") |> pull(wmean.DNAme)
#   w_test <- wilcox.test(x = WT, y = KO, paired = FALSE, na.rm = TRUE)
#   tmp <- data.frame("Cluster" = Current_cluster,
#                     "WT_n" = length(WT),
#                     "DKO_n" = length(KO),
#                     "WT_5mC" = round(mean(WT, na.rm = TRUE), digits = 2),
#                     "DKO_5mC" = round(mean(KO, na.rm = TRUE), digits = 2),
#                     "pval" = w_test$p.value)
#   output <- rbind(output, tmp)
# }
# output$paj <- p.adjust(output$pval, method = "bonferroni")
# ## Format before printing
# output$pval <- signif(output$pval, digits = 5)
# output$paj <- signif(output$paj, digits = 5)
#
# cat("5mC comparison per cluster (all peaks)\n",file = Name_log, sep = "", append = TRUE)
# suppressWarnings(write.table(x = output, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
# cat("\n",file = Name_log, sep = "", append = TRUE)
# rm(output, w_test, tmp, WT, KO)
# ## Repeat but for peaks in which both samples have DNAme measurement

## -----------------------------------
## -- Peaks with both measurements

cat("\nPeaks with 5mC measured in both samples\n",file = Name_log, sep = "", append = TRUE)
cat("n = ",nrow(res_both),"\n",file = Name_log, sep = "", append = TRUE)
write.table(x = data.frame(table(res_both$cluster)), file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = FALSE)
cat("\n",file = Name_log, sep = "\t", append = TRUE)

res_both_longer <- pivot_longer(data = res_both, cols = c(meth.wmean.WT, meth.wmean.KO), names_to = "genotype", values_to = "wmean.DNAme")
res_both_longer$genotype <- gsub(pattern = "meth.wmean.", replacement = "", x = res_both_longer$genotype)
res_both_longer$genotype <- gsub(pattern = "KO", replacement = "DKO", x = res_both_longer$genotype)
res_both_longer$genotype <- factor(res_both_longer$genotype, levels = c("WT", "DKO"))

## Boxplots
plot <- ggplot(res_both_longer, aes(x = cluster, y = wmean.DNAme, fill = genotype)) +
  geom_boxplot(outlier.shape = 21, outliers = FALSE, na.rm = TRUE) +
  scale_fill_manual(values = c("#266799", "#4c4c4c")) +
  labs(y = "5mC %", x = "Cluster") +
  coord_cartesian(ylim = c(0,100)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
plot
Name_plot <- paste("WT_",day,"_vs_DKO_",day,"_MATAC-5clusters_DNAme_boxplots_both.pdf", sep = "")
#ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 9, units = "cm")

## -- Calculate DNAme per cluster; perform wilcox test between WT and DKO
output_both <- data.frame("Cluster" = character(),
                          "WT_n" = numeric(),
                          "DKO_n" = numeric(),
                          "WT_5mC" = numeric(),
                          "DKO_5mC" = numeric(),
                          "pval" = numeric())
for(i in 1:5){
  Current_cluster <- c("c1", "c2", "c3", "c4", "c5")[i]
  WT <- res_both_longer |> filter(cluster == Current_cluster & genotype == "WT") |> pull(wmean.DNAme)
  KO <- res_both_longer |> filter(cluster == Current_cluster & genotype == "DKO") |> pull(wmean.DNAme)
  w_test <- wilcox.test(x = WT, y = KO, paired = FALSE, na.rm = TRUE)
  tmp <- data.frame("Cluster" = Current_cluster,
                    "WT_n" = length(WT),
                    "DKO_n" = length(KO),
                    "WT_5mC" = round(mean(WT, na.rm = TRUE), digits = 2),
                    "DKO_5mC" = round(mean(KO, na.rm = TRUE), digits = 2),
                    "pval" = w_test$p.value)
  output_both <- rbind(output_both, tmp)
}
output_both$paj <- p.adjust(output_both$pval, method = "bonferroni")
## Format before printing
output_both$pval <- signif(output_both$pval, digits = 5)
output_both$paj <- signif(output_both$paj, digits = 5)

cat("5mC comparison per cluster (both samples)\n",file = Name_log, sep = "", append = TRUE)
suppressWarnings(write.table(x = output_both, file = Name_log, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE))
cat("\n",file = Name_log, sep = "", append = TRUE)
rm(output_both, w_test, tmp, WT, KO)
