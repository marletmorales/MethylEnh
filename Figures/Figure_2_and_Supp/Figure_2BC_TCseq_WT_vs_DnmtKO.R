## Plot DA regions between WT and Dnmt KO
## Figure 2B,C and Supplemental Figure 2B,C,E
library(TCseq)
library(ggplot2)
library(dplyr)
library(tidyr)
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
#ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 12, height = 12, units = "cm")

Name_png <- paste("WT_",day,"_vs_TKO_",day,"_MATAC-5clusters_DA_volcano_color_cluster.png", sep = "")
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
DNAme_WT <- read.table(file = Name_WT, header = TRUE, sep = '\t')
## Remove unused columns
DNAme_WT <- DNAme_WT[, c("chr", "start", "end", "num.cpgs", "meth.wmean")]
## Rename columns by genotype and merge with results
colnames(DNAme_WT) <- c("chr", "start", "end", "num.cpgs.WT", "meth.wmean.WT")
res <- merge(x = res, y = DNAme_WT, by = c("chr", "start", "end"), all.x = TRUE)
rm(DNAme_WT)


## ------------------------------ ##
## -- Plot WT 5mC by direction -- ##
## ------------------------------ ##

## -- Comparing WT 5mC in down vs none vs up.

## Make a copy to remove NA (simplifies things downstream)
res_WT <- res
res_WT <- res_WT[!is.na(res_WT$meth.wmean.WT),]
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
#ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = plot, width = 6, height = 10, units = "cm")

## -- Calculate DNAme per direction; perform wilcox test
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
