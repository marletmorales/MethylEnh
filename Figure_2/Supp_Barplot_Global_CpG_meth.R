## 2026.02.23
## Barplot of global CpG metylation levels
library(ggplot2)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/Mean_DNAme_reports/")

CpG_meth <- read.table(file = "Mean_DNAme_Sample_summary.txt", header = TRUE, sep = '\t')
## Factorize
CpG_meth$Gen <- factor(CpG_meth$Gen, levels = c("WT", "DKO", "TKO"))
CpG_meth$Day <- factor(CpG_meth$Day, levels = c("D0", "D2", "D4"))
CpG_meth$Coverage <- factor(CpG_meth$Coverage, levels = c("cov5", "cov1"))

## Dot plot
ggplot(CpG_meth[CpG_meth$Coverage == "cov5",], aes(x = Day, y = CpG.meth, color = Gen, group = Gen, shape = Gen)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("#266799", "#4c4c4c")) +
  scale_shape_manual(values = c(16, 15)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_discrete(expand = c(0.15,0.15)) +
  scale_y_continuous(expand = c(0.02,0.02)) +
  labs(y = "Global 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

## Bar plot - side
p <- ggplot(CpG_meth[CpG_meth$Coverage == "cov5",], aes(x = Day, y = CpG.meth, fill = Gen, group = Gen)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_manual(values = c("#266799", "#4c4c4c")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  labs(y = "Global 5mC %") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p
# ggsave(filename = "Barplot_5mC_WT-DKO_D0-D2-D4_WGBS.pdf", plot = p, device = "pdf", width = 5.5, height = 8.5, units = "cm",
#        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_2/Supp_Barplot_Global_5mC_WGBS/")
