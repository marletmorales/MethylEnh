## 2026.01.06
## Compare 5hmC between 5hmCseq (whole-genome) and hMATAC
library(ggplot2)
library(dplyr)
library(bedtoolsr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/")

## -- Read cov5 files
## Whole genome (5hmCseq)
hmC_WG <- read.table(file = "Input_files/cov5/WT_D4_5hmCseq_mergedreps_trimV1_mm10_cov5_annotated.txt", header = TRUE, sep = '\t')
## Keep important columns
hmC_WG <- hmC_WG[,c("seqnames", "start", "end", "perc_meth", "c.meth", "c.unmeth", "annot.type", "coord")]
colnames(hmC_WG) <- c("chr", "start", "end", "perc_meth", "c.meth", "c.unmeth", "annot.type", "coord")
## CpGs WG: 147,908
nrow(hmC_WG)

## ATAC (hmATAC)
hmC_ATAC <- read.table(file = "Input_files/cov5/WT_D4_hMATAC_mergedreps_2507_trimV6_mm10_cov5.txt", header = FALSE, sep = '\t')
colnames(hmC_ATAC) <- c("chr", "start", "end", "perc_meth", "c.meth", "c.unmeth")
## Shift coordinates by 1 bp (WG coordinates were shifted to 1-based coordinate system by GRanges)
hmC_ATAC$start <- hmC_ATAC$start + 1
hmC_ATAC$end <- hmC_ATAC$end + 1
## CpGs ATAC: 3,355,292
nrow(hmC_ATAC)

## -- Filter out blacklisted regions
Blcklist <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Data/blacklist/mm10_blackList_ENCFF547MET.bed", header = FALSE, sep = '\t')
colnames(Blcklist) <- c("chr", "start", "end")

## WG
hmC_WG <- bt.intersect(a = hmC_WG, b = Blcklist, v = TRUE)
colnames(hmC_WG) <- c("chr", "start", "end", "perc_meth", "c.meth", "c.unmeth", "annot.type", "coord")
## CpGs WG: 145,875 (98.6%; 2033 CpGs removed)
nrow(hmC_WG)

## ATAC
hmC_ATAC <- bt.intersect(a = hmC_ATAC, b = Blcklist, v = TRUE)
colnames(hmC_ATAC) <- c("chr", "start", "end", "perc_meth", "c.meth", "c.unmeth")
## CpGs ATAC: 3,352,857 (99.9%; 2435 CpGs removed)
nrow(hmC_ATAC)

## -- Merge dataframes
hmC_shared <- merge(x = hmC_WG, y = hmC_ATAC, all = FALSE, by = c("chr", "start", "end"), suffixes = c(".WG", ".ATAC"))
## Rearrange columns
hmC_shared <- hmC_shared[,c("chr", "start", "end", "annot.type", "perc_meth.WG", "perc_meth.ATAC", "c.meth.WG", "c.unmeth.WG", "c.meth.ATAC", "c.unmeth.ATAC", "coord")]

## Calculate sequencing depth
hmC_shared <- hmC_shared |> mutate("Reads_WG" = c.meth.WG + c.unmeth.WG,
                                   "Reads_ATAC" = c.meth.ATAC + c.unmeth.ATAC,
                                   "FC_reads" = Reads_ATAC/Reads_WG)

## -- Stats

## CpGs WG: 145,875
nrow(hmC_WG)
## CpGs ATAC: 3,352,857
nrow(hmC_ATAC)
## CpGs in both datasets: 42,984
## 29.47% CpGs of WG
## 1.28% CpGs of ATAC
nrow(hmC_shared)
nrow(hmC_shared)/nrow(hmC_WG)*100
nrow(hmC_shared)/nrow(hmC_ATAC)*100

## 26,328 CpGs have 0% 5hmC in both datasets (61.25% of shared CpGs)
sum(hmC_shared$perc_meth.WG == 0 & hmC_shared$perc_meth.ATAC == 0)
sum(hmC_shared$perc_meth.WG == 0 & hmC_shared$perc_meth.ATAC == 0)/nrow(hmC_shared)*100

## Keep only CpGs that are not equal to 0 in both dataset
hmC_shared_not00 <- hmC_shared |> filter(!(perc_meth.WG == 0 & perc_meth.ATAC == 0))

## ------------------- ##
## -- Plotting time -- ##
## ------------------- ##

## ---- Sequencing depth

## -- Original: WG
ggplot(hmC_WG, aes(x = log2(c.meth + c.unmeth))) +
  geom_histogram(binwidth = 0.25, color = "gray40", fill = "gray80") +
  theme_bw()
median(hmC_WG$c.meth + hmC_WG$c.unmeth)
quantile((hmC_WG$c.meth + hmC_WG$c.unmeth), seq(0,1,0.1))
round(quantile((hmC_WG$c.meth + hmC_WG$c.unmeth), seq(0.9,1,0.01)))
sum((hmC_WG$c.meth + hmC_WG$c.unmeth) <= 6)/nrow(hmC_WG)
sum((hmC_WG$c.meth + hmC_WG$c.unmeth) <= 10)/nrow(hmC_WG)

## -- Original: ATAC
ggplot(hmC_ATAC, aes(x = log2(c.meth + c.unmeth))) +
  geom_histogram(binwidth = 0.25, color = "gray40", fill = "gray80") +
  theme_bw()
median(hmC_ATAC$c.meth + hmC_ATAC$c.unmeth)
quantile((hmC_ATAC$c.meth + hmC_ATAC$c.unmeth), seq(0,1,0.1))
sum((hmC_ATAC$c.meth + hmC_ATAC$c.unmeth) <= 6)/nrow(hmC_ATAC)
sum((hmC_ATAC$c.meth + hmC_ATAC$c.unmeth) <= 10)/nrow(hmC_ATAC)

## -- Shared: WG
ggplot(hmC_shared, aes(x = log2(Reads_WG))) +
  geom_histogram(binwidth = 0.25, color = "gray40", fill = "gray80") +
  theme_bw()
median(hmC_shared$Reads_WG)
quantile(hmC_shared$Reads_WG, seq(0,1,0.1))
round(quantile(hmC_shared$Reads_WG, seq(0.9,1,0.01)))
sum(hmC_shared$Reads_WG <= 6)/nrow(hmC_shared)
sum(hmC_shared$Reads_WG <= 10)/nrow(hmC_shared)

## -- Shared: ATAC
ggplot(hmC_shared, aes(x = log2(Reads_ATAC))) +
  geom_histogram(binwidth = 0.25, color = "gray40", fill = "gray80") +
  theme_bw()
median(hmC_shared$Reads_ATAC)
quantile(hmC_shared$Reads_ATAC, seq(0,1,0.1))
round(quantile(hmC_shared$Reads_ATAC, seq(0.9,1,0.01)))
sum(hmC_shared$Reads_ATAC <= 6)/nrow(hmC_shared)
sum(hmC_shared$Reads_ATAC <= 10)/nrow(hmC_shared)


## ---- 5hmC %

## -- All shared CpGs
## Scatter
ggplot(hmC_shared, aes(x = perc_meth.ATAC, y = perc_meth.WG)) +
  geom_point(alpha = 0.1,) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  theme_bw()

## Binned
ggplot(hmC_shared, aes(x = perc_meth.ATAC, y = perc_meth.WG)) +
  stat_bin2d(bins = 40) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,8000), oob = scales::squish) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "5hmC ATAC", y = "5hmC WG") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right")

## Histogram
# ggplot(hmC_shared, aes(x = perc_meth.ATAC - perc_meth.WG)) +
#   geom_histogram(color = "gray40", fill = "gray80", breaks = seq(-100,100,5)) +
#   theme_bw()


## 5hmC delta vs sequencing depth
## [This one is for the paper]
p <- ggplot(hmC_shared, aes(x = log2(Reads_WG), y = perc_meth.ATAC - perc_meth.WG, color = log2(Reads_ATAC))) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 1, size = 0.7) +
  scale_color_viridis_c() +
  labs(color = "log2(hMATAC r.)", x = "log2(5hmCseq reads)", y = "%5hmC (hM-ATAC - 5hmCseq)") +
  coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "right")
p
legend <- cowplot::ggdraw(cowplot::get_legend(p))
legend

## 5hmC delta vs sequencing depth - density plot
ggplot(hmC_shared, aes(x = log2(Reads_WG), y = perc_meth.ATAC - perc_meth.WG)) +
  stat_bin2d(bins = 40) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,3000), oob = scales::squish) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "right", panel.grid = element_blank())

## 5hmC delta vs sequencing depth - inversed
ggplot(hmC_shared, aes(x = log2(Reads_WG), y = log2(Reads_ATAC), color = (perc_meth.ATAC-perc_meth.WG))) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 1, size = 0.7) +
  scale_color_distiller(palette = "RdYlBu", limits = c(-10,10), oob = scales::squish) +
  labs(color = "d5hmC") +
  #coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "right", panel.grid = element_blank())


## -- CpGs not 0% 5hmC in both
## Scatter
ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC, y = perc_meth.WG)) +
  geom_point(alpha = 0.1,) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  theme_bw()

## Binned
ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC, y = perc_meth.WG)) +
  stat_bin2d(bins = 40) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,800), oob = scales::squish) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "5hmC ATAC", y = "5hmC WG") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right")

## 5hmC delta vs sequencing depth
ggplot(hmC_shared_not00, aes(x = log2(Reads_WG), y = perc_meth.ATAC - perc_meth.WG, color = log2(Reads_ATAC))) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 1, size = 0.7) +
  scale_color_viridis_c() +
  labs(color = "log2(ATAC r.)") +
  coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "none")

## 5hmC delta vs sequencing depth - density plot
ggplot(hmC_shared_not00, aes(x = log2(Reads_WG), y = perc_meth.ATAC - perc_meth.WG)) +
  stat_bin2d(bins = 30) +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,1000), oob = scales::squish) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "right", panel.grid = element_blank())

## Histogram
ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC - perc_meth.WG)) +
  geom_histogram(color = "gray40", fill = "gray80", breaks = seq(-100,100,2.5)) +
  geom_vline(xintercept = 0) +
  theme_bw()

sum(hmC_shared_not00$perc_meth.ATAC - hmC_shared_not00$perc_meth.WG == 0)/nrow(hmC_shared_not00)
sum(hmC_shared_not00$perc_meth.ATAC - hmC_shared_not00$perc_meth.WG > 0)/nrow(hmC_shared_not00)
sum(hmC_shared_not00$perc_meth.ATAC - hmC_shared_not00$perc_meth.WG < 0)/nrow(hmC_shared_not00)

## 5hmC delta vs sequencing depth - inversed
ggplot(hmC_shared_not00, aes(x = log2(Reads_WG), y = log2(Reads_ATAC), color = (perc_meth.ATAC-perc_meth.WG))) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 1, size = 0.7) +
  scale_color_distiller(palette = "RdYlBu", limits = c(-10,10), oob = scales::squish) +
  labs(color = "d5hmC") +
  #coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "right", panel.grid = element_blank())



# ## -- Difference in 5hmC delta vs FC sequencing depth
#
# ## Scatter
# ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC - perc_meth.WG, y = log2(FC_reads))) +
#   geom_point(alpha = 0.1, size = 0.8) +
#   coord_cartesian(xlim = c(-100,100)) +
#   labs(y = "log2([reads ATAC/WG])") +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   theme_bw() + theme(panel.grid = element_blank())
#
# ## Binned
# ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC - perc_meth.WG, y = log2(FC_reads))) +
#   stat_bin2d(bins = 40) +
#   scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0,400), oob = scales::squish) +
#   coord_cartesian(xlim = c(-100,100)) +
#   labs(y = "log2([reads ATAC/WG])") +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   theme_bw() + theme(panel.grid = element_blank())
#
# ## Scatter - y-axis depth WG
# ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC - perc_meth.WG, y = log2(Reads_WG))) +
#   geom_point(alpha = 0.1, size = 0.8) +
#   coord_cartesian(xlim = c(-100,100)) +
#   geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.2) +
#   geom_vline(xintercept = 0) +
#   theme_bw() + theme(panel.grid = element_blank())
#
# ## Scatter - y-axis depth ATAC
# ggplot(hmC_shared_not00, aes(x = perc_meth.ATAC - perc_meth.WG, y = log2(Reads_ATAC))) +
#   geom_point(alpha = 0.1, size = 0.8) +
#   coord_cartesian(xlim = c(-100,100)) +
#   geom_hline(yintercept = 8, linetype = "dashed", linewidth = 0.2) +
#   geom_vline(xintercept = 0) +
#   theme_bw() + theme(panel.grid = element_blank())
#

## -----------------------------
## Plot 5hmC distribution per genomic compartment

## Add bin of 5hmC % - WG
tmp_WG <- hmC_shared$perc_meth.WG
hmC_shared$perc_meth_bin.WG <- case_when(tmp_WG == 0 ~ "0",
                                         tmp_WG > 0 & tmp_WG <= 10 ~ "0-10",
                                         tmp_WG > 10 & tmp_WG <= 20 ~ "10-20",
                                         tmp_WG > 20 & tmp_WG <= 40 ~ "20-40",
                                         tmp_WG > 40 & tmp_WG <= 60 ~ "40-60",
                                         tmp_WG > 60 & tmp_WG <= 80 ~ "60-80",
                                         tmp_WG > 80 & tmp_WG <= 100 ~ "80-100",
                                         .default = NA)
hmC_shared <- hmC_shared |> relocate(perc_meth_bin.WG, .after = perc_meth.WG)
rm(tmp_WG)
## Factorize
hmC_shared$perc_meth_bin.WG <- factor(hmC_shared$perc_meth_bin.WG, levels = c("0", "0-10", "10-20", "20-40", "40-60", "60-80", "80-100"))

## Add bin of 5hmC % - ATAC
tmp_ATAC <- hmC_shared$perc_meth.ATAC
hmC_shared$perc_meth_bin.ATAC <- case_when(tmp_ATAC == 0 ~ "0",
                                           tmp_ATAC > 0 & tmp_ATAC <= 10 ~ "0-10",
                                           tmp_ATAC > 10 & tmp_ATAC <= 20 ~ "10-20",
                                           tmp_ATAC > 20 & tmp_ATAC <= 40 ~ "20-40",
                                           tmp_ATAC > 40 & tmp_ATAC <= 60 ~ "40-60",
                                           tmp_ATAC > 60 & tmp_ATAC <= 80 ~ "60-80",
                                           tmp_ATAC > 80 & tmp_ATAC <= 100 ~ "80-100",
                                           .default = NA)
hmC_shared <- hmC_shared |> relocate(perc_meth_bin.ATAC, .after = perc_meth.ATAC)
rm(tmp_ATAC)
## Factorize
hmC_shared$perc_meth_bin.ATAC <- factor(hmC_shared$perc_meth_bin.ATAC, levels = c("0", "0-10", "10-20", "20-40", "40-60", "60-80", "80-100"))

## -- Plot stacked bars
## Order annotation categories

Stats_summary <- hmC_shared |> group_by(annot.type) |> summarise("mean_WG" = mean(perc_meth.WG),
                                                                 "mean_ATAC" = mean(perc_meth.ATAC))
Stats_summary$mean_Both <- rowMeans(Stats_summary[,c("mean_WG", "mean_ATAC")])

New_annot_order <- Stats_summary |> arrange(mean_WG)
#New_annot_order <- Stats_summary |> arrange(mean_Both)
New_annot_order <- as.character(New_annot_order$annot.type)
hmC_shared$annot.type <- factor(hmC_shared$annot.type, levels = New_annot_order)

## WG
## [For paper?]
ggplot(hmC_shared, aes(x = annot.type, fill = perc_meth_bin.WG)) +
  geom_bar(stat="count", position = position_fill(reverse = TRUE)) +
  #scale_fill_discrete(palette = "YlGnBu") +
  scale_fill_manual(values = scales::pal_brewer(palette = "BuPu")(9)[3:9]) +
  #coord_cartesian(ylim = c(.4,1)) +
  labs(fill = "% 5hmC", y = "% CpGs", x = "Genomic category") +
  scale_x_discrete(labels = sub(pattern = "mm10_(genes_)*", replacement = "", x = levels(hmC_shared$annot.type))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.9), panel.grid = element_blank())

## Zoomed-in:
ggplot(hmC_shared, aes(x = annot.type, fill = perc_meth_bin.WG)) +
  geom_bar(stat="count", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = scales::pal_brewer(palette = "BuPu")(9)[3:9]) +
  coord_cartesian(ylim = c(.5,1)) +
  labs(fill = "% 5hmC", y = "% CpGs", x = "Genomic category") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_discrete(labels = sub(pattern = "mm10_(genes_)*", replacement = "", x = levels(hmC_shared$annot.type))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.9), panel.grid = element_blank())


## ATAC
ggplot(hmC_shared, aes(x = annot.type, fill = perc_meth_bin.ATAC)) +
  geom_bar(stat="count", position = position_fill(reverse = TRUE)) +
  #scale_fill_discrete(palette = "YlGnBu") +
  scale_fill_manual(values = scales::pal_brewer(palette = "BuPu")(9)[3:9]) +
  #coord_cartesian(ylim = c(.75,1)) +
  labs(fill = "% 5hmC", y = "% CpGs", x = "Genomic category") +
  scale_x_discrete(labels = sub(pattern = "mm10_(genes_)*", replacement = "", x = levels(hmC_shared$annot.type))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.9), panel.grid = element_blank())
