## Comparing 5hmC between 5hmCseq (whole-genome) and hMATAC (accessible-only) at single-CpG resolution
## Supplemental Figure 7C
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
#hmC_shared_not00 <- hmC_shared |> filter(!(perc_meth.WG == 0 & perc_meth.ATAC == 0))

## ------------------- ##
## -- Plotting time -- ##
## ------------------- ##

## 5hmC delta vs sequencing depth
p <- ggplot(hmC_shared, aes(x = log2(Reads_WG), y = perc_meth.ATAC - perc_meth.WG, color = log2(Reads_ATAC))) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 1, size = 0.7) +
  scale_color_viridis_c() +
  labs(color = "log2(hMATAC r.)", x = "log2(5hmCseq reads)", y = "delta %5hmC (hM-ATAC - 5hmCseq)") +
  coord_cartesian(ylim = c(-100,100)) +
  theme_bw() + theme(legend.position = "right")
p
legend <- cowplot::ggdraw(cowplot::get_legend(p))
legend

