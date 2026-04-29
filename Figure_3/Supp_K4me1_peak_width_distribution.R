## 21.11.2025
## Compare K4me1 RmProm width peaks
## Edit 17.02.2026: Merge reps into a single boxplot
## Edit 17.02.2026: Include mergedreps. Change from TCseq folder to Data/peaks folder
## Edit 26.03.2026: Inverted levels; WT > TKO > DKO
library(ggplot2)
library(stringr)
suppressMessages(library(dplyr))
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Data/peaks/K4me1_ChIPseeker/")
#setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT-DKO-TKO_K4me1_RmProm/peaks/")

## List of peaks
Peak_files <- list.files(path = ".", pattern = "RmProm.filtered.broadPeak", all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE, include.dirs = FALSE, no.. = TRUE)

All_peaks <- data.frame("chr" = character(),
                        "start" = numeric(),
                        "end" = numeric(),
                        "category" = character(),
                        "sample" = character())

for (i in 1:length(Peak_files)) {
  Current_peak_file <- Peak_files[i]
  Code <- sub(pattern = "_trimV6_bsmrk_mm10_mrkdup_RmProm.filtered.broadPeak", replacement = "", x = Current_peak_file)
  Peaks <- read.table(file = Current_peak_file, header = FALSE, sep = '\t', quote = "")
  colnames(Peaks) <- c("chr", "start", "end", "category")
  Peaks$sample <- Code
  All_peaks <- rbind(All_peaks, Peaks)
}
rm(Current_peak_file, i, Code, Peaks)

## Parse sample string
All_peaks$sample <- sub(pattern = "_K4me1_MCnT", replacement = "", x = All_peaks$sample)
Sample_parse <- as.data.frame(str_split_fixed(string = All_peaks$sample, pattern = "_", n = 3))
colnames(Sample_parse) <- c("Gen", "Day", "Rep")
## Parse replicate
Sample_parse$Rep <- sub(x = Sample_parse$Rep, pattern = "_.+", replacement = "", perl = TRUE)
All_peaks <- cbind(All_peaks, Sample_parse)
rm(Sample_parse)
## Replace "240409" for "rep1"
All_peaks$Rep <- gsub(pattern = "240409", replacement = "rep1", x = All_peaks$Rep)
## Clean
All_peaks$Gen <- factor(All_peaks$Gen, levels = c("WT", "TKO", "DKO"))
All_peaks$Day <- factor(All_peaks$Day, levels = c("D0", "D2", "D4"))
All_peaks$Rep <- factor(All_peaks$Rep, levels = c("rep1", "rep2", "mergedreps"))
All_peaks$width <- All_peaks$end - All_peaks$start + 1

## How many peaks per sample?
table(All_peaks$Day, All_peaks$Rep, All_peaks$Gen)

## -- Plotting time

## By genotype
ggplot(All_peaks, aes(x = Day, y = width, group = interaction(Day, Rep), fill = Gen)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = c(WT = "#266799", DKO = "#8c8c8c", TKO = "#F68532")) +
  labs(y = "width (bp)") +
  theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") +
  facet_wrap(~Gen)

## By day - individual replicates
p <- ggplot(All_peaks[All_peaks$Rep != "mergedreps",], aes(x = Gen, y = width/1000, group = interaction(Gen, Rep), fill = Gen)) +
  geom_boxplot(outlier.size = 0.8, outlier.shape = 21) +
  scale_fill_manual(values = c(WT = "#266799", DKO = "#8c8c8c", TKO = "#F68532")) +
  labs(y = "peak width (kb)", x = "Genotype") +
  theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") +
  facet_wrap(~Day)
p
## Export
# ggsave(filename = "Boxplot_K4me1_peak_width_WT-TKO-DKO_D0-D2-D4.pdf", plot = p, device = "pdf", width = 16, height = 10, units = "cm",
#        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_3/Supp_Boxplot_K4me1_peak_width_distribution/")

## By day - only mergedreps
p <- ggplot(All_peaks[All_peaks$Rep == "mergedreps",], aes(x = Gen, y = width/1000, fill = Gen)) +
  geom_boxplot(outlier.size = 0.8, outlier.shape = 21) +
  scale_fill_manual(values = c(WT = "#266799", DKO = "#8c8c8c", TKO = "#F68532")) +
  labs(y = "peak width (kb)", x = "Genotype") +
  theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") +
  facet_wrap(~Day)
p
## Export
# ggsave(filename = "Boxplot_K4me1_peak_width_WT-TKO-DKO_D0-D2-D4_mergedreps.pdf", plot = p, device = "pdf", width = 12, height = 10, units = "cm",
#        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_3/Supp_Boxplot_K4me1_peak_width_distribution/")


## -------------------
## -- Generate summary

Width_summary <- All_peaks |>
  filter(Rep == "mergedreps") |>
  group_by(Gen, Day) |>
  summarise("Num_peaks" = n(),
            "Mean_width" = round(mean(width)),
            "Q25_width" = round(quantile(width, 0.25)),
            "Q50_width" = round(quantile(width, 0.50)),
            "Q75_width" = round(quantile(width, 0.75)))
Width_summary <- Width_summary |> arrange(Day, Gen)
## Export
# write.table(x = Width_summary, append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE,
#             file = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_3/Supp_Boxplot_K4me1_peak_width_distribution/Info_mergedreps.txt")

## -- Calculate wilcox test

## Extract D0 values
WT <- All_peaks |> filter(Rep == "mergedreps" & Gen == "WT" & Day == "D0") |> pull(width)
DKO <- All_peaks |> filter(Rep == "mergedreps" & Gen == "DKO" & Day == "D0") |> pull(width)
TKO <- All_peaks |> filter(Rep == "mergedreps" & Gen == "TKO" & Day == "D0") |> pull(width)
## Two-sided wilcox
test_WT_DKO <- wilcox.test(x = WT, y = DKO)
test_WT_DKO$p.value
test_WT_TKO <- wilcox.test(x = WT, y = TKO)
test_WT_TKO$p.value

## Extract D2 values
WT <- All_peaks |> filter(Rep == "mergedreps" & Gen == "WT" & Day == "D2") |> pull(width)
DKO <- All_peaks |> filter(Rep == "mergedreps" & Gen == "DKO" & Day == "D2") |> pull(width)
TKO <- All_peaks |> filter(Rep == "mergedreps" & Gen == "TKO" & Day == "D2") |> pull(width)
## Two-sided wilcox
test_WT_DKO <- wilcox.test(x = WT, y = DKO)
test_WT_DKO$p.value
test_WT_TKO <- wilcox.test(x = WT, y = TKO)
test_WT_TKO$p.value

## Extract D4 values
WT <- All_peaks |> filter(Rep == "mergedreps" & Gen == "WT" & Day == "D4") |> pull(width)
DKO <- All_peaks |> filter(Rep == "mergedreps" & Gen == "DKO" & Day == "D4") |> pull(width)
TKO <- All_peaks |> filter(Rep == "mergedreps" & Gen == "TKO" & Day == "D4") |> pull(width)
## Two-sided wilcox
test_WT_DKO <- wilcox.test(x = WT, y = DKO)
test_WT_DKO$p.value
test_WT_TKO <- wilcox.test(x = WT, y = TKO)
test_WT_TKO$p.value

