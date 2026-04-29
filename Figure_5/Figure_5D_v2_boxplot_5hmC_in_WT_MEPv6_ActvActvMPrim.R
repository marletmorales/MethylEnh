## 2025.12.09
## Plotting 5hmC distribution depending on enhancer type
## Edit 2026.03.31: Added actvM category
library(ggplot2)
library(dplyr)
#library(tidyr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_Master_Enhancer_Peakome_v6/")

## Read data and merge
All_data <- read.table(file = "WT_D0-D2-D4_hMATAC__Coordinates_WT_Master_Enhancer_Peakome_v6.txt", header = TRUE, sep = '\t')
Ext_class <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_single-class_wActvM.txt", header = TRUE, sep = '\t')
All_data <- left_join(x = All_data, y = Ext_class)
rm(Ext_class)

## Factorize enhancer classes
All_data$D0_class <- factor(All_data$D0_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D2_class <- factor(All_data$D2_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D4_class <- factor(All_data$D4_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D0_class_ext <- factor(All_data$D0_class_ext, levels = c("actv", "actvM", "prim", "inac", "othr"))
All_data$D2_class_ext <- factor(All_data$D2_class_ext, levels = c("actv", "actvM", "prim", "inac", "othr"))
All_data$D4_class_ext <- factor(All_data$D4_class_ext, levels = c("actv", "actvM", "prim", "inac", "othr"))

## Pivot
All_data_long <- tidyr::pivot_longer(data = All_data,
                                     cols = c("meth.wmean_D0_hMATAC", "meth.wmean_D2_hMATAC", "meth.wmean_D4_hMATAC"),
                                     names_to = "Day", values_to = "meth.wmean_hMATAC")
All_data_long <- All_data_long |> filter(!is.na(meth.wmean_hMATAC))

## Parse day
All_data_long$Day <- gsub(pattern = "meth.wmean_", replacement = "", x = All_data_long$Day)
All_data_long$Day <- gsub(pattern = "_hMATAC", replacement = "", x = All_data_long$Day)
All_data_long$Day <- factor(All_data_long$Day, levels = c("D0", "D2", "D4"))

## Add enhancer class as a variable
All_data_long <- All_data_long |> mutate(Class = case_when(
  Day == "D0" ~ D0_class,
  Day == "D2" ~ D2_class,
  Day == "D4" ~ D4_class,
  .default = NA))
All_data_long <- All_data_long |> mutate(Class_ext = case_when(
  Day == "D0" ~ D0_class_ext,
  Day == "D2" ~ D2_class_ext,
  Day == "D4" ~ D4_class_ext,
  .default = NA))

## For those enhancers that are classified as actvM in D0, replace Class_ext with NA (we are omitting this category at D0)
All_data_long$Class_ext[All_data_long$Day == "D0" & All_data_long$Class_ext == "actvM"] <- NA

## Factorize
All_data_long$Class <- factor(All_data_long$Class, levels = c("actv", "prim", "inac", "othr"))
All_data_long$Class_ext <- factor(All_data_long$Class_ext, levels = c("actv", "actvM", "prim", "inac", "othr"))


## ---- Using Class --
## Boxplot w outliers
ggplot(All_data_long[All_data_long$Class %in% c("actv", "prim"),], aes(x = Day, y = meth.wmean_hMATAC, fill = Class)) +
  geom_boxplot(na.rm = TRUE, outliers = TRUE, outlier.shape = 21, outlier.size = 1) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B")) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

## Boxplot wo outliers
p <- ggplot(All_data_long[All_data_long$Class %in% c("actv", "prim"),], aes(x = Day, y = meth.wmean_hMATAC, fill = Class)) +
  geom_boxplot(na.rm = TRUE, outliers = FALSE) +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B")) +
  labs(y = "5hmC %") +
  coord_cartesian(ylim = c(0,30)) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p
# ggsave(filename = "WT_D0-D2-D4_5hmC_hMATAC_per_enhancer_type.pdf", plot = p, device = "pdf", width = 9, height = 9, units = "cm",
#        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_5/Panel-D_Boxplots_5hmC_active_vs_primed_enhancers/")


## ---- Using Class_ext --
## Boxplot w outliers
ggplot(All_data_long[All_data_long$Class_ext %in% c("actv", "actvM", "prim"),], aes(x = Day, y = meth.wmean_hMATAC, fill = Class_ext)) +
  geom_boxplot(na.rm = TRUE, outliers = TRUE, outlier.shape = 21, outlier.size = 1) +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B")) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

## Boxplot wo outliers
p <- ggplot(All_data_long[All_data_long$Class_ext %in% c("actv", "actvM", "prim"),], aes(x = Day, y = meth.wmean_hMATAC, fill = Class_ext)) +
  geom_boxplot(na.rm = TRUE, outliers = FALSE) +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B")) +
  labs(y = "5hmC %") +
  coord_cartesian(ylim = c(0,30)) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p
# ggsave(filename = "WT_D0-D2-D4_5hmC_hMATAC_per_enhancer_type_wActvM.pdf", plot = p, device = "pdf", width = 9, height = 9, units = "cm",
#        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_5/Panel-D_Boxplots_5hmC_active_vs_primed_enhancers/")


## -- Statistics -----------------------------

## ---- Calculate 5hmC per ehancer class; perform wilcox test
## Summarize results
Summary_hmC <- All_data_long |> group_by(Day, Class_ext) |>
  filter(Class_ext %in% c("actv", "actvM", "prim")) |>
  summarise("Count" = n(),
            "Mean_5hmC" = round(mean(meth.wmean_hMATAC), digits = 2),
            "Q25_5hmC" = round(quantile(meth.wmean_hMATAC, 0.25), digits = 2),
            "Q50_5hmC" = round(quantile(meth.wmean_hMATAC, 0.50), digits = 2),
            "Q75_5hmC" = round(quantile(meth.wmean_hMATAC, 0.75), digits = 2))

## ---- Wilcox testing
## -- D0
hmC_actv <- All_data_long |> filter(Day == "D0" & Class_ext == "actv") |> pull(meth.wmean_hMATAC)
hmC_prim <- All_data_long |> filter(Day == "D0" & Class_ext == "prim") |> pull(meth.wmean_hMATAC)

t.a_p <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "two.sided")

t_Output <- data.frame("Day" = "D0",
                       "GroupA" = "actv",
                       "GroupB" = "prim",
                       "A_n" = length(hmC_actv),
                       "B_n" = length(hmC_prim),
                       "A_median_5hmC" = round(median(hmC_actv), digits = 4),
                       "B_median_5hmC" = round(median(hmC_prim), digits = 4),
                       "pval" = t.a_p$p.value)

## -- D2
hmC_actv <- All_data_long |> filter(Day == "D2" & Class_ext == "actv") |> pull(meth.wmean_hMATAC)
hmC_actvM <- All_data_long |> filter(Day == "D2" & Class_ext == "actvM") |> pull(meth.wmean_hMATAC)
hmC_prim <- All_data_long |> filter(Day == "D2" & Class_ext == "prim") |> pull(meth.wmean_hMATAC)

t.a_aM <- wilcox.test(x = hmC_actv, y = hmC_actvM, alternative = "two.sided")
t.aM_p <- wilcox.test(x = hmC_actvM, y = hmC_prim, alternative = "two.sided")
t.a_p <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "two.sided")

t_Output <- rbind(t_Output, data.frame("Day" = "D2",
                                       "GroupA" = c("actv", "actvM", "actv"),
                                       "GroupB" = c("actvM", "prim", "prim"),
                                       "A_n" = c(length(hmC_actv), length(hmC_actvM), length(hmC_actv)),
                                       "B_n" = c(length(hmC_actvM), length(hmC_prim), length(hmC_prim)),
                                       "A_median_5hmC" = round(c(median(hmC_actv), median(hmC_actvM), median(hmC_actv)), digits = 4),
                                       "B_median_5hmC" = round(c(median(hmC_actvM), median(hmC_prim), median(hmC_prim)), digits = 4),
                                       "pval" = c(t.a_aM$p.value, t.aM_p$p.value, t.a_p$p.value)))

## -- D4
hmC_actv <- All_data_long |> filter(Day == "D4" & Class_ext == "actv") |> pull(meth.wmean_hMATAC)
hmC_actvM <- All_data_long |> filter(Day == "D4" & Class_ext == "actvM") |> pull(meth.wmean_hMATAC)
hmC_prim <- All_data_long |> filter(Day == "D4" & Class_ext == "prim") |> pull(meth.wmean_hMATAC)

t.a_aM <- wilcox.test(x = hmC_actv, y = hmC_actvM, alternative = "two.sided")
t.aM_p <- wilcox.test(x = hmC_actvM, y = hmC_prim, alternative = "two.sided")
t.a_p <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "two.sided")

t_Output <- rbind(t_Output, data.frame("Day" = "D4",
                                       "GroupA" = c("actv", "actvM", "actv"),
                                       "GroupB" = c("actvM", "prim", "prim"),
                                       "A_n" = c(length(hmC_actv), length(hmC_actvM), length(hmC_actv)),
                                       "B_n" = c(length(hmC_actvM), length(hmC_prim), length(hmC_prim)),
                                       "A_median_5hmC" = round(c(median(hmC_actv), median(hmC_actvM), median(hmC_actv)), digits = 4),
                                       "B_median_5hmC" = round(c(median(hmC_actvM), median(hmC_prim), median(hmC_prim)), digits = 4),
                                       "pval" = c(t.a_aM$p.value, t.aM_p$p.value, t.a_p$p.value)))

## -- Export results
Name_file <- "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_5/Panel-D_Boxplots_5hmC_active_vs_primed_enhancers/Info_wActvM.txt"
# cat("Summary stats\n", file = Name_file)
# write.table(x = Summary_hmC, file = Name_file, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE)
# cat("\n\nWilcox test\n", file = Name_file, append = TRUE)
# write.table(x = t_Output, file = Name_file, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE)
