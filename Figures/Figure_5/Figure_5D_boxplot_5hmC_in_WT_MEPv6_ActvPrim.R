## 09.12.2025
## Plotting 5hmC distribution depending on enhancer type
library(ggplot2)
library(dplyr)
#library(tidyr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_Master_Enhancer_Peakome_v6/")

## Read data
All_data <- read.table(file = "WT_D0-D2-D4_hMATAC__Coordinates_WT_Master_Enhancer_Peakome_v6.txt", header = TRUE, sep = '\t')
All_data$D0_class <- factor(All_data$D0_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D2_class <- factor(All_data$D2_class, levels = c("actv", "prim", "inac", "othr"))
All_data$D4_class <- factor(All_data$D4_class, levels = c("actv", "prim", "inac", "othr"))

## Pivot
All_data_long <- tidyr::pivot_longer(data = All_data,
                                     cols = c("meth.wmean_D0_hMATAC", "meth.wmean_D2_hMATAC", "meth.wmean_D4_hMATAC"),
                                     names_to = "Day", values_to = "meth.wmean_hMATAC")
All_data_long <- All_data_long |> filter(!is.na(meth.wmean_hMATAC))

## Parse day
All_data_long$Day <- gsub(pattern = "meth.wmean_", replacement = "", x = All_data_long$Day)
All_data_long$Day <- gsub(pattern = "_hMATAC", replacement = "", x = All_data_long$Day)
All_data_long$Day <- factor(All_data_long$Day, levels = c("D0", "D2", "D4"))

## Add enhancer class
All_data_long$Class <- NA
All_data_long$Class[All_data_long$Day == "D0"] <- as.character(All_data_long$D0_class[All_data_long$Day == "D0"])
All_data_long$Class[All_data_long$Day == "D2"] <- as.character(All_data_long$D2_class[All_data_long$Day == "D2"])
All_data_long$Class[All_data_long$Day == "D4"] <- as.character(All_data_long$D4_class[All_data_long$Day == "D4"])
All_data_long$Class <- factor(All_data_long$Class, levels = c("actv", "prim", "inac", "othr"))

# ## Boxplot w outliers
# ggplot(All_data_long[All_data_long$Class %in% c("actv", "prim"),], aes(x = Day, y = meth.wmean_hMATAC, fill = Class)) +
#   geom_boxplot(na.rm = TRUE, outliers = TRUE, outlier.shape = 21, outlier.size = 1) +
#   scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbf"), drop = FALSE) +
#   theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

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


## -- Statistics -----------------------------

## -- Calculate 5hmC per ehancer class; perform wilcox test

## Summarize results
Summary_hmC <- All_data_long |> group_by(Day, Class) |>
  filter(Class %in% c("actv", "prim")) |>
  summarise("Count" = n(),
            "Mean_5hmC" = round(mean(meth.wmean_hMATAC), digits = 4),
            "Q25_5hmC" = round(quantile(meth.wmean_hMATAC, 0.25), digits = 4),
            "Q50_5hmC" = round(quantile(meth.wmean_hMATAC, 0.50), digits = 4),
            "Q75_5hmC" = round(quantile(meth.wmean_hMATAC, 0.75), digits = 4))

## -- Wilcox testing

## D0
hmC_actv <- All_data_long |> filter(Day == "D0" & Class == "actv") |> pull(meth.wmean_hMATAC)
hmC_prim <- All_data_long |> filter(Day == "D0" & Class == "prim") |> pull(meth.wmean_hMATAC)
T_two <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "two.sided")
T_greater <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "greater")
T_less <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "less")
T_results <- data.frame("Day" = "D0",
                        "Class" = c("actv", "prim"),
                        "Pval.TwoSided" = c(NA, T_two$p.value),
                        "Pval.Greater" = c(NA, T_greater$p.value),
                        "Pval.Less" = c(NA, T_less$p.value))

## D2
hmC_actv <- All_data_long |> filter(Day == "D2" & Class == "actv") |> pull(meth.wmean_hMATAC)
hmC_prim <- All_data_long |> filter(Day == "D2" & Class == "prim") |> pull(meth.wmean_hMATAC)
T_two <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "two.sided")
T_greater <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "greater")
T_less <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "less")
T_results <- rbind(T_results, data.frame("Day" = "D2",
                                         "Class" = c("actv", "prim"),
                                         "Pval.TwoSided" = c(NA, T_two$p.value),
                                         "Pval.Greater" = c(NA, T_greater$p.value),
                                         "Pval.Less" = c(NA, T_less$p.value)))

## D4
hmC_actv <- All_data_long |> filter(Day == "D4" & Class == "actv") |> pull(meth.wmean_hMATAC)
hmC_prim <- All_data_long |> filter(Day == "D4" & Class == "prim") |> pull(meth.wmean_hMATAC)
T_two <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "two.sided")
T_greater <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "greater")
T_less <- wilcox.test(x = hmC_actv, y = hmC_prim, alternative = "less")
T_results <- rbind(T_results, data.frame("Day" = "D4",
                                         "Class" = c("actv", "prim"),
                                         "Pval.TwoSided" = c(NA, T_two$p.value),
                                         "Pval.Greater" = c(NA, T_greater$p.value),
                                         "Pval.Less" = c(NA, T_less$p.value)))

## Append results to summary stats
Summary_hmC <- merge(Summary_hmC, T_results)
## Sort
Summary_hmC <- Summary_hmC |> arrange(Day, Class)

## Export
# write.table(x = Summary_hmC, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE,
#             file = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_5/Panel-D_Boxplots_5hmC_active_vs_primed_enhancers/Info.txt")
