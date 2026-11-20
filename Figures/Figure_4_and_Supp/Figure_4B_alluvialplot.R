## Alluvial plot of enhancer dynamics: category swapping through differentiation (mESCs -> EpiLCs -> ExEpiLCs)
## Figure 4B
## Edit 2026.11.17: Changed input file. Changed color coding. Polished generation of summary.
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(lamisc)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/")

## -- Import all enhancers
Enhancers <- read.table(file = "WT_Master_Enhancer_Peakome_v6_single-class_wActvM_v2.txt", header = TRUE, sep = '\t')
Enhancers$D0_class <- factor(Enhancers$D0_class, levels = c("actv", "prim", "inac", "othr"))
Enhancers$D2_class <- factor(Enhancers$D2_class, levels = c("actv", "prim", "inac", "othr"))
Enhancers$D4_class <- factor(Enhancers$D4_class, levels = c("actv", "prim", "inac", "othr"))

## -- Enhancer flow between days
## Simplify the data.frame by counting the frequency of each combination of states (4x4x4 - 1; 63 columns instead of >59K)
Cats <- count_duplicates(Enhancers[, c("D0_class", "D2_class", "D4_class")])
Cats <- arrange(Cats, D0_class, D2_class, D4_class)

## Plot flow
Color_scale <- c("actv" = "#b1c848",
                 "actvM" = "#2cb11b",
                 "prim" = "#1bb6af",
                 "inac" = "#485d7b",
                 "othr" = "#d7d7d7")

p <- ggplot(data = Cats, aes(axis1 = D0_class, axis2 = D2_class, axis3 = D4_class, y = dupe_count)) +
  geom_alluvium(aes(fill = D0_class), alpha = 0.85, curve_type = "cubic") +
  scale_fill_manual(values = Color_scale) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_continuous(breaks = 1:3, labels = c("D0", "D2", "D4")) +
  theme_void() + theme(legend.position = "none")
p

## Visualize color-blind
#colorblindr::cvd_grid(p)

## Export plot
# ggsave(filename = "WT_Master_Enhancer_Peakome_v6_Flux_v2.pdf", plot = p, width = 10, height = 10, units = "cm",
#        path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Paper_additional_files/Figures_Revision_1/Figure_4/Panel-B_Alluvial_plots/")

## -- Create summary table
## Pivot to long
Enh_summary <- Enhancers |>
  select(chr:D4_class) |>
  tidyr::pivot_longer(cols = c("D0_class", "D2_class", "D4_class"), names_to = "Day", values_to = "Class")
Enh_summary$Day <- sub(pattern = "_class", replacement = "", x = Enh_summary$Day)
Enh_summary$Day <- factor(Enh_summary$Day, levels = c("D0", "D2", "D4"))
## Group and summarize
Enh_summary <- Enh_summary |> group_by(Day, Class) |> summarise("Count" = n())
## Pivot to wide
Enh_summary <- Enh_summary |> tidyr::pivot_wider(names_from = Day, values_from = Count)
Enh_summary

## Export
# write.table(x = Enh_summary, sep = "\t", eol = "\n", na = "", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE,
#             file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Paper_additional_files/Figures_Revision_1/Figure_4/Panel-B_Alluvial_plots/Info.txt")
