## Alluvial plot of enhancer dynamics: category swapping through differentiation (mESCs -> EpiLCs -> ExEpiLCs)
## Figure 4B
library(ggplot2)
library(ggalluvial)
library(dplyr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/")
gen <- "WT"

## -- Import all enhancers
Name_file <- paste(gen,"_Enhancer_peakome/",gen,"_Master_Enhancer_Peakome_v6_single-class.bed", sep = "")
Enhancers <- read.table(file = Name_file, header = FALSE, sep = '\t')
colnames(Enhancers) <- c("chr", "start", "end", "D0_class", "D2_class", "D4_class")
Enhancers$D0_class <- factor(Enhancers$D0_class, levels = c("actv", "prim", "inac", "othr"))
Enhancers$D2_class <- factor(Enhancers$D2_class, levels = c("actv", "prim", "inac", "othr"))
Enhancers$D4_class <- factor(Enhancers$D4_class, levels = c("actv", "prim", "inac", "othr"))

## -- Enhancer flow between days
## Simplify the data.frame by counting the frequency of each combination of states (4x4x4 - 1; 63 columns instead of >59K)
Cats <- lamisc::count_duplicates(Enhancers[, c("D0_class", "D2_class", "D4_class")])
Cats <- arrange(Cats, D0_class, D2_class, D4_class)

## Plot flow
p <- ggplot(data = Cats, aes(axis1 = D0_class, axis2 = D2_class, axis3 = D4_class, y = dupe_count)) +
  geom_alluvium(aes(fill = D0_class), alpha = 0.85, curve_type = "cubic") +
  scale_fill_manual(values = c("#E5614E", "#F0AC4B", "#77BEF0", "#bfbfbfe6")) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_continuous(breaks = 1:3, labels = c("D0", "D2", "D4")) +
  theme_void() + theme(legend.position = "none")
p

## Export plot
Name_plot <- paste(gen,"_Master_Enhancer_Peakome_v6_Flux.pdf", sep = "")
# ggsave(filename = Name_plot, plot = p, device = "pdf", width = 12, height = 10, units = "cm")

## -- Create summary table
## Pivot to long
Enhancers_long <- Enhancers |> tidyr::pivot_longer(cols = c("D0_class", "D2_class", "D4_class"), names_to = "Day", values_to = "Class")
Enhancers_long$Day <- factor(Enhancers_long$Day, levels = c("D0_class", "D2_class", "D4_class"))
## Group and summarize
Enhancers_summary <- Enhancers_long |> group_by(Day, Class) |> summarise("Count" = n())
## Pivot to wide
Enhancers_summary <- Enhancers_summary |> tidyr::pivot_wider(names_from = Day, values_from = Count)
Enhancers_summary
