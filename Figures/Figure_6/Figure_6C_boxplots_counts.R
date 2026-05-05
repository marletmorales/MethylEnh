## Plotting normalized read counts for genes associated to each type of enhancer
## Figure 6C and Supplemental Figure 8D
library(dplyr)
library(forcats)
library(ggplot2)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/")

## -- Functions ----
Factorize_enhancers <- function(Class_df){
  if(class(Class_df) == "data.frame") {
    ## Assess that the enhancer categories are the only values in input
    if(all(apply(Class_df, 2, function(x) all(x %in% c("actv", "actvM", "prim", "inac", "othr"))))){
      ## Figure out whether we are using 3-class enhancers or 4-class (includes actvM)
      if(sum(Class_df == "actvM") > 0) { ## We are in 4-class
        for(i in 1:ncol(Class_df)) {
          Class_df[,i] <- factor(Class_df[,i], levels = c("actv", "actvM", "prim", "inac", "othr"))
        }
      } else { ## We are in 3-class
        for(i in 1:ncol(Class_df)) {
          Class_df[,i] <- factor(Class_df[,i], levels = c("actv", "prim", "inac", "othr"))
        }
      }
      return(Class_df)
    } else {
      print("Input contains strings outside of enhancer classes (actv,prim,inac,othr). Returning original input")
      return(Class_df)
    }
  } else {
    print("Error. Input must be a data.frame")
  }
}

## -- Part 1: Classify genes to a single enhancer category ----

## -- Read table
Region_gene_parsed <- read.table(file = "WT_Master_Enhancer_Peakome_v6_annotated_TargetGene_wExpression_wActvM.txt", header = TRUE, sep = '\t')
## Remove unused columns
Region_gene_parsed <- Region_gene_parsed[,!colnames(Region_gene_parsed) %in% c("peak_id", "dist_sense", "dist_bp", "WT_D7_avg", "DKO_D7_avg", "TKO_D7_avg")]
## -- Factorize enhancer classes
## 3-class
Region_gene_parsed[,c("D0_class", "D2_class", "D4_class")] <- Factorize_enhancers(Region_gene_parsed[,c("D0_class", "D2_class", "D4_class")])
## 4-class
Region_gene_parsed[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")] <- Factorize_enhancers(Region_gene_parsed[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")])

## -- Classify each gene based on which enhancers is it associated to
## If a gene is linked to multiple enhancers, keep only enhancers of a single class with the following priority:
## actvM = actv > prim > inac

## -- For D0
## DO NOT do actvM (the category is too small and doesn't make sense at D0)
D0_actv <- Region_gene_parsed |> filter(D0_class_ext == "actv") |> pull(gene_symbol) |> unique()
D0_prim <- Region_gene_parsed |> filter(D0_class_ext == "prim") |> pull(gene_symbol) |> unique()
D0_inac <- Region_gene_parsed |> filter(D0_class_ext == "inac") |> pull(gene_symbol) |> unique()
## If already in actv, drop from prim
D0_prim <- D0_prim[!D0_prim %in% c(D0_actv)]
## If already in actv/prim, drop from inac
D0_inac <- D0_inac[!D0_inac %in% c(D0_actv, D0_prim)]
## Convert to df
D0_actv <- data.frame("gene_symbol" = D0_actv, "Day" = "D0", "Class_ext" = "actv")
D0_prim <- data.frame("gene_symbol" = D0_prim, "Day" = "D0", "Class_ext" = "prim")
D0_inac <- data.frame("gene_symbol" = D0_inac, "Day" = "D0", "Class_ext" = "inac")
D0_all <- rbind(D0_actv, D0_prim, D0_inac)
## Drop empty gene
D0_all <- D0_all |> filter(gene_symbol != "")
## Add WT expression
D0_all <- merge(x = D0_all, y = distinct(Region_gene_parsed[,c("gene_symbol", "WT_D0_avg")]), by = "gene_symbol", all.x = TRUE, all.y = FALSE)
colnames(D0_all) <- c("gene_symbol", "Day", "Class_ext", "WT_avg")
## Clear workspace
rm(D0_actv, D0_prim, D0_inac)

## -- For D2
D2_actv <- Region_gene_parsed |> filter(D2_class_ext == "actv") |> pull(gene_symbol) |> unique()
D2_actvM <- Region_gene_parsed |> filter(D2_class_ext == "actvM") |> pull(gene_symbol) |> unique()
D2_prim <- Region_gene_parsed |> filter(D2_class_ext == "prim") |> pull(gene_symbol) |> unique()
D2_inac <- Region_gene_parsed |> filter(D2_class_ext == "inac") |> pull(gene_symbol) |> unique()
## If already in actv/actvM, drop from prim
D2_prim <- D2_prim[!D2_prim %in% c(D2_actv, D2_actvM)]
## If already in actv/actvM/prim, drop from inac
D2_inac <- D2_inac[!D2_inac %in% c(D2_actv, D2_actvM, D2_prim)]
## Convert to df
D2_actv <- data.frame("gene_symbol" = D2_actv, "Day" = "D2", "Class_ext" = "actv")
D2_actvM <- data.frame("gene_symbol" = D2_actvM, "Day" = "D2", "Class_ext" = "actvM")
D2_prim <- data.frame("gene_symbol" = D2_prim, "Day" = "D2", "Class_ext" = "prim")
D2_inac <- data.frame("gene_symbol" = D2_inac, "Day" = "D2", "Class_ext" = "inac")
D2_all <- rbind(D2_actv, D2_actvM, D2_prim, D2_inac)
## Drop empty gene
D2_all <- D2_all |> filter(gene_symbol != "")
## Add WT expression
D2_all <- merge(x = D2_all, y = distinct(Region_gene_parsed[,c("gene_symbol", "WT_D2_avg")]), by = "gene_symbol", all.x = TRUE, all.y = FALSE)
colnames(D2_all) <- c("gene_symbol", "Day", "Class_ext", "WT_avg")
## Clear workspace
rm(D2_actv, D2_actvM, D2_prim, D2_inac)

## -- For D4
D4_actv <- Region_gene_parsed |> filter(D4_class_ext == "actv") |> pull(gene_symbol) |> unique()
D4_actvM <- Region_gene_parsed |> filter(D4_class_ext == "actvM") |> pull(gene_symbol) |> unique()
D4_prim <- Region_gene_parsed |> filter(D4_class_ext == "prim") |> pull(gene_symbol) |> unique()
D4_inac <- Region_gene_parsed |> filter(D4_class_ext == "inac") |> pull(gene_symbol) |> unique()
## If already in actv/actvM, drop from prim
D4_prim <- D4_prim[!D4_prim %in% c(D4_actv, D4_actvM)]
## If already in actv/actvM/prim, drop from inac
D4_inac <- D4_inac[!D4_inac %in% c(D4_actv, D4_actvM, D4_prim)]
## Convert to df
D4_actv <- data.frame("gene_symbol" = D4_actv, "Day" = "D4", "Class_ext" = "actv")
D4_actvM <- data.frame("gene_symbol" = D4_actvM, "Day" = "D4", "Class_ext" = "actvM")
D4_prim <- data.frame("gene_symbol" = D4_prim, "Day" = "D4", "Class_ext" = "prim")
D4_inac <- data.frame("gene_symbol" = D4_inac, "Day" = "D4", "Class_ext" = "inac")
D4_all <- rbind(D4_actv, D4_actvM, D4_prim, D4_inac)
## Drop empty gene
D4_all <- D4_all |> filter(gene_symbol != "")
## Add WT expression
D4_all <- merge(x = D4_all, y = distinct(Region_gene_parsed[,c("gene_symbol", "WT_D4_avg")]), by = "gene_symbol", all.x = TRUE, all.y = FALSE)
colnames(D4_all) <- c("gene_symbol", "Day", "Class_ext", "WT_avg")
## Clear workspace
rm(D4_actv, D4_actvM, D4_prim, D4_inac)


## -- Convert into new table
All_days <- rbind(D0_all, D2_all, D4_all)
All_days$Day <- factor(All_days$Day, levels = c("D0", "D2", "D4"))
All_days$Class_ext <- factor(All_days$Class_ext, levels = c("actv", "actvM", "prim", "inac"))

#rm(D0_all, D2_all, D4_all)

## -- Plot time!

Plot_folder <- "Figure_6/Panel-C_Boxplot_WT_GeneExp_by_EnhancerCat/"

## -- Plot all days
ggplot(All_days, aes(x = Class_ext, y = log2(WT_avg+1), fill = Class_ext)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B", "#77BEF0"), drop = FALSE) +
  labs(x = "Enhancer categories", y = "log2(WT+1)") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) +
  facet_wrap(~Day)

## -- D0
p0b <- ggplot(All_days[All_days$Day == "D0",], aes(x = Class_ext, y = log2(WT_avg+1), fill = Class_ext)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B", "#77BEF0"), drop = FALSE) +
  labs(x = "Enhancer categories", y = "log2(WT+1)") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p0b
#ggsave(filename = "WT_D0_EnhCats_Boxplot.pdf", plot = p0b, device = "pdf", path = Plot_folder, width = 8, height = 8, units = "cm")


## -- D2
p2b <- ggplot(All_days[All_days$Day == "D2",], aes(x = Class_ext, y = log2(WT_avg+1), fill = Class_ext)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B", "#77BEF0"), drop = FALSE) +
  labs(x = "Enhancer categories", y = "log2(WT+1)") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p2b
#ggsave(filename = "WT_D2_EnhCats_Boxplot.pdf", plot = p2b, device = "pdf", path = Plot_folder, width = 8, height = 8, units = "cm")


## -- D4
p4b <- ggplot(All_days[All_days$Day == "D4",], aes(x = Class_ext, y = log2(WT_avg+1), fill = Class_ext)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B", "#77BEF0"), drop = FALSE) +
  labs(x = "Enhancer categories", y = "log2(WT+1)") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
p4b
ggsave(filename = "WT_D4_EnhCats_Boxplot.pdf", plot = p4b, device = "pdf", path = Plot_folder, width = 8, height = 8, units = "cm")


## -- Summary
Cats_summary <- All_days |> group_by(Day, Class_ext) |> summarise("Num_genes" = n(),
                                                                  "Mean_exp" = round(mean(WT_avg), digits = 2),
                                                                  "Q25_exp" = round(quantile(WT_avg, 0.25), digits = 2),
                                                                  "Q50_exp" = round(quantile(WT_avg, 0.50), digits = 2),
                                                                  "Q75_exp" = round(quantile(WT_avg, 0.75), digits = 2))
Cats_summary

## Export table
#write.table(x = Cats_summary, file = paste0(Plot_folder, "Info.txt"), append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)

