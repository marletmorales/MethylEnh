## 16.02.2026
## Exporting pluripotency markers expression
## Edit 2026.03.19: Changing order WT-TKO-DKO
library(ggplot2)
library(dplyr)
library(tidyr)
library(DESeq2)
setwd(dir = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/htseq_count/RNAseq_v2/")

## -- Read experimental design --
sampleTable <- read.table(file = "ExperimentDesign.txt", header = TRUE, sep = '\t')
sampleTable$condition <- factor(sampleTable$condition, levels = unique(sampleTable$condition))

## -- Read counts --
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition)

## -- Filter low-count genes --
## Remove genes in which there are no counts for any sample
keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq <- ddsHTSeq[keep,]
rm(keep)

## Calculate model fitting
dds <- DESeq(ddsHTSeq)

## -----------------------------
## -- Read markers of interest
MOI <- c("Dnmt1",
         "Dnmt3a",
         "Dnmt3b",
         "Tet1",
         "Tet2",
         "Tet3",
         "Prdm14",
         "Esrrb",
         "Zfp42",
         "Nanog",
         "Dnmt3l",
         "Pou3f1",
         "Otx2",
         "Pou5f1",
         "Sox2",
         "T")


## Extract (normalized) expression of our markers
Counts_MOI <- counts(dds, normalized=TRUE)[MOI, ]

#tmp <- rownames(counts(dds, normalized=TRUE))

Counts_MOI <- as.data.frame(Counts_MOI) |> tibble::rownames_to_column(var="gene_name")

## Convert to tidy format
Tidy_counts <- pivot_longer(data = Counts_MOI,
                            cols = colnames(Counts_MOI)[2:ncol(Counts_MOI)],
                            names_to = "Sample",
                            values_to = "Counts")
## Add metadata
MD <- stringr::str_split_fixed(string = Tidy_counts$Sample, pattern = "_", n = 3)
colnames(MD) <- c("Genotype", "Day", "Rep")
Tidy_counts <- cbind(Tidy_counts, MD)
Tidy_counts <- Tidy_counts[,c("Sample", "Genotype", "Day", "Rep", "gene_name", "Counts")]
rm(MD)

## Specify factors
Tidy_counts$Genotype <- factor(Tidy_counts$Genotype, levels = c("WT","TKO","DKO"))
Tidy_counts$Day <- factor(Tidy_counts$Day, levels = c("D0","D2","D4","D7"))
Tidy_counts$gene_name <- factor(Tidy_counts$gene_name, levels = MOI)

## Calculate mean expression for line fit
Mean_fit <- Tidy_counts |> group_by(Genotype, Day, gene_name) |> summarise("Counts" = mean(Counts))

## ------------------------- ##
## -- By individual genes -- ##
## -----------------=------- ##

Folder_plot <- "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_6/Supp_Pluripotency_markers/"

## Export dimensions: 380 x (80n + 20)

## -- Gene deletions
## Raw counts

for(i in 1:length(MOI)) {
  Current_gene <- MOI[i]
  Sub_table <- Tidy_counts |> filter(gene_name == Current_gene)
  Sub_mean <- Mean_fit |> filter(gene_name == Current_gene)

  p <- ggplot(Sub_table, aes(x = Day, y = Counts, color = Genotype, group = Genotype)) +
    geom_point(alpha = 0.8) +
    geom_line(data = Sub_mean) +
    scale_color_manual(values = c(WT = "#266799", DKO = "#4c4c4c", TKO = "#f68532")) +
    labs(y = "Counts") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_blank()) +
    facet_grid(rows = vars(gene_name), cols = vars(Genotype), scales = "free_y")
  Name_plot <- paste0("RNAseq_expression_WT-KOs_D0-D7_",Current_gene,".pdf")
  ggsave(filename = Name_plot, plot = p, device = "pdf", path = Folder_plot, width = 10, height = 4, units = "cm")
}


