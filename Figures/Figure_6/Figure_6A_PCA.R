## RNAseq analysis: PCA and and individual gene plots
## Figure 6A and Supplemental Figure 8A,C
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
## Inital number of genes
n_orig <- dim(ddsHTSeq)[1]
cat("Genes counted: ",n_orig, sep = "")

## Remove genes in which there are no counts for any sample
keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq <- ddsHTSeq[keep,]
n_new <- dim(ddsHTSeq)[1]
cat("Genes retained: ",n_new," (",round(100*n_new/n_orig, digits = 2), "%)", sep = "")
rm(n_orig, n_new, keep)

## -- Calculate model fitting --
dds <- DESeq(ddsHTSeq)

## ----------- ##
## --  PCA  -- ##
## ----------- ##

# The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST
# when there is a wide range of sequencing depth across samples.

## -- Stabilize variance --
rld <- rlog(dds, blind = TRUE)
rld <- plotPCA(rld, intgroup=c("Genotype", "Day"), returnData=TRUE)

## Add factors (make sure KO order matches colors)
rld$Genotype <- factor(rld$Genotype, levels = c("WT","DKO","TKO"))
rld$Day <- factor(rld$Day, levels = c("D0","D2","D4","D7"))
rld$condition <- factor(rld$condition, levels = unique(rld$condition))

PCvar_rld <- round(100 * attr(rld, "percentVar"))

## -- Plot rld --
p <- ggplot(rld, aes(x = PC1, y = PC2, color = condition, shape = Day)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("#78B1C5","#3B90AF","#266799","#154B75",
                                "#a6a6a6","#8c8c8c","#4c4c4c","#0d0d0d",
                                "#FBCC82","#F9BA4E","#f68532","#DF6A16")) +
  scale_shape_manual(values = c(8,16,15,17)) +
  xlab(paste0("PC1: ",PCvar_rld[1],"% variance")) +
  ylab(paste0("PC2: ",PCvar_rld[2],"% variance")) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
p
#ggsave(filename = "PCA_rld_top500.pdf", plot = p, device = "pdf", path = "../../../Plots/htseq_count/RNAseq_v2/", width = 10, height = 10, units = "cm")

## ------------------- ##
## --  Gene counts  -- ##
## ------------------- ##

## -- Specify markers of interest --
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

## -- Extract (normalized) expression of markers --
Counts_MOI <- counts(dds, normalized=TRUE)[MOI, ]
Counts_MOI <- as.data.frame(Counts_MOI) |> tibble::rownames_to_column(var="gene_name")

## -- Convert to tidy format --
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

## -- Calculate mean expression for line fit --
Mean_fit <- Tidy_counts |> group_by(Genotype, Day, gene_name) |> summarise("Counts" = mean(Counts))

## -- Plot and export --
Folder_plot <- "Figure_6/Supp_Pluripotency_markers/"

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
    facet_grid(rows = vars(gene_name), cols = vars(Genotype))
  Name_plot <- paste0("RNAseq_expression_WT-KOs_D0-D7_",Current_gene,".pdf")
  #ggsave(filename = Name_plot, plot = p, device = "pdf", path = Folder_plot, width = 10, height = 4, units = "cm")
}

