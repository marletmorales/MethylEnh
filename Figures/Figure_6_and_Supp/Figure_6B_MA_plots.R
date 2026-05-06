## RNAseq differential expression (DE) with DESeq2
## Figure 6B and and Supplemental Figure 8B
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(DESeq2))
#setwd(dir = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/htseq_count/RNAseq_v2/")
setwd(dir = "/Users/marle/Documents/Greenberg_Lab/RNAseq_v2/")

## -------------------------------------------
## -- Specify contrast for Diff Exp analysis

## Make all contrasts
G1 <- c("WT_D0", "WT_D2", "WT_D4", "WT_D7",
        "WT_D0", "WT_D2", "WT_D4", "WT_D7")
G2 <- c("DKO_D0", "DKO_D2", "DKO_D4", "DKO_D7",
        "TKO_D0", "TKO_D2", "TKO_D4", "TKO_D7")
G1_G2 <- data.frame("Group_1" = G1, "Group_2" = G2)

## Iterate over every contrast combination
for(i in 1:nrow(G1_G2)){

  group_1 <- G1_G2$Group_1[i]
  group_2 <- G1_G2$Group_2[i]
  cat("Comparing ",group_1," vs ",group_2,"\n",sep = "")

  ## Specify log file and Plot_folder
  Name_log <- paste("log_files/",group_1,"_vs_",group_2,"_RNAseq_stats.txt", sep = "")
  cat("Comparing ",group_1," vs ",group_2,"\n\n",file = Name_log, sep = "", append = FALSE)
  Plot_folder <- "../../../Plots/htseq_count/RNAseq_v2/"

  ## -----------------
  ## -- Read samples

  ## -- Read experimental design --
  sampleTable <- read.table(file = "ExperimentDesign.txt", header = TRUE, sep = '\t')
  sampleTable$condition <- factor(sampleTable$condition, levels = unique(sampleTable$condition))

  ## -- Read counts --
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition)

  ## -- Filter low-count genes --
  ## Inital number of genes
  n_orig <- dim(ddsHTSeq)[1]
  cat("Genes quantified:\t",n_orig, "\n",file = Name_log, sep = "", append = TRUE)

  ## Remove genes in which there are no counts for any sample
  keep <- rowSums(counts(ddsHTSeq)) > 0
  ddsHTSeq <- ddsHTSeq[keep,]
  n_new <- dim(ddsHTSeq)[1]
  cat("Pass low-count filter:\t",n_new," (",round(100*n_new/n_orig, digits = 2), "%)\n\n", file = Name_log, sep = "", append = TRUE)
  rm(n_orig, n_new, keep)

  ## -- Model fitting
  ddsHTSeq <- DESeq(ddsHTSeq)
  ## Export size factors
  tmp <- sizeFactors(ddsHTSeq)
  tmp <- sort(tmp)
  cat("Size factors range\n", file = Name_log, sep = "", append = TRUE)
  cat("Min: ",names(tmp)[1]," (",tmp[1],")\n",file = Name_log, sep = "", append = TRUE)
  cat("Max: ",names(tmp)[length(tmp)]," (",tmp[length(tmp)],")\n\n",file = Name_log, sep = "", append = TRUE)
  rm(tmp)

  ## ----------------------------- ##
  ## -- Differential Expression -- ##
  ## ----------------------------- ##

  res <- results(ddsHTSeq, alpha = 0.05, contrast = c("condition",group_2,group_1))

  ## Extract stats and remove rows with padj = NA (they were not used for Bonferroni correction either way)
  res_df <- as.data.frame(res)
  res_df <- res_df |> filter(!is.na(padj))
  cat("\nRemaining genes:\t",nrow(res_df),"\n",file = Name_log, sep = "", append = TRUE)

  ## Mark differentially-expressed
  res_df$direction <- "none"
  res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "up"
  res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "down"
  res_df$direction <- factor(res_df$direction, levels = c("down", "none", "up"))
  res_df <- res_df |> tibble::rownames_to_column(var = "gene_name")

  tmp_none <- sum(res_df$direction == "none")
  tmp_up <- sum(res_df$direction == "up")
  tmp_down <- sum(res_df$direction == "down")
  cat("Not Diff Expr:\t\t",tmp_none," (",round(tmp_none*100/nrow(res_df), 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("Up in ",group_2,":\t\t",tmp_up," (",round(tmp_up*100/nrow(res_df), 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("Down in ",group_2,":\t\t",tmp_down," (",round(tmp_down*100/nrow(res_df), 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  rm(tmp_none, tmp_up, tmp_down)

  ## Export in nice format
  tmp <- res_df
  tmp[,c("baseMean", "log2FoldChange", "lfcSE", "stat")] <- round(tmp[,c("baseMean", "log2FoldChange", "lfcSE", "stat")], digits = 4)
  tmp[,c("pvalue", "padj")] <- signif(tmp[,c("pvalue","padj")], digits = 5)
  Name_file <- paste("results/",group_1,"_vs_",group_2,"_RNAseq_DESeq2_DE.txt", sep = "")
  write.table(x = tmp, file = Name_file, append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
  rm(tmp)

  ## -- Plots

  ## Re-store res to keep low_count genes
  res_full <- as.data.frame(res)
  res_full$direction <- "none"
  res_full$direction[res_full$padj < 0.05 & res_full$log2FoldChange > 1] <- "up"
  res_full$direction[res_full$padj < 0.05 & res_full$log2FoldChange < -1] <- "down"
  res_full$direction <- factor(res_full$direction, levels = c("down", "none", "up"))
  res_full <- res_full |> tibble::rownames_to_column(var = "gene_name")
  lim_fc <- max(abs(res_full$log2FoldChange))

  ## -- Volcano
  p <- ggplot(res_full, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
    geom_point(size = 1, alpha = 0.4, na.rm = TRUE) +
    scale_color_manual(values = c("#00BFC4","gray30","#F8766D")) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
    geom_vline(xintercept = c(-1,1), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
    coord_cartesian(xlim = c(-lim_fc, lim_fc)) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")

  Name_plot <- paste("png/",group_1,"_vs_",group_2,"_logFC_padj_volcano.png", sep = "")
  #ggsave(filename = Name_plot, device = "png", path = Plot_folder, plot = p, width = 12, height = 12, units = "cm")
  
  ## -- MA plot
  up_col <- "#266799"
  if(grepl(pattern = "DKO", x = group_2)){ up_col <- "#4c4c4c"
  } else if(grepl(pattern = "TKO", x = group_2)){ up_col <- "#f68532" }

  p <- ggplot(res_full, aes(x = log10(baseMean), y = log2FoldChange, color = direction)) +
    geom_point(size = 0.5, alpha = 1) +
    scale_color_manual(values = c("#266799","gray80",up_col)) +
    coord_cartesian(ylim = c(-lim_fc, lim_fc)) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")

  Name_plot <- paste("png/",group_1,"_vs_",group_2,"_meanExp_logFC_MA.png", sep = "")
  #ggsave(filename = Name_plot, device = "png", path = Plot_folder, plot = p, width = 12, height = 12, units = "cm")
}