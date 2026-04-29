## 2026.03.20
## Compare DE peaks, color dots by MATAC cluster
## Based on 2025_10_13_TCseq_plot_K27ac.R
## Implement for DiffBind
library(ggplot2)
library(dplyr)
library(stringr)
library(bedtoolsr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/DiffBind/")

## -- Run automatically for all comparisons, both marks.

Genotypes <- c("TKO")
#Genotypes <- c("DKO", "TKO")
Days <- c("D0", "D2", "D4")
Marks <- c("K4me1", "K27ac")

gen1 <- "WT"

for(i in 1:length(Genotypes)) {
  gen2 <- Genotypes[i]
  for(j in 1:length(Days)) {
    day1 <- Days[j]
    day2 <- Days[j]
    for(k in 1:length(Marks)) {
      mark <- Marks[k]

      ## -- Initial arguments --
      group1 <- paste0(gen1,"_",day1)
      group2 <- paste0(gen2,"_",day2)
      cat("Analyzing: ",group1,"_vs_",group2," [",mark,"]\n", sep = "")

      ## --------------------------- ##
      ## -- Read DiffBind results -- ##
      ## --------------------------- ##

      ## Specify folder
      if(mark == "K4me1") {
        mark_folder <- "DiffBind_WT-DKO-TKO_K4me1_RmProm/results"
        #} else if (mark == "MATAC") {
        #  mark_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT-DKO-TKO_MATAC/results"
      } else if (mark == "K27ac") {
        mark_folder <- "DiffBind_WT-DKO-TKO_K27ac_i_MATAC_RmProm/results"
      }

      ## Specify file name
      mark_name <- list.files(path = mark_folder, full.names = FALSE)
      mark_name <- mark_name[startsWith(x = mark_name, prefix = paste0("DiffBind_",group1,"_vs_",group2,"_Fullpeaks_DESeq2"))]
      mark_name <- mark_name[!grepl(pattern = "AllPeaks", x = mark_name, ignore.case = FALSE)]
      if(length(mark_name) == 0) {
        cat("No DiffBind-DESeq2 results found for: ",group1,"_vs_",group2," [",mark,"]", sep = "")
        q(save = "ask")
      } else if(length(mark_name) > 1) {
        cat("Multiple DiffBind-DESeq2 results found for: ",group1,"_vs_",group2," [",mark,"]", sep = "")
        q(save = "ask")
      }
      ## If only one file found
      DiffBind_mark <- read.table(file = paste(mark_folder, mark_name, sep = "/"), header = TRUE, sep = '\t')
      ## Add coordinate
      DiffBind_mark$coord <- str_c(str_c(DiffBind_mark$chr, DiffBind_mark$start, sep = ":"), DiffBind_mark$end, sep = "-")
      DiffBind_mark <- DiffBind_mark |> relocate(coord, .after = end)

      ## Identify differentially-enriched
      DiffBind_mark$direction <- "none"
      DiffBind_mark$direction[DiffBind_mark$FDR < 0.05 & DiffBind_mark$logFC > 1] <- "up"
      DiffBind_mark$direction[DiffBind_mark$FDR < 0.05 & DiffBind_mark$logFC < -1] <- "down"
      DiffBind_mark$direction <- factor(DiffBind_mark$direction, levels = c("down", "none", "up"))
      table(DiffBind_mark$direction)

      ## ------------------------- ##
      ## -- Read MATAC clusters -- ##
      ## ------------------------- ##

      peaks_clust <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT_MATAC_v1/Dynamic/K5/bed/All_clusters_labeled.bed", header = FALSE, sep = '\t')
      colnames(peaks_clust) <- c("chr", "start", "end", "cluster")
      static_clust <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT_MATAC_v1/Dynamic/K5/bed/Static_peaks.bed3", header = FALSE, sep = '\t')
      colnames(static_clust) <- c("chr", "start", "end")
      static_clust$cluster <- "static"
      peaks_clust <- rbind(peaks_clust, static_clust)
      rm(static_clust)

      ## -- Annotate DiffBind peaks
      DB_ov <- bt.intersect(a = DiffBind_mark, b = peaks_clust, c = TRUE)
      colnames(DB_ov) <- c(colnames(DiffBind_mark), "num_ov")
      table(DB_ov$num_ov)
      # 0:   2,916
      # 1:  25,455
      # 2:     947
      # 3:      46
      # 4:       4

      ## Treat no-overlaps, single-overlaps, and multiple-overlaps independently.
      DB_NoOv <- DB_ov[DB_ov$num_ov == 0,]
      DB_NoOv <- DB_NoOv[, colnames(DB_NoOv) != "num_ov"]
      DB_SingleOv <- DB_ov[DB_ov$num_ov == 1,]
      DB_SingleOv <- DB_SingleOv[, colnames(DB_SingleOv) != "num_ov"]
      DB_MultiOv <- DB_ov[DB_ov$num_ov > 1,]
      DB_MultiOv <- DB_MultiOv[, colnames(DB_MultiOv) != "num_ov"]
      rm(DB_ov)

      ## Anotate those with overlap
      DB_SingleOv <- bt.intersect(a = DB_SingleOv, b = peaks_clust, wo = TRUE)
      DB_MultiOv <- bt.intersect(a = DB_MultiOv, b = peaks_clust, wo = TRUE)
      colnames(DB_SingleOv) <- c(colnames(DiffBind_mark), "MATAC.chr", "MATAC.start", "MATAC.end", "MATAC.cluster", "bp_ov")
      colnames(DB_MultiOv) <- c(colnames(DiffBind_mark), "MATAC.chr", "MATAC.start", "MATAC.end", "MATAC.cluster", "bp_ov")

      table(table(DB_MultiOv$coord))
      ## Should correspond to the above values

      ## For no overlap, add fake columns.
      DB_NoOv$MATAC.chr <- NA
      DB_NoOv$MATAC.start <- NA
      DB_NoOv$MATAC.end <- NA
      DB_NoOv$MATAC.cluster <- "ATACneg"
      DB_NoOv$bp_ov <- 0

      ## For multi-overlap, iterate over every coordinate and chose the loci with the biggest overlap
      ## Factor cluster for sorting
      DB_MultiOv$MATAC.cluster <- factor(DB_MultiOv$MATAC.cluster, levels = c("c1", "c2", "c3", "c4", "c5", "static"))
      DB_MultiOv_final <- DB_MultiOv[0,]
      Coord_multi_ov <- unique(DB_MultiOv$coord)

      for(i in 1:length(Coord_multi_ov)){
        Sub_coord <- DB_MultiOv[DB_MultiOv$coord == Coord_multi_ov[i],] |>
          arrange(desc(bp_ov), MATAC.cluster) |> slice_head(n = 1)
        DB_MultiOv_final <- rbind(DB_MultiOv_final, Sub_coord)
      }
      rm(i, Sub_coord, Coord_multi_ov)

      ## Combine all the overlaps
      DB_final <- rbind(DB_NoOv, DB_SingleOv, DB_MultiOv_final) |> arrange(chr, start, end)
      ## Consider exporting table

      rm(DB_NoOv, DB_SingleOv, DB_MultiOv, DB_MultiOv_final)

      ## --------------- ##
      ## -- Plot time -- ##
      ## --------------- ##

      ## Prepare table for plotting
      DB_final$direction <- factor(DB_final$direction, levels = c("down", "none", "up"))
      DB_final$MATAC.cluster <- factor(DB_final$MATAC.cluster, levels = c("c1", "c2", "c3", "c4", "c5", "static", "ATACneg"))

      ## Make column for color
      DB_final$Color <- DB_final$MATAC.cluster
      ## Increase factor levels
      DB_final$Color <- forcats::fct_expand(DB_final$Color, "Not_DE")
      ## Change color for those peaks that are not differential
      DB_final$Color[DB_final$direction == "none"] <- "Not_DE"

      # table(DB_final$Color)
      # table(DB_final[,c("MATAC.cluster", "Color")])

      if (gen2 == "DKO") { label2 <- "Tet KO"
      } else if (gen2 == "TKO") { label2 <- "Dnmt KO"
      }

      Name_folder <- gsub(pattern = "/results", replacement = "", x = mark_folder)

      ## -- Color by up/down
      p <- ggplot(DB_final, aes(x = logFC, y = -log10(FDR), color = direction)) +
        geom_point(size = 0.8, alpha = 0.4) +
        scale_color_manual(values = c(none = "gray50", up = "#F8766D", down = "#00BFC4")) +
        geom_hline(yintercept = -log10(0.05), linewidth = 0.2, linetype = "dashed", color = 'gray30') +
        geom_vline(xintercept = c(-1,1), linewidth = 0.2, linetype = "dashed", color = 'gray30') +
        labs(x = paste0("logFC(",label2,"/",gen1,") ",mark)) +
        theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
      p
      Name_file <- paste0(group1,"_vs_",group2,"_",mark,"_DE_volcano.png")
      ggsave(filename = Name_file, device = "png", plot = p, width = 12, height = 12, units = "cm",
             path = paste0("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/DiffBind/",Name_folder,"/"))

      ## Colored by MATAC cluster
      p <- ggplot(DB_final, aes(x = logFC, y = -log10(FDR), color = Color)) +
        geom_point(size = 0.8, alpha = 0.4) +
        scale_color_manual(values = c("#E4E319", "#90D643", "#2D908C", "#5A4791", "#571B67","#95B6D6","gray80","gray20"), drop = FALSE) +
        geom_hline(yintercept = -log10(0.05), linewidth = 0.2, linetype = "dashed", color = 'gray30') +
        geom_vline(xintercept = c(-1,1), linewidth = 0.2, linetype = "dashed", color = 'gray30') +
        labs(x = paste0("logFC(",label2,"/",gen1,") ",mark)) +
        theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
      p
      Name_file <- paste0(group1,"_vs_",group2,"_",mark,"_DE_volcano_byCluster.png")
      ggsave(filename = Name_file, device = "png", plot = p, width = 12, height = 12, units = "cm",
             path = paste0("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/DiffBind/",Name_folder,"/"))

      ## -- Create summary table -- ##
      tot <- nrow(DB_final)
      not_de <- sum(DB_final$direction == "none")
      tot_up <- sum(DB_final$direction == "up")
      tot_down <- sum(DB_final$direction == "down")

      Clust_table <- as.data.frame(table(DB_final$direction, DB_final$MATAC.cluster))
      colnames(Clust_table) <- c("direction", "cluster", "freq")
      ## Relevel direction
      Clust_table$direction <- forcats::fct_relevel(Clust_table$direction, c("down", "up", "none"))
      Clust_table <- arrange(Clust_table, direction, cluster)
      Clust_table$code <- paste(Clust_table$direction, Clust_table$cluster, sep = "_")
      Clust_table <- Clust_table[,c("code", "freq")]

      Sum_table <- data.frame("Category" = c("Total", "Not_DE", "Up", "Down", Clust_table$code),
                              "Value" = c(tot, not_de, tot_up, tot_down, Clust_table$freq))

      write.table(x = Sum_table, file = paste0(Name_folder,"/MATAC_clust_labeling/Category_DEpeaks_",group1,"_vs_",group2,"_",mark,".txt"),
                  append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
    }
  }
}
