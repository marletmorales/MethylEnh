## 09.03.2026
## Approaching the project backwards. Start with enhancers, then go into DiffBind results
## Edit from 2026_03_03_Link_Enhancers-GeneExp_with_TCseq_v3.R
## Edit from 2026_03_04_Inverse_Link_Enhancers-GeneExp_with_TCseq_addDNAme.R
## Edit 2026.03.09: Merge DKO and TKO scripts into one
## Edit 2026.03.13: Add p-value adjustment for pval_greater and pval_less
## Edit 2026.03.26: Export table for Gene Ontologies (middle of "Part 8: Intersect Enhancers with RNAseq data")
library(dplyr)
library(stringr)
library(bedtoolsr)
library(ggplot2)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Set_operations/Enhancers_intersect_v6/WT_Enhancer_peakome/")

## -- Initial arguments ----
group1 <- "WT_D4"
group2 <- "DKO_D4"
mark <- "K27ac"
gen1 <- str_split_1(string = group1, pattern = "_")[1]
day1 <- str_split_1(string = group1, pattern = "_")[2]
gen2 <- str_split_1(string = group2, pattern = "_")[1]
day2 <- str_split_1(string = group2, pattern = "_")[2]
## ----
## Create log file
# Name_log <- paste0("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/MEP_DiffBind_RNAseq/Enhancers_intersect_v6/WT_Enhancer_peakome/",
#                    group1,"_vs_",group2,"__",mark,"_DNAme_RNAseq_log_v1.txt")
#cat("Analyzing: ",group1,"_vs_",group2," [",mark,"]\n\n", sep = "", file = Name_log, append = FALSE)
## ----

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
Annotate_TCseq <- function(WT_data, Mark_label, Mark_data) {
  ## Initialize output
  Final_output <- cbind(WT_data[0,], Mark_data[0,], data.frame("bp_ov" = numeric()))
  ## Do not add WT. to the beginning of the columns (at least for this script)
  colnames(Final_output) <- c(colnames(WT_data),
                              paste(Mark_label, colnames(Mark_data), sep = "."),
                              paste(Mark_label, "bp_ov", sep = "."))
  # colnames(Final_output) <- c(paste("WT", colnames(WT_data), sep = "."),
  #                             paste(Mark_label, colnames(Mark_data), sep = "."),
  #                             paste(Mark_label, "bp_ov", sep = "."))

  ## Extract number of overlaps per WT_data loci
  ovs <- bt.intersect(a = WT_data, b = Mark_data, c = TRUE)
  colnames(ovs) <- c(colnames(WT_data), "num_ov")
  #cat("Intersects: (WT_data ∩ Mark_data)\n")
  print(table(ovs$num_ov))

  ## Process no overlaps
  res_no_ov <- ovs[ovs$num_ov == 0,]
  if(nrow(res_no_ov) > 0) {
    res_no_ov <- res_no_ov[, colnames(res_no_ov) != "num_ov"]
    #colnames(res_no_ov) <- paste("WT", colnames(res_no_ov), sep = ".")
    ## Append empty columns
    tmp <- matrix(data = NA, nrow = nrow(res_no_ov), ncol = ncol(Mark_data)+1)
    colnames(tmp) <- c(paste(Mark_label, colnames(Mark_data), sep = "."),
                       paste(Mark_label, "bp_ov", sep = "."))
    tmp <- as.data.frame(tmp)
    res_no_ov <- cbind(res_no_ov, tmp)
    ## Specify direction as "none"
    res_no_ov[, paste(Mark_label, "direction", sep = ".")] <- "none"
    ## Bind to final output
    Final_output <- rbind(Final_output, res_no_ov)
  }

  ## Process single overlaps
  res_single_ov <- ovs[ovs$num_ov == 1,]
  if(nrow(res_single_ov) > 0) {
    res_single_ov <- res_single_ov[, colnames(res_single_ov) != "num_ov"]
    res_single_ov <- bt.intersect(a = res_single_ov, b = Mark_data, wo = TRUE)
    ## Do not add WT. to the beginning of the columns (at least for this script)
    colnames(res_single_ov) <- c(colnames(WT_data),
                                 paste(Mark_label, colnames(Mark_data), sep = "."),
                                 paste(Mark_label, "bp_ov", sep = "."))
    # colnames(res_single_ov) <- c(paste("WT", colnames(WT_data), sep = "."),
    #                              paste(Mark_label, colnames(Mark_data), sep = "."),
    #                              paste(Mark_label, "bp_ov", sep = "."))
    ## Bind to final output
    Final_output <- rbind(Final_output, res_single_ov)
  }

  ## Process multi overlaps
  res_mult_ov <- ovs[ovs$num_ov > 1,]
  if(nrow(res_mult_ov) > 0) {
    res_mult_ov <- res_mult_ov[, colnames(res_mult_ov) != "num_ov"]
    res_mult_ov <- bt.intersect(a = res_mult_ov, b = Mark_data, wo = TRUE)
    ## Do not add WT. to the beginning of the columns (at least for this script)
    colnames(res_mult_ov) <- c(colnames(WT_data),
                               paste(Mark_label, colnames(Mark_data), sep = "."), "bp_ov")
    # colnames(res_mult_ov) <- c(paste("WT", colnames(WT_data), sep = "."),
    #                            paste(Mark_label, colnames(Mark_data), sep = "."), "bp_ov")
    ## Iterate over every coordinate and chose the loci with the biggest overlap
    Coord_mult_ov <- unique(res_mult_ov$coord)
    #Coord_mult_ov <- unique(res_mult_ov$WT.coord)
    res_mult_ov_single <- res_mult_ov[0,]
    for(i in 1:length(Coord_mult_ov)){
      Sub_coord <- res_mult_ov[res_mult_ov$coord == Coord_mult_ov[i],]
      #Sub_coord <- res_mult_ov[res_mult_ov$WT.coord == Coord_mult_ov[i],]
      max_bp <- max(Sub_coord$bp_ov)
      Sub_coord <- Sub_coord[Sub_coord$bp_ov == max_bp,]
      Sub_coord <- Sub_coord[1,] # in case mutliple loci have the max bp
      res_mult_ov_single <- rbind(res_mult_ov_single, Sub_coord)
    }
    ## Bind to final output
    ## Do not add WT. to the beginning of the columns (at least for this script)
    colnames(res_mult_ov_single) <- c(colnames(WT_data),
                                      paste(Mark_label, colnames(Mark_data), sep = "."),
                                      paste(Mark_label, "bp_ov", sep = "."))
    # colnames(res_mult_ov_single) <- c(paste("WT", colnames(WT_data), sep = "."),
    #                                   paste(Mark_label, colnames(Mark_data), sep = "."),
    #                                   paste(Mark_label, "bp_ov", sep = "."))
    Final_output <- rbind(Final_output, res_mult_ov_single)
  }
  ## Clean and return
  Final_output <- arrange(Final_output, chr, start, end)
  #Final_output <- arrange(Final_output, WT.chr, WT.start, WT.end)
  return(Final_output)
}
## ----

## -- Part 1: Read Enhancer-GeneExp data ----
#cat("== Reading enhancers ==\n\n", sep = "", file = Name_log, append = TRUE)

## Read original, then drop gene symbol and counts
Enh_GeneExp_orig <- read.table(file = "WT_Master_Enhancer_Peakome_v6_annotated_TargetGene_wExpression_wActvM.txt", header = TRUE, sep = '\t', na.strings = "")
Enh_GeneExp_orig$coord <- str_c(str_c(Enh_GeneExp_orig$chr, Enh_GeneExp_orig$start, sep = ":"), Enh_GeneExp_orig$end, sep = "-")
Enh_GeneExp_orig <- Enh_GeneExp_orig |> relocate(coord, .after = end)

## Drop gene symbols and counts
Enh <- Enh_GeneExp_orig[,c("chr", "start", "end", "coord", "D0_class", "D2_class", "D4_class", "peak_id", "D0_class_ext", "D2_class_ext", "D4_class_ext")] |> distinct()

## -- Factorize enhancer classes
## 3-class
Enh[,c("D0_class", "D2_class", "D4_class")] <- Factorize_enhancers(Enh[,c("D0_class", "D2_class", "D4_class")])
## 4-class
Enh[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")] <- Factorize_enhancers(Enh[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")])

## -- Keep only those enhancers that are not 'othr'
#n_orig <- length(unique(Enh$peak_id))
index1 <- Enh[,paste(day1,"_class",sep = "")] != "othr"
index2 <- Enh[,paste(day2,"_class",sep = "")] != "othr"
Enh <- Enh[index1|index2,]
#n_tot <- length(unique(Enh$peak_id))
#cat("Number of total enhancers: ",n_orig,"\n", sep = "", file = Name_log, append = TRUE)
#cat("Number of enhancers considered for ",group1,"_vs_",group2," comparison: ",n_tot,"\n", sep = "", file = Name_log, append = TRUE)
## Print enhancer classification
# tmp <- Enh[,c("chr","start","end","coord","peak_id",paste0(day1,"_class_ext"))] |> distinct()
# tmp <- table(tmp[,paste0(day1,"_class_ext")], exclude = "othr")
# tmp <- data.frame("Class" = names(tmp),
#                   "Value" = as.vector(tmp))
#write.table(x = tmp, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
#cat("\n", sep = "", file = Name_log, append = TRUE)
rm(index1,index2, tmp)

## ------------------------------------------
## ------------------------------------------

## -- Part 2: Import information on DNAme ----
DNAme_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/Enhancers_intersect_v3"
DNAme_files <- list.files(path = DNAme_folder, full.names = FALSE, pattern = "meanDNAme_.+__Coordinates_WT_Master_Enhancer_Peakome_v3.*\\.txt")
Group1_files <- DNAme_files[grepl(pattern = paste0("meanDNAme_",group1,"_"), x = DNAme_files)]
Group2_files <- DNAme_files[grepl(pattern = paste0("meanDNAme_",group2,"_"), x = DNAme_files)]

## Make empty table to append information to
Final_DNAme <- read.table(file = paste0(DNAme_folder,"/",Group1_files[1]), header = TRUE, sep = "\t")[0,]
Final_DNAme <- cbind(Final_DNAme, data.frame("gen" = character(), "mark" = character()))

## Read files from group1 (should be WT)
if(gen1 == "WT") {
  if(length(Group1_files) > 0) {
    for(i in 1:length(Group1_files)) {
      tmp <- read.table(file = paste0(DNAme_folder,"/",Group1_files[i]), header = TRUE, sep = "\t")
      tmp$gen <- group1
      tmp$mark <- case_when(grepl(pattern = "K4me1", x = Group1_files[i]) ~ "K4me1",
                            grepl(pattern = "K27ac", x = Group1_files[i]) ~ "K27ac",
                            grepl(pattern = "MATAC", x = Group1_files[i]) ~ "MATAC",
                            grepl(pattern = "WGBS", x = Group1_files[i]) ~ "WGBS")
      Final_DNAme <- rbind(Final_DNAme, tmp)
    }
    rm(tmp)
  }
} else {
  #cat("Group 1 is different from WT; DNAme comparisons may be innaccurate.\n")
  #cat("Group 1 is different from WT; DNAme comparisons may be innaccurate.\n", file = Name_log, append = TRUE)
}

## Read files from group2 (DKO or TKO)
if(gen2 != "TKO") {
  if(length(Group2_files) > 0) {
    for(i in 1:length(Group2_files)) {
      tmp <- read.table(file = paste0(DNAme_folder,"/",Group2_files[i]), header = TRUE, sep = "\t")
      tmp$gen <- group2
      tmp$mark <- case_when(grepl(pattern = "K4me1", x = Group2_files[i]) ~ "K4me1",
                            grepl(pattern = "K27ac", x = Group2_files[i]) ~ "K27ac",
                            grepl(pattern = "MATAC", x = Group2_files[i]) ~ "MATAC",
                            grepl(pattern = "WGBS", x = Group2_files[i]) ~ "WGBS")
      Final_DNAme <- rbind(Final_DNAme, tmp)
    }
    rm(tmp)
  }
} else {
  #cat("Group 2 is equal to TKO; skipping DNAme analysis for this genotype (using only Group1).\n\n")
  #cat("Group 2 is equal to TKO; skipping DNAme analysis for this genotype (using only Group1).\n\n", file = Name_log, append = TRUE)
}

## Add coordinate to DNAme levels
Final_DNAme$coord <- str_c(str_c(Final_DNAme$chr, Final_DNAme$start, sep = ":"), Final_DNAme$end, sep = "-")
Final_DNAme <- Final_DNAme |> relocate(coord, .after = end)
rm(DNAme_files, DNAme_folder, Group1_files, Group2_files, i)

## ------------------------------------------
## ------------------------------------------

## -- Part 3: Add DNAme information to Enhancers ----
Final_DNAme <- Final_DNAme |> filter(coord %in% Enh$coord)
Final_DNAme <- Final_DNAme |> tidyr::pivot_wider(names_from = c(gen, mark),
                                                 values_from = num.cpgs:meth.wmean,
                                                 names_glue = "{gen}_{mark}_{.value}")
## Append to Enhancers
Enh <- merge(x = Enh, y = Final_DNAme, by = c("chr", "start", "end", "coord"), all.x = TRUE, all.y = FALSE)
rm(Final_DNAme)

## ------------------------------------------
## ------------------------------------------

## -- Part 4: Read mark-specific DiffBind results ----
#cat("\n== Reading chrom mark results ==\n\n", sep = "", file = Name_log, append = TRUE)

## Specify folder
if(mark == "K4me1") {
  mark_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/DiffBind/DiffBind_WT-DKO-TKO_K4me1_RmProm/results"
  #} else if (mark == "MATAC") {
  #  mark_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/TCseq/WT-DKO-TKO_MATAC/results"
} else if (mark == "K27ac") {
  mark_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/DiffBind/DiffBind_WT-DKO-TKO_K27ac_i_MATAC_RmProm/results"
}

## Specify file name
mark_name <- list.files(path = mark_folder, full.names = FALSE)
mark_name <- mark_name[startsWith(x = mark_name, prefix = paste0("DiffBind_",group1,"_vs_",group2,"_Fullpeaks_DESeq2"))]
mark_name <- mark_name[!grepl(pattern = "AllPeaks", x = mark_name, ignore.case = FALSE)]
if(length(mark_name) == 0) {
  #cat("No DiffBind-DESeq2 results found for: ",group1,"_vs_",group2," [",mark,"]\nExiting", sep = "")
  q(save = "no")
} else if(length(mark_name) > 1) {
  #cat("Multiple DiffBind-DESeq2 results found for: ",group1,"_vs_",group2," [",mark,"]\nExiting", sep = "")
  q(save = "no")
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
#DiffBind_mark_diff <- DiffBind_mark |> filter(direction != "none")

## -- Stats per dataset
#n_tot <- nrow(DiffBind_mark)
#n_dif <- nrow(DiffBind_mark_diff)
#cat("Number of tested ",mark," regions: ",n_tot,"\n", sep = "", file = Name_log, append = TRUE)
#cat("Number of differential ",mark," regions: ",n_dif," (",round(n_dif/n_tot*100, 2)," %)\n", sep = "", file = Name_log, append = TRUE)
#cat("Down:\t",sum(DiffBind_mark_diff$direction == "down"),"\n", sep = "", file = Name_log, append = TRUE)
#cat("Up:\t",sum(DiffBind_mark_diff$direction == "up"),"\n\n", sep = "", file = Name_log, append = TRUE)

## How many of the total assessed regions intersect with enhancers?
#cat("\n== Intersecting enhancers with chrom mark ==\n\n", sep = "", file = Name_log, append = TRUE)
#n_int <- nrow(bt.intersect(a = DiffBind_mark, b = Enh, u = TRUE))
#cat("Number of tested ",mark," that intersects with enhancers: ",n_int," (",round(n_int/n_tot*100, 2)," %)\n", sep = "", file = Name_log, append = TRUE)

## How many of the DE regions intersect with enhancers?
#n_int_diff <- nrow(bt.intersect(a = DiffBind_mark_diff, b = Enh, u = TRUE))
#cat("Number of differential ",mark," that intersects with enhancers: ",n_int_diff," (",
#    round(n_int_diff/n_dif*100, 2)," % of differential; ",
#    round(n_int_diff/n_int*100, 2)," % of Intersect)\n", sep = "", file = Name_log, append = TRUE)
## Direction of the DE_intersect
#table_int_diff <- bt.intersect(a = DiffBind_mark_diff, b = Enh, u = TRUE)
#cat("Down:\t",sum(table_int_diff[,ncol(table_int_diff)] == "down"),"\n", sep = "", file = Name_log, append = TRUE)
#cat("Up:\t",sum(table_int_diff[,ncol(table_int_diff)] == "up"),"\n\n", sep = "", file = Name_log, append = TRUE)

## Of the X tested regions that intersect with enhancers, how many enhancers are represented?
#tmp_Enh <- Enh[,c("chr","start","end","peak_id")] |> distinct()
#n_int_rev <- nrow(bt.intersect(a = tmp_Enh, b = DiffBind_mark, u = TRUE))
#cat("Number of enhancers that intersect with tested ",mark,": ",n_int_rev," (",round(n_int_rev/nrow(tmp_Enh)*100, 2)," %)\n", sep = "", file = Name_log, append = TRUE)

## How many enhancers intersect with DE regions?
#n_int_diff_rev <- nrow(bt.intersect(a = tmp_Enh, b = DiffBind_mark_diff, u = TRUE))
#cat("Number of enhancers that intersect with differential ",mark,": ",n_int_diff_rev," (",round(n_int_diff_rev/n_int_rev*100, 2)," % of Intersect)\n\n", sep = "", file = Name_log, append = TRUE)

rm(table_int_diff, tmp_Enh, n_dif, n_int, n_orig, n_int_diff, n_int_diff_rev, n_int_rev, n_tot, mark_folder, mark_name, DiffBind_mark_diff)

## ------------------------------------------
## ------------------------------------------

## -- Part 5: Merge Enhancers and DiffBind results ----
Enh_int <- bt.intersect(a = Enh, b = DiffBind_mark, u = TRUE)
colnames(Enh_int) <- colnames(Enh)

## -- Factorize enhancer classes
## 3-class
Enh_int[,c("D0_class", "D2_class", "D4_class")] <- Factorize_enhancers(Enh_int[,c("D0_class", "D2_class", "D4_class")])
## 4-class
Enh_int[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")] <- Factorize_enhancers(Enh_int[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")])

## -- Stats
#cat("Classification of intersected enhancers:\n", sep = "", file = Name_log, append = TRUE)
# tmp <- Enh_int[,c("chr","start","end","coord","peak_id",paste0(day1,"_class_ext"))] |> distinct()
# tmp <- table(tmp[,paste0(day1,"_class_ext")], exclude = "othr")
# tmp <- data.frame("Class" = names(tmp),
#                   "Value" = as.vector(tmp))
#write.table(x = tmp, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
#cat("\n", sep = "", file = Name_log, append = TRUE)

## -- Merge intersected enhancers with DiffBind metrics
Enh_int <- Annotate_TCseq(WT_data = Enh_int, Mark_label = "DiffBind", Mark_data = DiffBind_mark)
## Factorize enhancer classes
Enh_int[,c("D0_class", "D2_class", "D4_class")] <- Factorize_enhancers(Enh_int[,c("D0_class", "D2_class", "D4_class")])
Enh_int[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")] <- Factorize_enhancers(Enh_int[,c("D0_class_ext", "D2_class_ext", "D4_class_ext")])
## Factorize direction
Enh_int$DiffBind.direction <- factor(Enh_int$DiffBind.direction, levels = c("down", "none", "up"))

## -- Export number of DE observations per enhancer class
#cat("Number and direction of DE enhancers:\n", sep = "", file = Name_log, append = TRUE)
# tmp_class <- Enh_int |> group_by(.data[[paste0(day1,"_class_ext")]]) |> summarise("down" = sum(DiffBind.direction == "down"),
#                                                                                   "none" = sum(DiffBind.direction == "none"),
#                                                                                   "up" = sum(DiffBind.direction == "up"),
#                                                                                   "Total" = n())
# tmp_Tot <- Enh_int |> summarise("down" = sum(DiffBind.direction == "down"),
#                                 "none" = sum(DiffBind.direction == "none"),
#                                 "up" = sum(DiffBind.direction == "up"),
#                                 "Total" = n())
#cat("\tdown\tnone\tup\tTotal\n", sep = "", file = Name_log, append = TRUE)
#write.table(x = tmp_class, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
#cat("----------------------------------------\n", sep = "", file = Name_log, append = TRUE)
#cat("Total\t", sep = "", file = Name_log, append = TRUE)
#write.table(x = tmp_Tot, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
#cat("\n", sep = "", file = Name_log, append = TRUE)

## Clear workspace
rm(Enh, DiffBind_mark, tmp, tmp_class, tmp_Tot)

## ------------------------------------------
## ------------------------------------------

## -- Part 6: Plots Enhancers <=> DiffBind ----

## Specify plots folder
#Plot_folder <- paste0("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Plots/MEP_DiffBind_RNAseq/Enhancers_intersect_v6/WT_Enhancer_peakome/",group1,"_vs_",group2,"_",mark,"_v2/")
gen_colors <- c("WT" = "#266799",
                "DKO" = "#4C4C4C",
                "TKO" = "#F68532")

## -- Distribution of DiffBind.logFC by enhancer class
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_distribution.png")
p <- ggplot(Enh_int, aes(x = DiffBind.logFC, fill = .data[[paste0(day1, "_class_ext")]])) +
  geom_density(alpha = 0.5) +
  facet_wrap({{paste0(day1,"_class_ext")}}, ncol = 1) +
  scale_fill_manual(values = c("#E5614E", "#97C124", "#F0AC4B", "#77BEF0")) +
  geom_vline(xintercept = c(-1,0,1), color = "gray20", linewidth = 0.25, linetype = "dashed") +
  labs(x = paste0("log2FC(",mark,")")) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 10, height = 13, units = "cm", create.dir = TRUE)
rm(p)

## -- Differential Enrichment: DiffBind.logFC vs FDR; facet by enhancer class
## All enhancer classes
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_FDR_byEnhClassExt_full.png")
p <- ggplot(Enh_int, aes(x = DiffBind.logFC, y = -log10(DiffBind.FDR), color = DiffBind.direction)) +
  geom_point(na.rm = TRUE, size = 0.85, alpha = 0.7) +
  scale_color_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
  labs(x = paste0("log2FC(",mark,")"), y = "-log10(FDR)") +
  geom_hline(yintercept = -log10(0.05), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}} )
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
rm(p)

## Only actv, actvM, and prim enhancers
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_FDR_byEnhClassExt_aap.png")
p <- ggplot(Enh_int[Enh_int[,paste0(day1,"_class_ext")] %in% c("actv", "actvM", "prim"),],
            aes(x = DiffBind.logFC, y = -log10(DiffBind.FDR), color = DiffBind.direction)) +
  geom_point(na.rm = TRUE, size = 0.85, alpha = 0.7) +
  scale_color_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
  labs(x = paste0("log2FC(",mark,")"), y = "-log10(FDR)") +
  geom_hline(yintercept = -log10(0.05), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 30, height = 10, units = "cm", create.dir = TRUE)
rm(p)

## --------------------------------------------------------------------- ##
## -- For paper: Only actv and prim enhancers (no actvM) when day==D0 -- ##
## --------------------------------------------------------------------- ##
if(day1 == "D0"){
  p <- ggplot(Enh_int[Enh_int[,paste0(day1,"_class_ext")] %in% c("actv", "prim"),],
              aes(x = DiffBind.logFC, y = -log10(DiffBind.FDR), color = DiffBind.direction)) +
    geom_point(na.rm = TRUE, size = 0.85, alpha = 0.7) +
    scale_color_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
    labs(x = paste0("log2FC(",mark,")"), y = "-log10(FDR)") +
    geom_hline(yintercept = -log10(0.05), color = "gray0", linewidth = 0.2, linetype = "dashed") +
    geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
  p
  Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_FDR_byEnhClassExt_ap.png")
  # ggsave(filename = Name_plot, plot = p, device = "png", width = 20, height = 10, units = "cm", create.dir = FALSE,
  #        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_4/Supp-X_Volcano_plot_K27ac_DE/")
  rm(p)
}

## ------------------------------------------
## ------------------------------------------

## -- Part 7: Analyse DNAme ----
#cat("\n== Adding DNAme information ==\n\n", sep = "", file = Name_log, append = TRUE)

## Make a copy of Enh_int for pivoting
Enh_int_pivot <- Enh_int |> tidyr::pivot_longer(cols = contains("_meth.wmean"),
                                                names_to = c("Group", "Mark"),
                                                values_to = "meth.wmean",
                                                names_pattern = "^(.+)_([^_]+)_meth\\.wmean$")
## Drop missing meth.wmean
Enh_int_pivot <- Enh_int_pivot |> filter(!is.na(meth.wmean))
## Factorize columns
Enh_int_pivot$Group <- factor(Enh_int_pivot$Group, levels = c(group1, group2))
Enh_int_pivot$Mark <- factor(Enh_int_pivot$Mark, levels = c("K4me1", "MATAC", "K27ac", "WGBS"))

## -- Plot DNAme for each genotype, for each mark, faceted by enhancer class
## All enhancer classes
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAme_byMark_byEnhClassExt_full.png")
p <- ggplot(Enh_int_pivot, aes(x = DiffBind.direction, y = meth.wmean, fill = Group)) +
  geom_boxplot(alpha = 0.9, outliers = FALSE) +
  scale_fill_manual(values = as.vector(gen_colors[c(gen1,gen2)])) +
  labs(x = paste0(mark," DE"), y = "DNAme %") +
  coord_cartesian(ylim = c(0,100)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
rm(p)

## Only actv, actvM, and prim enhancers; no K27ac 5mC
tmp_aap <- Enh_int_pivot |> filter(Mark != "K27ac" & .data[[paste0(day1, "_class_ext")]] %in% c("actv", "actvM", "prim"))
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAme_byMark_byEnhClassExt_aap.png")
p <- ggplot(tmp_aap, aes(x = DiffBind.direction, y = meth.wmean, fill = Group)) +
  geom_boxplot(alpha = 0.9, outliers = FALSE) +
  scale_fill_manual(values = as.vector(gen_colors[c(gen1,gen2)])) +
  labs(x = paste0(mark," DE"), y = "DNAme %") +
  coord_cartesian(ylim = c(0,100)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 16, units = "cm", create.dir = TRUE)
rm(p, tmp_aap)
## ----

## -- Stats for DNAme WT ----
## Only useful if comparing WT vs TKO

# if(gen1=="WT" & gen2=="TKO") {
#   Enh_classes <- c("actv", "actvM", "prim")
#   Mark_classes <- c("K4me1", "MATAC", "K27ac", "WGBS")
#   t_Output <- data.frame("Enh_class" = character(),
#                          "Mark" = character(),
#                          "GroupA" = character(),
#                          "GroupB" = character(),
#                          "A_n" = numeric(),
#                          "B_n" = numeric(),
#                          "A_mean_5mC" = numeric(),
#                          "B_mean_5mC" = numeric(),
#                          "A_median_5mC" = numeric(),
#                          "B_median_5mC" = numeric(),
#                          "pval" = numeric())
#   for(i in 1:length(Enh_classes)) {
#     for(j in 1:length(Mark_classes)) {
#       ## Pull DNAme values
#       Down <- Enh_int_pivot |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                         Mark == Mark_classes[j] & DiffBind.direction == "down" & Group == group1) |> pull(meth.wmean)
#       None <- Enh_int_pivot |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                         Mark == Mark_classes[j] & DiffBind.direction == "none" & Group == group1) |> pull(meth.wmean)
#       Up <- Enh_int_pivot |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                       Mark == Mark_classes[j] & DiffBind.direction == "up" & Group == group1) |> pull(meth.wmean)
#
#       ## Test comparisons -- minimum sample size 6
#       if(min(length(Down), length(None), length(Up)) >= 6) {
#         t_DN <- wilcox.test(x = Down, y = None)
#         t_UN <- wilcox.test(x = Up, y = None)
#         t_DU <- wilcox.test(x = Down, y = Up)
#
#         ## Append information to output
#         t_Output <- rbind(t_Output, data.frame("Enh_class" = Enh_classes[i],
#                                                "Mark" = Mark_classes[j],
#                                                "GroupA" = c("down", "up", "down"),
#                                                "GroupB" = c("none", "none", "up"),
#                                                "A_n" = c(length(Down), length(Up), length(Down)),
#                                                "B_n" = c(length(None), length(None), length(Up)),
#                                                "A_mean_5mC" = round(c(mean(Down), mean(Up), mean(Down)), digits = 2),
#                                                "B_mean_5mC" = round(c(mean(None), mean(None), mean(Up)), digits = 2),
#                                                "A_median_5mC" = round(c(median(Down), median(Up), median(Down)), digits = 2),
#                                                "B_median_5mC" = round(c(median(None), median(None), median(Up)), digits = 2),
#                                                "pval" = c(t_DN$p.value, t_UN$p.value, t_DU$p.value)))
#       }
#     }
#   }
#   ## Clean output
#   t_Output$Enh_class <- factor(t_Output$Enh_class, levels = Enh_classes)
#   t_Output$Mark <- factor(t_Output$Mark, levels = Mark_classes)
#   t_Output$GroupA <- factor(t_Output$GroupA, levels = c("down", "none", "up"))
#   t_Output$GroupB <- factor(t_Output$GroupB, levels = c("down", "none", "up"))
#   t_Output <- t_Output |> arrange(Enh_class, GroupA, GroupB, Mark)
#   t_Output <- t_Output |> relocate(Mark, .after = GroupB)
#   t_Output$paj <- p.adjust(t_Output$pval, method = "bonferroni")
#   t_Output$sig <- t_Output$paj < 0.05
#   ## Round to four significant decimals
#   t_Output$pval <- signif(t_Output$pval, digits = 5)
#   t_Output$paj <- signif(t_Output$paj, digits = 5)
#
#   #cat("Calculating two-sided wilcox test (5mC ",group1,")\n\n", sep = "", file = Name_log, append = TRUE)
#   #write.table(x = t_Output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
#   #cat("\n", sep = "", file = Name_log, append = TRUE)
#   ## Clean workspace
#   rm(Down, None, Up, t_DN, t_UN, t_DU, t_Output, i, j)
# }
## ----

## -- Stats for DNAme delta ----
## Only useful if working with two DNAme datasets (WT and/or DKO)

if(gen1!="TKO" & gen2!="TKO") {
  ## Create column for easier operation
  Enh_int_pivot <- Enh_int_pivot |> tidyr::pivot_wider(names_from = Group,
                                                       names_glue = "{Group}_meth.wmean",
                                                       values_from = meth.wmean)
  Enh_int_pivot <- Enh_int_pivot |> mutate("DNAme_delta" = .data[[paste0(group2,"_meth.wmean")]] - .data[[paste0(group1,"_meth.wmean")]])

  ## -- Plot d5mC, for each mark, faceted by enhancer class
  ## All enhancer classes
  #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_byMark_byEnhClassExt_full.png")
  p <- ggplot(Enh_int_pivot, aes(x = DiffBind.direction, y = DNAme_delta, fill = DiffBind.direction)) +
    geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
    geom_boxplot(outliers = FALSE, na.rm = TRUE) +
    scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
    labs(x = paste0(mark," DE"), y = paste0("DNAme delta (",group2," - ",group1,")")) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
  p
  #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
  rm(p)

  ## Only actv, actvM, and prim enhancers; no K27ac 5mC
  tmp_aap <- Enh_int_pivot |> filter(Mark != "K27ac" & .data[[paste0(day1, "_class_ext")]] %in% c("actv", "actvM", "prim"))
  #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_byMark_byEnhClassExt_aap.png")
  p <- ggplot(tmp_aap, aes(x = DiffBind.direction, y = DNAme_delta, fill = DiffBind.direction)) +
    geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
    geom_boxplot(outliers = FALSE, na.rm = TRUE) +
    scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
    labs(x = paste0(mark," DE"), y = paste0("DNAme delta (",group2," - ",group1,")")) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
  p
  #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 16, units = "cm", create.dir = TRUE)
  rm(p, tmp_aap)

  ## -- For paper: Only MATAC
  if(day1 %in% c("D2", "D4")){
    ## Only actv, actvM, and prim enhancers
    tmp_aap <- Enh_int_pivot |> filter(Mark == "MATAC" & .data[[paste0(day1, "_class_ext")]] %in% c("actv", "actvM", "prim"))
    p <- ggplot(tmp_aap, aes(x = DiffBind.direction, y = DNAme_delta, fill = DiffBind.direction)) +
      geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
      geom_boxplot(outliers = FALSE, na.rm = TRUE) +
      scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
      labs(x = paste0(mark," DE"), y = "d5mC %") +
      theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
      facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
    p
    Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_MATAC_byEnhClassExt_aap.pdf")
    # ggsave(filename = Name_plot, plot = p, device = "pdf", width = 15, height = 6, units = "cm", create.dir = FALSE,
    #        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_4/Supp-X_Volcano_plot_K27ac_DE/")
    rm(p, tmp_aap)
  } else if(day1 == "D0") {
    ## Only actv and prim enhancers
    tmp_aap <- Enh_int_pivot |> filter(Mark == "MATAC" & .data[[paste0(day1, "_class_ext")]] %in% c("actv", "prim"))
    p <- ggplot(tmp_aap, aes(x = DiffBind.direction, y = DNAme_delta, fill = DiffBind.direction)) +
      geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
      geom_boxplot(outliers = FALSE, na.rm = TRUE) +
      scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
      labs(x = paste0(mark," DE"), y = "d5mC %") +
      theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
      facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
    p
    Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_MATAC_byEnhClassExt_ap.pdf")
    # ggsave(filename = Name_plot, plot = p, device = "pdf", width = 11, height = 6, units = "cm", create.dir = FALSE,
    #        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_4/Supp-X_Volcano_plot_K27ac_DE/")
    rm(p, tmp_aap)
  }

  # ## Export stats for d5mC
  # Enh_classes <- c("actv", "actvM", "prim")
  # Mark_classes <- c("K4me1", "MATAC", "K27ac", "WGBS")
  # t_Output <- data.frame("Enh_class" = character(),
  #                        "Mark" = character(),
  #                        "GroupA" = character(),
  #                        "GroupB" = character(),
  #                        "A_n" = numeric(),
  #                        "B_n" = numeric(),
  #                        "A_mean_d5mC" = numeric(),
  #                        "B_mean_d5mC" = numeric(),
  #                        "A_median_d5mC" = numeric(),
  #                        "B_median_d5mC" = numeric(),
  #                        "pval" = numeric())
  # for(i in 1:length(Enh_classes)) {
  #   for(j in 1:length(Mark_classes)) {
  #     ## Pull DNAme_delta values
  #     Down <- Enh_int_pivot |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
  #                                       Mark == Mark_classes[j] & DiffBind.direction == "down" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
  #     None <- Enh_int_pivot |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
  #                                       Mark == Mark_classes[j] & DiffBind.direction == "none" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
  #     Up <- Enh_int_pivot |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
  #                                     Mark == Mark_classes[j] & DiffBind.direction == "up" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
  #
  #     ## Test comparisons -- minimum sample size 6
  #     if(min(length(Down), length(None), length(Up)) >= 6) {
  #       t_DN <- wilcox.test(x = Down, y = None)
  #       t_UN <- wilcox.test(x = Up, y = None)
  #       t_DU <- wilcox.test(x = Down, y = Up)
  #
  #       ## Append information to output
  #       t_Output <- rbind(t_Output, data.frame("Enh_class" = Enh_classes[i],
  #                                              "Mark" = Mark_classes[j],
  #                                              "GroupA" = c("down", "up", "down"),
  #                                              "GroupB" = c("none", "none", "up"),
  #                                              "A_n" = c(length(Down), length(Up), length(Down)),
  #                                              "B_n" = c(length(None), length(None), length(Up)),
  #                                              "A_mean_d5mC" = round(c(mean(Down), mean(Up), mean(Down)), digits = 2),
  #                                              "B_mean_d5mC" = round(c(mean(None), mean(None), mean(Up)), digits = 2),
  #                                              "A_median_d5mC" = round(c(median(Down), median(Up), median(Down)), digits = 2),
  #                                              "B_median_d5mC" = round(c(median(None), median(None), median(Up)), digits = 2),
  #                                              "pval" = c(t_DN$p.value, t_UN$p.value, t_DU$p.value)))
  #     }
  #   }
  # }
  # ## Clean output
  # t_Output$Enh_class <- factor(t_Output$Enh_class, levels = Enh_classes)
  # t_Output$Mark <- factor(t_Output$Mark, levels = Mark_classes)
  # t_Output$GroupA <- factor(t_Output$GroupA, levels = c("down", "none", "up"))
  # t_Output$GroupB <- factor(t_Output$GroupB, levels = c("down", "none", "up"))
  # t_Output <- t_Output |> arrange(Enh_class, GroupA, GroupB, Mark)
  # t_Output <- t_Output |> relocate(Mark, .after = GroupB)
  # t_Output$paj <- p.adjust(t_Output$pval, method = "bonferroni")
  # t_Output$sig <- t_Output$paj < 0.05
  # ## Round to four significant decimals
  # t_Output$pval <- signif(t_Output$pval, digits = 5)
  # t_Output$paj <- signif(t_Output$paj, digits = 5)
  #
  # #cat("Calculating two-sided wilcox test (d5mC ",group2,"-⁠",group1,")\n\n", sep = "", file = Name_log, append = TRUE)
  # #write.table(x = t_Output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
  # #cat("\n", sep = "", file = Name_log, append = TRUE)
  # ## Clean workspace
  # rm(Down, None, Up, t_DN, t_UN, t_DU, t_Output, i, j)
}
rm(Enh_int_pivot)
## ----

## ------------------------------------------
## ------------------------------------------

## -- Part 8: Intersect Enhancers with RNAseq data ----
#cat("\n== Reading RNAseq information ==\n\n", sep = "", file = Name_log, append = TRUE)

## Retrieve genes symbols and counts from original df
Enh_int <- merge(x = Enh_int, y = Enh_GeneExp_orig,
                 by = intersect(colnames(Enh_int), colnames(Enh_GeneExp_orig)),
                 all.x = TRUE, all.y = FALSE)
## Re-sort by chromosomal position
Enh_int <- Enh_int |> arrange(chr, start, end)
rm(Enh_GeneExp_orig)

## Read DESeq2 output to get p-values
RNAseq_name <- list.files(path = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/htseq_count/RNAseq_v2/results/", full.names = FALSE)
RNAseq_name <- RNAseq_name[startsWith(x = RNAseq_name, prefix = paste(group1,"_vs_",group2,sep = ""))]

if(length(RNAseq_name) == 0) {
  #cat("No RNAseq DESeq2 results found for: ",group1,"_vs_",group2,"\nExiting", sep = "")
  q(save = "no")
} else if(length(RNAseq_name) > 1) {
  #cat("Multiple RNAseq DESeq2 results found for: ",group1,"_vs_",group2,"\nExiting", sep = "")
  q(save = "no")
}
## If only one file found, read
RNAseq_res <- read.table(file = paste("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/htseq_count/RNAseq_v2/results/", RNAseq_name, sep = ""), header = TRUE, sep = '\t')
## Clean table
colnames(RNAseq_res) <- c("gene_symbol", "RNA.baseMean", "RNA.log2FoldChange", "RNA.lfcSE", "RNA.stat", "RNA.pvalue", "RNA.padj", "RNA.direction")

## Add results to table
Enh_int <- merge(x = Enh_int, y = RNAseq_res, by = "gene_symbol", all.x = TRUE, all.y = FALSE)
Enh_int$RNA.direction <- factor(Enh_int$RNA.direction, levels = c("down", "none", "up"))
Enh_int <- Enh_int |> relocate(gene_symbol, .before = dist_sense)
Enh_int <- Enh_int |> arrange(chr, start, end)


## ~~ Little Break
## ---------------------------------------------------------------------------- ##
## -- For Gene Ontologies: Export genes associated to enhancers with DE mark -- ##
## ---------------------------------------------------------------------------- ##

## Do before getting rid of genes without quantified expression
GO_df <- Enh_int |> filter(D4_class_ext %in% c("actv", "actvM", "prim") & DiffBind.direction != "none" & !is.na(gene_symbol))
GO_df <- GO_df |> select(chr, start, end, D4_class_ext, DiffBind.direction, gene_symbol)

Enh_class <- c("actv", "actvM", "prim")
DB_direction <- c("up", "down")

#Output_folder <- "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/GeneOntologies/ROIs"

for(i in 1:length(Enh_class)) {
  for(j in 1:length(DB_direction)) {
    Sub_GO <- GO_df |> filter(D4_class_ext == Enh_class[i] & DiffBind.direction == DB_direction[j])
    Sub_GO <- Sub_GO |> select(gene_symbol) |> distinct() |> arrange(gene_symbol)

    Name_file <- paste0(group1,"_vs_",group2,"_",mark,"_",Enh_class[i],"_",DB_direction[j],"_AssociatedGenes.txt")
    #write.table(x = Sub_GO, file = paste0(Output_folder,"/",Name_file), append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
  }
}
rm(GO_df, Enh_class, DB_direction, Output_folder, i, j, Sub_GO, Name_file)
## -----------------
## ~~ End of Break


## -- Stats
## Total genes
#cat("Number of total expressed genes: ",length(unique(RNAseq_res$gene_symbol)),"\n", sep = "", file = Name_log, append = TRUE)
#cat("Down:\t",sum(RNAseq_res$RNA.direction == "down"),"\n", sep = "", file = Name_log, append = TRUE)
#cat("Up:\t",sum(RNAseq_res$RNA.direction == "up"),"\n\n", sep = "", file = Name_log, append = TRUE)

## Original mapping of enhancers to genes
#cat("Number of genes linked to intersected enhancers: ",length(unique(Enh_int$gene_symbol)),"\n\n", sep = "", file = Name_log, append = TRUE)

## Mapped genes that are expressed
## Drop genes without counts in both samples (logFC==NA); will also drop enhancers without an associated gene
Enh_int <- Enh_int |> filter(!is.na(RNA.log2FoldChange))
#tmp <- Enh_int[,c("gene_symbol", "RNA.direction")] |> distinct()
#cat("Number of expressed, linked genes: ",nrow(tmp),"\n", sep = "", file = Name_log, append = TRUE)
#cat("Down:\t",sum(tmp$RNA.direction == "down"),"\n", sep = "", file = Name_log, append = TRUE)
#cat("Up:\t",sum(tmp$RNA.direction == "up"),"\n\n", sep = "", file = Name_log, append = TRUE)
#rm(tmp)

## From the enhancer point of view
# tmp_class <- Enh_int[,c("chr","start","end","coord","peak_id",paste0(day1,"_class_ext"),"DiffBind.direction")] |> distinct() |>
#   group_by(.data[[paste0(day1,"_class_ext")]]) |> summarise("down" = sum(DiffBind.direction == "down"),
#                                                             "none" = sum(DiffBind.direction == "none"),
#                                                             "up" = sum(DiffBind.direction == "up"),
#                                                             "Total" = n())
# tmp_Tot <- Enh_int[,c("chr","start","end","coord","peak_id",paste0(day1,"_class_ext"),"DiffBind.direction")] |> distinct() |>
#   summarise("down" = sum(DiffBind.direction == "down"),
#             "none" = sum(DiffBind.direction == "none"),
#             "up" = sum(DiffBind.direction == "up"),
#             "Total" = n())
#cat("Number of intersected enhancers linked to expressed genes: (",length(unique(Enh_int$peak_id)),")\n", sep = "", file = Name_log, append = TRUE)
#cat("\tdown\tnone\tup\tTotal\n", sep = "", file = Name_log, append = TRUE)
#write.table(x = tmp_class, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
#cat("----------------------------------------\n", sep = "", file = Name_log, append = TRUE)
#cat("Total\t", sep = "", file = Name_log, append = TRUE)
#write.table(x = tmp_Tot, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = FALSE)
#cat("\n", sep = "", file = Name_log, append = TRUE)

rm(RNAseq_res, RNAseq_name, tmp_class, tmp_Tot)

## ------------------------------------------
## ------------------------------------------

## -- Part 9: Plots Enhancers <=> DiffBind <=> RNAseq ----

## Make a df copy specifically for plotting. Re-sort direction factors in DiffBind and RNA
Enh_int_plot <- Enh_int
Enh_int_plot$DiffBind.direction <- Enh_int_plot$DiffBind.direction |> forcats::fct_relevel(c("none", "up", "down"))
Enh_int_plot$RNA.direction <- Enh_int_plot$RNA.direction |> forcats::fct_relevel(c("none", "up", "down"))
Enh_int_plot <- Enh_int_plot |> arrange(DiffBind.direction, RNA.direction)

## ________________________________________________________
## __ Enhancer POV >> Assess RNA differential expression __
## ________________________________________________________

## -- DiffBind.logFC vs RNA.logFC; color by DiffBind DE; facet by enhancer class
## All enhancer classes
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_RNAseqFC_byEnhClassExt_full.png")
p <- ggplot(Enh_int_plot, aes(x = DiffBind.logFC, y = RNA.log2FoldChange, color = DiffBind.direction, fill = DiffBind.direction, shape = RNA.direction)) +
  geom_point(na.rm = TRUE, size = 1.2, stroke = 0.35) +
  scale_fill_manual(values = c("gray70", as.vector(gen_colors[c(gen2, gen1)]))) +
  scale_color_manual(values = c("gray70", "gray20", "gray20")) +
  scale_shape_manual(values = c(21,24,25)) +
  labs(x = paste0("log2FC(",mark,")"), y = "log2FC(RNA)") +
  geom_hline(yintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 2)
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
rm(p)

## Only actv, actvM, and prim enhancers
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_RNAseqFC_byEnhClassExt_aap.png")
p <- ggplot(Enh_int_plot[Enh_int_plot[,paste0(day1,"_class_ext")] %in% c("actv", "actvM", "prim"),],
            aes(x = DiffBind.logFC, y = RNA.log2FoldChange, color = DiffBind.direction, fill = DiffBind.direction, shape = RNA.direction)) +
  geom_point(na.rm = TRUE, size = 1.2, stroke = 0.35) +
  scale_fill_manual(values = c("gray70", as.vector(gen_colors[c(gen2, gen1)]))) +
  scale_color_manual(values = c("gray70", "gray20", "gray20")) +
  scale_shape_manual(values = c(21,24,25)) +
  labs(x = paste0("log2FC(",mark,")"), y = "log2FC(RNA)") +
  geom_hline(yintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 30, height = 10, units = "cm", create.dir = TRUE)
rm(p)

## ________

## -- DiffBind.logFC vs RNA.logFC; doing the inverse --> color by RNA DE; facet by enhancer class
## All enhancer classes
Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_RNAseqFC_byEnhClassExt_full_inverse.png")
p <- ggplot(Enh_int_plot, aes(x = DiffBind.logFC, y = RNA.log2FoldChange, color = RNA.direction, fill = RNA.direction, shape = DiffBind.direction)) +
  geom_point(na.rm = TRUE, size = 1.2, stroke = 0.35) +
  scale_fill_manual(values = c("gray70", as.vector(gen_colors[c(gen2, gen1)]))) +
  scale_color_manual(values = c("gray70", "gray20", "gray20")) +
  scale_shape_manual(values = c(21,24,25)) +
  labs(x = paste0("log2FC(",mark,")"), y = "log2FC(RNA)") +
  geom_hline(yintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 2)
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
rm(p)

## Only actv, actvM, and prim enhancers
Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindFC_vs_RNAseqFC_byEnhClassExt_aap_inverse.png")
p <- ggplot(Enh_int_plot[Enh_int_plot[,paste0(day1,"_class_ext")] %in% c("actv", "actvM", "prim"),],
       aes(x = DiffBind.logFC, y = RNA.log2FoldChange, color = RNA.direction, fill = RNA.direction, shape = DiffBind.direction)) +
  geom_point(na.rm = TRUE, size = 1.2, stroke = 0.35) +
  scale_fill_manual(values = c("gray70", as.vector(gen_colors[c(gen2, gen1)]))) +
  scale_color_manual(values = c("gray70", "gray20", "gray20")) +
  scale_shape_manual(values = c(21,24,25)) +
  labs(x = paste0("log2FC(",mark,")"), y = "log2FC(RNA)") +
  geom_hline(yintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 30, height = 10, units = "cm", create.dir = TRUE)
rm(p)

## ________

## -- DiffBind.DE vs RNA.logFC; facet by enhancer class, omit inac
## Consider only unique observations
tmp_int <- Enh_int |> filter(.data[[paste0(day1,"_class_ext")]] != "inac") |> select(.data[[paste0(day1, "_class_ext")]], DiffBind.direction, gene_symbol, RNA.log2FoldChange, RNA.direction) |> distinct()
## Plot
#Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindDE_vs_RNAseqFC_byEnhClassExt_aap.png")
p <- ggplot(tmp_int, aes(x = DiffBind.direction, y = RNA.log2FoldChange, fill = DiffBind.direction)) +
  geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
  labs(x = paste0(mark," DE"), y = "log2FC(RNA)") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
p
#ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 10, units = "cm", create.dir = TRUE)

## ------------------------------ ##
## -- For paper, export in PDF -- ##
## ------------------------------ ##

if(day1 %in% c("D2", "D4")){
  ## Only export
  Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindDE_vs_RNAseqFC_byEnhClassExt_aap.pdf")
  # ggsave(filename = Name_plot, plot = p, device = "pdf", width = 20, height = 10, units = "cm", create.dir = FALSE,
  #        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_6/Supp-X_Boxplot_K27acDE_vs_RNAexp_by_EnhancerCat/")
  rm(p)
} else if (day1 == "D0") {
  ## Omit actvM; re-do plot
  tmp_int <- Enh_int |> filter(.data[[paste0(day1,"_class_ext")]] %in% c("actv", "prim")) |> select(.data[[paste0(day1, "_class_ext")]], DiffBind.direction, gene_symbol, RNA.log2FoldChange, RNA.direction) |> distinct()
  p <- ggplot(tmp_int, aes(x = DiffBind.direction, y = RNA.log2FoldChange, fill = DiffBind.direction)) +
    geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
    geom_boxplot(outliers = FALSE) +
    scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
    labs(x = paste0(mark," DE"), y = "log2FC(RNA)") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
  p
  ## Export as pdf
  Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DiffBindDE_vs_RNAseqFC_byEnhClassExt_ap.pdf")
  # ggsave(filename = Name_plot, plot = p, device = "pdf", width = 14, height = 10, units = "cm", create.dir = FALSE,
  #        path = "/Users/user/Dropbox/Greenberg Lab Sharing/Collaborations/MethylEnh/Manuscript/Figures/Figure_6/Supp-X_Boxplot_K27acDE_vs_RNAexp_by_EnhancerCat/")
  rm(p)
}


# ## -- Export stats
# t_Output <- data.frame("Enh_class" = character(),
#                        "GroupA" = character(),
#                        "GroupB" = character(),
#                        "A_n" = numeric(),
#                        "B_n" = numeric(),
#                        "A_mean_logFC" = numeric(),
#                        "B_mean_logFC" = numeric(),
#                        "A_median_logFC" = numeric(),
#                        "B_median_logFC" = numeric(),
#                        "pval" = numeric())
# for(i in 1:length(Enh_classes)) {
#   ## Pull RNA.logFC
#   Down <- tmp_int |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] & DiffBind.direction == "down") |> pull(RNA.log2FoldChange)
#   None <- tmp_int |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] & DiffBind.direction == "none") |> pull(RNA.log2FoldChange)
#   Up <- tmp_int |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] & DiffBind.direction == "up") |> pull(RNA.log2FoldChange)
#
#   ## Test comparisons -- minimum sample size 6
#   if(min(length(Down), length(None), length(Up)) >= 6) {
#     t_DN <- wilcox.test(x = Down, y = None)
#     t_UN <- wilcox.test(x = Up, y = None)
#     t_DU <- wilcox.test(x = Down, y = Up)
#
#     ## Append information to output
#     t_Output <- rbind(t_Output, data.frame("Enh_class" = Enh_classes[i],
#                                            "GroupA" = c("down", "up", "down"),
#                                            "GroupB" = c("none", "none", "up"),
#                                            "A_n" = c(length(Down), length(Up), length(Down)),
#                                            "B_n" = c(length(None), length(None), length(Up)),
#                                            "A_mean_logFC" = round(c(mean(Down), mean(Up), mean(Down)), digits = 4),
#                                            "B_mean_logFC" = round(c(mean(None), mean(None), mean(Up)), digits = 4),
#                                            "A_median_logFC" = round(c(median(Down), median(Up), median(Down)), digits = 4),
#                                            "B_median_logFC" = round(c(median(None), median(None), median(Up)), digits = 4),
#                                            "pval" = c(t_DN$p.value, t_UN$p.value, t_DU$p.value)))
#   }
# }
# ## Clean output
# t_Output$Enh_class <- factor(t_Output$Enh_class, levels = Enh_classes)
# t_Output$GroupA <- factor(t_Output$GroupA, levels = c("down", "none", "up"))
# t_Output$GroupB <- factor(t_Output$GroupB, levels = c("down", "none", "up"))
# t_Output <- t_Output |> arrange(Enh_class, GroupA, GroupB)
# t_Output$paj <- p.adjust(t_Output$pval, method = "bonferroni")
# t_Output$sig <- t_Output$paj < 0.05
# ## Round to four significant decimals
# t_Output$pval <- signif(t_Output$pval, digits = 5)
# t_Output$paj <- signif(t_Output$paj, digits = 5)
#
# #cat("\n== Enhancer POV: Analyzing differential RNAseq ==\n\n", sep = "", file = Name_log, append = TRUE)
# #cat("Calculating two-sided wilcox test of RNA.logFC(",group2,"/⁠",group1,") per ",mark," direction per enhancer class\n\n", sep = "", file = Name_log, append = TRUE)
# #write.table(x = t_Output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# ## Clean workspace
# rm(tmp_int, Down, None, Up, t_DN, t_UN, t_DU, t_Output, i)
#
# ## ________________________________________________________
# ## __ RNA POV >> Assess Enhancer differential Enrichment __
# ## ________________________________________________________
#
# ## Approach it from the opposite direction.
# ## Dif Exp genes, to which enhancers are they linked?
# #cat("\n== RNA POV: Analyzing differential chrom mark ==\n\n", sep = "", file = Name_log, append = TRUE)
# #cat("Enhancer classes associated to DE genes:\n", sep = "", file = Name_log, append = TRUE)
# tmp <- table(Enh_int[,c("RNA.direction", paste0(day1, "_class_ext"))], exclude = "othr")
# #cat("\t", sep = "", file = Name_log, append = TRUE)
# #write.table(x = tmp, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = TRUE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
# rm(tmp)
#
# ## -- RNA.logFC vs DiffBind.logFC; color by RNA DE; facet by enhancer class
# ## Transposed version of plot above
# Enh_int_plot <- Enh_int_plot |> arrange(RNA.direction, DiffBind.direction)
#
# ## All enhancer classes
# #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_RNAseqFC_vs_DiffBindFC_byEnhClassExt_full.png")
# p <- ggplot(Enh_int_plot, aes(x = RNA.log2FoldChange, y = DiffBind.logFC, color = RNA.direction, fill = RNA.direction, shape = DiffBind.direction)) +
#   geom_point(na.rm = TRUE, size = 1.2, stroke = 0.35) +
#   scale_fill_manual(values = c("gray70", as.vector(gen_colors[c(gen2, gen1)]))) +
#   scale_color_manual(values = c("gray70", "gray20", "gray20")) +
#   scale_shape_manual(values = c(21,24,25)) +
#   labs(x = "log2FC(RNA)", y = paste0("log2FC(",mark,")")) +
#   geom_hline(yintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
#   geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#   facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 2)
# p
# #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
# rm(p)
#
# ## Only actv, actvM, and prim enhancers
# #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_RNAseqFC_vs_DiffBindFC_byEnhClassExt_aap.png")
# p <- ggplot(Enh_int_plot[Enh_int_plot[,paste0(day1,"_class_ext")] %in% c("actv", "actvM", "prim"),],
#             aes(x = RNA.log2FoldChange, y = DiffBind.logFC, color = RNA.direction, fill = RNA.direction, shape = DiffBind.direction)) +
#   geom_point(na.rm = TRUE, size = 1.2, stroke = 0.35) +
#   scale_fill_manual(values = c("gray70", as.vector(gen_colors[c(gen2, gen1)]))) +
#   scale_color_manual(values = c("gray70", "gray20", "gray20")) +
#   scale_shape_manual(values = c(21,24,25)) +
#   labs(x = "log2FC(RNA)", y = paste0("log2FC(",mark,")")) +
#   geom_hline(yintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
#   geom_vline(xintercept = c(-1,1), color = "gray0", linewidth = 0.2, linetype = "dashed") +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#   facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
# p
# #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 30, height = 10, units = "cm", create.dir = TRUE)
# rm(p)
#
# rm(Enh_int_plot)
#
# ## ------------------------------------------ ##
# ## -- logFC(K27ac) for each gene direction -- ##
#
# ## -- RNA.DE vs DiffBind.logFC; color by RNA DE; facet by enhancer class
# ## All enhancer classes together
# #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_RNAseqDE_vs_DiffBindFC_full.png")
# p <- ggplot(Enh_int, aes(x = RNA.direction, y = DiffBind.logFC, fill = RNA.direction)) +
#   geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
#   geom_boxplot(outliers = FALSE) +
#   scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
#   #scale_x_discrete(limits = c("up", "none", "down")) +
#   labs(y = paste0("logFC(",mark,")")) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
# p
# #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 8, height = 10, units = "cm", create.dir = TRUE)
# rm(p)
#
# ## -- Statistics
# Up <- Enh_int |> filter(RNA.direction == "up") |> pull(DiffBind.logFC)
# None <- Enh_int |> filter(RNA.direction == "none") |> pull(DiffBind.logFC)
# Down <- Enh_int |> filter(RNA.direction == "down") |> pull(DiffBind.logFC)
#
# ## Define empty df before filling them
# Comp_UN <- data.frame("GroupA" = character(),
#                       "GroupB" = character(),
#                       "A_n" = numeric(),
#                       "B_n" = numeric(),
#                       "A_mean_logFC" = numeric(),
#                       "B_mean_logFC" = numeric(),
#                       "A_median_logFC" = numeric(),
#                       "B_median_logFC" = numeric(),
#                       "pval_greater" = numeric(),
#                       "pval_less" = numeric())
# Comp_DN <- Comp_UN
# Comp_UD <- Comp_UN
#
# ## Calculate p-vals
# if(min(length(Up), length(None)) >= 6) {
#   t.UN.greater <- wilcox.test(x = Up, y = None, alternative = "greater")
#   t.UN.less <- wilcox.test(x = Up, y = None, alternative = "less")
#   Comp_UN <- data.frame("GroupA" = "up",
#                         "GroupB" = "none",
#                         "A_n" = length(Up),
#                         "B_n" = length(None),
#                         "A_mean_logFC" = round(mean(Up), digits = 4),
#                         "B_mean_logFC" = round(mean(None), digits = 4),
#                         "A_median_logFC" = round(median(Up), digits = 4),
#                         "B_median_logFC" = round(median(None), digits = 4),
#                         "pval_greater" = t.UN.greater$p.value,
#                         "pval_less" = t.UN.less$p.value)
# }
# if(min(length(Down), length(None)) >= 6) {
#   t.DN.greater <- wilcox.test(x = Down, y = None, alternative = "greater")
#   t.DN.less <- wilcox.test(x = Down, y = None, alternative = "less")
#   Comp_DN <- data.frame("GroupA" = "down",
#                         "GroupB" = "none",
#                         "A_n" = length(Down),
#                         "B_n" = length(None),
#                         "A_mean_logFC" = round(mean(Down), digits = 4),
#                         "B_mean_logFC" = round(mean(None), digits = 4),
#                         "A_median_logFC" = round(median(Down), digits = 4),
#                         "B_median_logFC" = round(median(None), digits = 4),
#                         "pval_greater" = t.DN.greater$p.value,
#                         "pval_less" = t.DN.less$p.value)
# }
# if(min(length(Up), length(Down)) >= 6) {
#   t.UD.greater <- wilcox.test(x = Up, y = Down, alternative = "greater")
#   t.UD.less <- wilcox.test(x = Up, y = Down, alternative = "less")
#   Comp_UD <- data.frame("GroupA" = "up",
#                         "GroupB" = "down",
#                         "A_n" = length(Up),
#                         "B_n" = length(Down),
#                         "A_mean_logFC" = round(mean(Up), digits = 4),
#                         "B_mean_logFC" = round(mean(Down), digits = 4),
#                         "A_median_logFC" = round(median(Up), digits = 4),
#                         "B_median_logFC" = round(median(Down), digits = 4),
#                         "pval_greater" = t.UD.greater$p.value,
#                         "pval_less" = t.UD.less$p.value)
# }
# Comp_output <- rbind(Comp_UN, Comp_DN, Comp_UD)
# Comp_output$paj_greater <- p.adjust(Comp_output$pval_greater, method = "bonferroni")
# Comp_output$paj_less <- p.adjust(Comp_output$pval_less, method = "bonferroni")
# Comp_output$sig_greater <- Comp_output$paj_greater < 0.05
# Comp_output$sig_less <- Comp_output$paj_less < 0.05
# ## Round to four significant decimals
# Comp_output$pval_greater <- signif(Comp_output$pval_greater, digits = 5)
# Comp_output$pval_less <- signif(Comp_output$pval_less, digits = 5)
# Comp_output$paj_greater <- signif(Comp_output$paj_greater, digits = 5)
# Comp_output$paj_less <- signif(Comp_output$paj_less, digits = 5)
#
# #cat("Calculating wilcox test of logFC(",mark,") per RNA direction\n\n", sep = "", file = Name_log, append = TRUE)
# #write.table(x = Comp_output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# rm(Up, None, Down, t.UN.greater, t.UN.less, t.DN.greater, t.DN.less, t.UD.greater, t.UD.less, Comp_UN, Comp_DN, Comp_UD, Comp_output)
#
# ## ------------------------------------------------------------------- ##
# ## -- logFC(K27ac) for each gene direction; facet by enhancer class -- ##
#
# ## -- RNA.DE vs DiffBind.logFC; color by RNA DE; facet by enhancer class
# #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_RNAseqDE_vs_DiffBindFC_byEnhClassExt_aap.png")
# p <- ggplot(Enh_int[Enh_int[,paste0(day1,"_class_ext")] != "inac",], aes(x = RNA.direction, y = DiffBind.logFC, fill = RNA.direction)) +
#   geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
#   geom_boxplot(outliers = FALSE) +
#   scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
#   #scale_x_discrete(limits = c("up", "none", "down")) +
#   labs(y = paste0("logFC(",mark,")")) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#   facet_wrap( {{paste(day1,"_class_ext",sep = "")}}, nrow = 1)
# p
# #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 16, height = 10, units = "cm", create.dir = TRUE)
# rm(p)
#
# ## -- Statistics by enhancer class
# Comp_output <- data.frame("Enh_class" = character(),
#                           "GroupA" = character(),
#                           "GroupB" = character(),
#                           "A_n" = numeric(),
#                           "B_n" = numeric(),
#                           "A_mean_logFC" = numeric(),
#                           "B_mean_logFC" = numeric(),
#                           "A_median_logFC" = numeric(),
#                           "B_median_logFC" = numeric(),
#                           "pval_greater" = numeric(),
#                           "pval_less" = numeric())
#
# for(i in 1:length(Enh_classes)) {
#   ## Pull DiffBind logFC
#   Up <- Enh_int |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] & RNA.direction == "up") |> pull(DiffBind.logFC)
#   None <- Enh_int |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] & RNA.direction == "none") |> pull(DiffBind.logFC)
#   Down <- Enh_int |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] & RNA.direction == "down") |> pull(DiffBind.logFC)
#
#   ## Define empty df before filling them
#   Comp_UN <- data.frame("Enh_class" = character(),
#                         "GroupA" = character(),
#                         "GroupB" = character(),
#                         "A_n" = numeric(),
#                         "B_n" = numeric(),
#                         "A_mean_logFC" = numeric(),
#                         "B_mean_logFC" = numeric(),
#                         "A_median_logFC" = numeric(),
#                         "B_median_logFC" = numeric(),
#                         "pval_greater" = numeric(),
#                         "pval_less" = numeric())
#   Comp_DN <- Comp_UN
#   Comp_UD <- Comp_UN
#
#   ## Calculate p-vals
#   if(min(length(Up), length(None)) >= 6) {
#     t.UN.greater <- wilcox.test(x = Up, y = None, alternative = "greater")
#     t.UN.less <- wilcox.test(x = Up, y = None, alternative = "less")
#     Comp_UN <- data.frame("Enh_class" = Enh_classes[i],
#                           "GroupA" = "up",
#                           "GroupB" = "none",
#                           "A_n" = length(Up),
#                           "B_n" = length(None),
#                           "A_mean_logFC" = round(mean(Up), digits = 4),
#                           "B_mean_logFC" = round(mean(None), digits = 4),
#                           "A_median_logFC" = round(median(Up), digits = 4),
#                           "B_median_logFC" = round(median(None), digits = 4),
#                           "pval_greater" = t.UN.greater$p.value,
#                           "pval_less" = t.UN.less$p.value)
#   }
#   if(min(length(Down), length(None)) >= 6) {
#     t.DN.greater <- wilcox.test(x = Down, y = None, alternative = "greater")
#     t.DN.less <- wilcox.test(x = Down, y = None, alternative = "less")
#     Comp_DN <- data.frame("Enh_class" = Enh_classes[i],
#                           "GroupA" = "down",
#                           "GroupB" = "none",
#                           "A_n" = length(Down),
#                           "B_n" = length(None),
#                           "A_mean_logFC" = round(mean(Down), digits = 4),
#                           "B_mean_logFC" = round(mean(None), digits = 4),
#                           "A_median_logFC" = round(median(Down), digits = 4),
#                           "B_median_logFC" = round(median(None), digits = 4),
#                           "pval_greater" = t.DN.greater$p.value,
#                           "pval_less" = t.DN.less$p.value)
#   }
#   if(min(length(Up), length(Down)) >= 6) {
#     t.UD.greater <- wilcox.test(x = Up, y = Down, alternative = "greater")
#     t.UD.less <- wilcox.test(x = Up, y = Down, alternative = "less")
#     Comp_UD <- data.frame("Enh_class" = Enh_classes[i],
#                           "GroupA" = "up",
#                           "GroupB" = "down",
#                           "A_n" = length(Up),
#                           "B_n" = length(Down),
#                           "A_mean_logFC" = round(mean(Up), digits = 4),
#                           "B_mean_logFC" = round(mean(Down), digits = 4),
#                           "A_median_logFC" = round(median(Up), digits = 4),
#                           "B_median_logFC" = round(median(Down), digits = 4),
#                           "pval_greater" = t.UD.greater$p.value,
#                           "pval_less" = t.UD.less$p.value)
#   }
#   Comp_output <- rbind(Comp_output, Comp_UN, Comp_DN, Comp_UD)
# }
# Comp_output$paj_greater <- p.adjust(Comp_output$pval_greater, method = "bonferroni")
# Comp_output$paj_less <- p.adjust(Comp_output$pval_less, method = "bonferroni")
# Comp_output$sig_greater <- Comp_output$paj_greater < 0.05
# Comp_output$sig_less <- Comp_output$paj_less < 0.05
# ## Round to four significant decimals
# Comp_output$pval_greater <- signif(Comp_output$pval_greater, digits = 5)
# Comp_output$pval_less <- signif(Comp_output$pval_less, digits = 5)
# Comp_output$paj_greater <- signif(Comp_output$paj_greater, digits = 5)
# Comp_output$paj_less <- signif(Comp_output$paj_less, digits = 5)
#
# #cat("Calculating wilcox test of logFC(",mark,") per RNA direction per enhancer class\n\n", sep = "", file = Name_log, append = TRUE)
# #write.table(x = Comp_output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# ## Clear workspace
# rm(Up, None, Down, t.UN.greater, t.UN.less, t.DN.greater, t.DN.less, t.UD.greater, t.UD.less, Comp_UN, Comp_DN, Comp_UD, Comp_output, i)
#
# ## ------------------------------------------
# ## ------------------------------------------
#
# ## -- Part 10: Analyze relationship between DiffBind and RNAseq data ----
# #cat("\n== Enhancer POV ∩ RNA POV ==\n\n", sep = "", file = Name_log, append = TRUE)
#
# tmp <- table(Enh_int$DiffBind.direction,
#              Enh_int$RNA.direction)
# rownames(tmp) <- paste(mark, rownames(tmp), sep = ".")
# colnames(tmp) <- paste("RNA", colnames(tmp), sep = ".")
# #cat("Relationship DiffBind DE <-> RNA DE\n\n",sep = "", file = Name_log, append = TRUE)
# #cat("\t",sep = "", file = Name_log, append = TRUE)
# #write.table(x = tmp, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = TRUE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
# rm(tmp)
#
# #cat("Relationship DiffBind DE <-> RNA DE per enhancer class\n\n",sep = "", file = Name_log, append = TRUE)
#
# ## -- Active
# tmp_actv <- Enh_int |>
#   filter(.data[[paste0(day1,"_class_ext")]] == "actv") |>
#   select(DiffBind.direction, RNA.direction) |>
#   table()
# rownames(tmp_actv) <- paste(mark, rownames(tmp_actv), sep = ".")
# colnames(tmp_actv) <- paste("RNA", colnames(tmp_actv), sep = ".")
# #cat("[actv]\t",sep = "", file = Name_log, append = TRUE)
# #write.table(x = tmp_actv, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = TRUE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# ## -- Active Methylated
# tmp_actvM <- Enh_int |>
#   filter(.data[[paste0(day1,"_class_ext")]] == "actvM") |>
#   select(DiffBind.direction, RNA.direction) |>
#   table()
# rownames(tmp_actvM) <- paste(mark, rownames(tmp_actvM), sep = ".")
# colnames(tmp_actvM) <- paste("RNA", colnames(tmp_actvM), sep = ".")
# #cat("[actvM]\t",sep = "", file = Name_log, append = TRUE)
# #write.table(x = tmp_actvM, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = TRUE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# ## -- Primed
# tmp_prim <- Enh_int |>
#   filter(.data[[paste0(day1,"_class_ext")]] == "prim") |>
#   select(DiffBind.direction, RNA.direction) |>
#   table()
# rownames(tmp_prim) <- paste(mark, rownames(tmp_prim), sep = ".")
# colnames(tmp_prim) <- paste("RNA", colnames(tmp_prim), sep = ".")
# #cat("[prim]\t",sep = "", file = Name_log, append = TRUE)
# #write.table(x = tmp_prim, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = TRUE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# ## -- Inaccessible
# tmp_inac <- Enh_int |>
#   filter(.data[[paste0(day1,"_class_ext")]] == "inac") |>
#   select(DiffBind.direction, RNA.direction) |>
#   table()
# rownames(tmp_inac) <- paste(mark, rownames(tmp_inac), sep = ".")
# colnames(tmp_inac) <- paste("RNA", colnames(tmp_inac), sep = ".")
# #cat("[inac]\t",sep = "", file = Name_log, append = TRUE)
# #write.table(x = tmp_inac, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = TRUE, col.names = TRUE)
# #cat("\n", sep = "", file = Name_log, append = TRUE)
#
# rm(tmp_actv, tmp_actvM, tmp_prim, tmp_inac)
#
# ## ------------------------------------------
# ## ------------------------------------------
#
# ## -- Part 11: Analyze relationship between DNAme and RNAseq data ----
# #cat("\n== Analyzing relationship between DNAme and RNAseq ==\n\n", sep = "", file = Name_log, append = TRUE)
#
# ## Make a copy of Enh_int for pivoting
# Enh_int_pivot <- Enh_int |> tidyr::pivot_longer(cols = contains("_meth.wmean"),
#                                                 names_to = c("Group", "Mark"),
#                                                 values_to = "meth.wmean",
#                                                 names_pattern = "^(.+)_([^_]+)_meth\\.wmean$")
# ## Drop missing meth.wmean
# Enh_int_pivot <- Enh_int_pivot |> filter(!is.na(meth.wmean))
# ## Factorize columns
# Enh_int_pivot$Group <- factor(Enh_int_pivot$Group, levels = c(group1, group2))
# Enh_int_pivot$Mark <- factor(Enh_int_pivot$Mark, levels = c("K4me1", "MATAC", "K27ac", "WGBS"))
#
# ## Remove duplicated counts (same enhancer linked to multiple genes with the same direction)
# Enh_int_pivot_dedup <- Enh_int_pivot |> select(!gene_symbol:RNA.padj) |> distinct()
#
# ## -- Plot DNAme for each genotype, for each mark, faceted by enhancer class (by RNAseq.DE)
# ## All enhancer classes
# #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAme_byMark_byEnhClassExt_byRNAseqDE_full.png")
# p <- ggplot(Enh_int_pivot_dedup, aes(x = RNA.direction, y = meth.wmean, fill = Group)) +
#   geom_boxplot(alpha = 0.9, outliers = FALSE) +
#   scale_fill_manual(values = as.vector(gen_colors[c(gen1,gen2)])) +
#   labs(x = "RNA DE", y = "DNAme %") +
#   coord_cartesian(ylim = c(0,100)) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#   facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
# p
# #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
# rm(p)
#
# ## Only actv, actvM, and prim enhancers; no K27ac 5mC
# tmp_aap <- Enh_int_pivot_dedup |> filter(Mark != "K27ac" & .data[[paste0(day1, "_class_ext")]] %in% c("actv", "actvM", "prim"))
# #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAme_byMark_byEnhClassExt_byRNAseqDE_aap.png")
# p <- ggplot(tmp_aap, aes(x = RNA.direction, y = meth.wmean, fill = Group)) +
#   geom_boxplot(alpha = 0.9, outliers = FALSE) +
#   scale_fill_manual(values = as.vector(gen_colors[c(gen1,gen2)])) +
#   labs(x = "RNA DE", y = "DNAme %") +
#   coord_cartesian(ylim = c(0,100)) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#   facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
# p
# #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 16, units = "cm", create.dir = TRUE)
# rm(p, tmp_aap)
# ## ----
#
# ## -- Stats for DNAme WT of the above plot ----
# ## Only useful if comparing WT vs TKO
#
# if(gen1=="WT" & gen2=="TKO") {
#   t_Output <- data.frame("Enh_class" = character(),
#                          "Mark" = character(),
#                          "GroupA" = character(),
#                          "GroupB" = character(),
#                          "A_n" = numeric(),
#                          "B_n" = numeric(),
#                          "A_mean_5mC" = numeric(),
#                          "B_mean_5mC" = numeric(),
#                          "A_median_5mC" = numeric(),
#                          "B_median_5mC" = numeric(),
#                          "pval" = numeric())
#   for(i in 1:length(Enh_classes)) {
#     for(j in 1:length(Mark_classes)) {
#       ## Pull DNAme values
#       Down <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                               Mark == Mark_classes[j] & RNA.direction == "down" & Group == group1) |> pull(meth.wmean)
#       None <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                               Mark == Mark_classes[j] & RNA.direction == "none" & Group == group1) |> pull(meth.wmean)
#       Up <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                             Mark == Mark_classes[j] & RNA.direction == "up" & Group == group1) |> pull(meth.wmean)
#
#       ## Test comparisons -- minimum sample size 6
#       if(min(length(Down), length(None), length(Up)) >= 6) {
#         t_DN <- wilcox.test(x = Down, y = None)
#         t_UN <- wilcox.test(x = Up, y = None)
#         t_DU <- wilcox.test(x = Down, y = Up)
#
#         ## Append information to output
#         t_Output <- rbind(t_Output, data.frame("Enh_class" = Enh_classes[i],
#                                                "Mark" = Mark_classes[j],
#                                                "GroupA" = c("down", "up", "down"),
#                                                "GroupB" = c("none", "none", "up"),
#                                                "A_n" = c(length(Down), length(Up), length(Down)),
#                                                "B_n" = c(length(None), length(None), length(Up)),
#                                                "A_mean_5mC" = round(c(mean(Down), mean(Up), mean(Down)), digits = 2),
#                                                "B_mean_5mC" = round(c(mean(None), mean(None), mean(Up)), digits = 2),
#                                                "A_median_5mC" = round(c(median(Down), median(Up), median(Down)), digits = 2),
#                                                "B_median_5mC" = round(c(median(None), median(None), median(Up)), digits = 2),
#                                                "pval" = c(t_DN$p.value, t_UN$p.value, t_DU$p.value)))
#       }
#     }
#   }
#   ## Clean output
#   t_Output$Enh_class <- factor(t_Output$Enh_class, levels = Enh_classes)
#   t_Output$Mark <- factor(t_Output$Mark, levels = Mark_classes)
#   t_Output$GroupA <- factor(t_Output$GroupA, levels = c("down", "none", "up"))
#   t_Output$GroupB <- factor(t_Output$GroupB, levels = c("down", "none", "up"))
#   t_Output <- t_Output |> arrange(Enh_class, GroupA, GroupB, Mark)
#   t_Output <- t_Output |> relocate(Mark, .after = GroupB)
#   t_Output$paj <- p.adjust(t_Output$pval, method = "bonferroni")
#   t_Output$sig <- t_Output$paj < 0.05
#   ## Round to four significant decimals
#   t_Output$pval <- signif(t_Output$pval, digits = 5)
#   t_Output$paj <- signif(t_Output$paj, digits = 5)
#
#   #cat("Calculating two-sided wilcox test (5mC ",group1,") per RNA direction per enhancer class\n\n", sep = "", file = Name_log, append = TRUE)
#   #write.table(x = t_Output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
#   #cat("\n", sep = "", file = Name_log, append = TRUE)
#   ## Clean workspace
#   rm(Down, None, Up, t_DN, t_UN, t_DU, t_Output, i, j)
# }
# ## ----
#
# ## -- Stats for DNAme delta of the above plot----
# ## Only useful if working with two DNAme datasets (WT and/or DKO)
#
# if(gen1!="TKO" & gen2!="TKO") {
#   ## Create column for easier operation
#   Enh_int_pivot_dedup <- Enh_int_pivot_dedup |> tidyr::pivot_wider(names_from = Group,
#                                                                    names_glue = "{Group}_meth.wmean",
#                                                                    values_from = meth.wmean)
#   Enh_int_pivot_dedup <- Enh_int_pivot_dedup |> mutate("DNAme_delta" = .data[[paste0(group2,"_meth.wmean")]] - .data[[paste0(group1,"_meth.wmean")]])
#
#   ## -- Plot d5mC, for each mark, faceted by enhancer class
#   ## All enhancer classes
#   #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_byMark_byEnhClassExt_byRNAseqDE_full.png")
#   p <- ggplot(Enh_int_pivot_dedup, aes(x = RNA.direction, y = DNAme_delta, fill = RNA.direction)) +
#     geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
#     geom_boxplot(outliers = FALSE, na.rm = TRUE) +
#     scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
#     labs(x = "RNA DE", y = paste0("DNAme delta (",group2," - ",group1,")")) +
#     theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#     facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
#   p
#   #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 18, units = "cm", create.dir = TRUE)
#   rm(p)
#
#   ## Only actv, actvM, and prim enhancers; no K27ac 5mC
#   tmp_aap <- Enh_int_pivot_dedup |> filter(Mark != "K27ac" & .data[[paste0(day1, "_class_ext")]] %in% c("actv", "actvM", "prim"))
#   #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_byMark_byEnhClassExt_byRNAseqDE_aap.png")
#   p <- ggplot(tmp_aap, aes(x = RNA.direction, y = DNAme_delta, fill = RNA.direction)) +
#     geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
#     geom_boxplot(outliers = FALSE, na.rm = TRUE) +
#     scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
#     labs(x = "RNA DE", y = paste0("DNAme delta (",group2," - ",group1,")")) +
#     theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#     facet_grid(rows = vars(Mark), cols = vars(.data[[paste0(day1, "_class_ext")]]))
#   p
#   #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 20, height = 16, units = "cm", create.dir = TRUE)
#   rm(p, tmp_aap)
#
#   ## Export stats for d5mC
#   t_Output <- data.frame("Enh_class" = character(),
#                          "Mark" = character(),
#                          "GroupA" = character(),
#                          "GroupB" = character(),
#                          "A_n" = numeric(),
#                          "B_n" = numeric(),
#                          "A_mean_d5mC" = numeric(),
#                          "B_mean_d5mC" = numeric(),
#                          "A_median_d5mC" = numeric(),
#                          "B_median_d5mC" = numeric(),
#                          "pval" = numeric())
#   for(i in 1:length(Enh_classes)) {
#     for(j in 1:length(Mark_classes)) {
#       ## Pull DNAme_delta values
#       Down <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                               Mark == Mark_classes[j] & RNA.direction == "down" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
#       None <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                               Mark == Mark_classes[j] & RNA.direction == "none" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
#       Up <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] == Enh_classes[i] &
#                                             Mark == Mark_classes[j] & RNA.direction == "up" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
#       ## Test comparisons -- minimum sample size 6
#       if(min(length(Down), length(None), length(Up)) >= 6) {
#         t_DN <- wilcox.test(x = Down, y = None)
#         t_UN <- wilcox.test(x = Up, y = None)
#         t_DU <- wilcox.test(x = Down, y = Up)
#
#         ## Append information to output
#         t_Output <- rbind(t_Output, data.frame("Enh_class" = Enh_classes[i],
#                                                "Mark" = Mark_classes[j],
#                                                "GroupA" = c("down", "up", "down"),
#                                                "GroupB" = c("none", "none", "up"),
#                                                "A_n" = c(length(Down), length(Up), length(Down)),
#                                                "B_n" = c(length(None), length(None), length(Up)),
#                                                "A_mean_d5mC" = round(c(mean(Down), mean(Up), mean(Down)), digits = 2),
#                                                "B_mean_d5mC" = round(c(mean(None), mean(None), mean(Up)), digits = 2),
#                                                "A_median_d5mC" = round(c(median(Down), median(Up), median(Down)), digits = 2),
#                                                "B_median_d5mC" = round(c(median(None), median(None), median(Up)), digits = 2),
#                                                "pval" = c(t_DN$p.value, t_UN$p.value, t_DU$p.value)))
#       }
#     }
#   }
#   ## Clean output
#   t_Output$Enh_class <- factor(t_Output$Enh_class, levels = Enh_classes)
#   t_Output$Mark <- factor(t_Output$Mark, levels = Mark_classes)
#   t_Output$GroupA <- factor(t_Output$GroupA, levels = c("down", "none", "up"))
#   t_Output$GroupB <- factor(t_Output$GroupB, levels = c("down", "none", "up"))
#   t_Output <- t_Output |> arrange(Enh_class, GroupA, GroupB, Mark)
#   t_Output <- t_Output |> relocate(Mark, .after = GroupB)
#   t_Output$paj <- p.adjust(t_Output$pval, method = "bonferroni")
#   t_Output$sig <- t_Output$paj < 0.05
#   ## Round to four significant decimals
#   t_Output$pval <- signif(t_Output$pval, digits = 5)
#   t_Output$paj <- signif(t_Output$paj, digits = 5)
#
#   #cat("Calculating two-sided wilcox test (d5mC ",group2,"-⁠",group1,") per RNA direction per enhancer class\n\n", sep = "", file = Name_log, append = TRUE)
#   #write.table(x = t_Output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
#   #cat("\n", sep = "", file = Name_log, append = TRUE)
#   ## Clean workspace
#   rm(Down, None, Up, t_DN, t_UN, t_DU, t_Output, i, j)
#
#   ## -- Analyze DNAme delta without separating for enhancer class
#   ## Only actv, actvM, and prim enhancers
#   tmp_aap <- Enh_int_pivot_dedup |> filter(.data[[paste0(day1, "_class_ext")]] %in% c("actv", "actvM", "prim"))
#   #Name_plot <- paste0(group1,"_vs_",group2,"_",mark,"_DNAmeDelta_byMark_byRNAseqDE_aap.png")
#   p <- ggplot(tmp_aap, aes(x = RNA.direction, y = DNAme_delta, fill = RNA.direction)) +
#     geom_hline(yintercept = 0, color = "gray0", linewidth = 0.2, linetype = "dashed") +
#     geom_boxplot(outliers = FALSE, na.rm = TRUE) +
#     scale_fill_manual(values = as.vector(c(gen_colors[gen1],"gray70",gen_colors[gen2]))) +
#     labs(x = "RNA DE", y = paste0("DNAme delta (",group2," - ",group1,")")) +
#     theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
#     facet_grid(rows = vars(Mark))
#   p
#   #ggsave(filename = Name_plot, plot = p, device = "png", path = Plot_folder, width = 7, height = 24, units = "cm", create.dir = TRUE)
#   rm(p)
#
#   ## -- Export stats
#   t_Output <- data.frame("Mark" = character(),
#                          "GroupA" = character(),
#                          "GroupB" = character(),
#                          "A_n" = numeric(),
#                          "B_n" = numeric(),
#                          "A_mean_d5mC" = numeric(),
#                          "B_mean_d5mC" = numeric(),
#                          "A_median_d5mC" = numeric(),
#                          "B_median_d5mC" = numeric(),
#                          "pval" = numeric())
#   for(j in 1:length(Mark_classes)) {
#     ## Pull DNAme_delta values
#     Down <- tmp_aap |> filter(Mark == Mark_classes[j] & RNA.direction == "down" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
#     None <- tmp_aap |> filter(Mark == Mark_classes[j] & RNA.direction == "none" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
#     Up <- tmp_aap |> filter(Mark == Mark_classes[j] & RNA.direction == "up" & !is.na(DNAme_delta)) |> pull(DNAme_delta)
#     ## Test comparisons -- minimum sample size 6
#     if(min(length(Down), length(None), length(Up)) >= 6) {
#       t_DN <- wilcox.test(x = Down, y = None)
#       t_UN <- wilcox.test(x = Up, y = None)
#       t_DU <- wilcox.test(x = Down, y = Up)
#
#       ## Append information to output
#       t_Output <- rbind(t_Output, data.frame("Mark" = Mark_classes[j],
#                                              "GroupA" = c("down", "up", "down"),
#                                              "GroupB" = c("none", "none", "up"),
#                                              "A_n" = c(length(Down), length(Up), length(Down)),
#                                              "B_n" = c(length(None), length(None), length(Up)),
#                                              "A_mean_d5mC" = round(c(mean(Down), mean(Up), mean(Down)), digits = 2),
#                                              "B_mean_d5mC" = round(c(mean(None), mean(None), mean(Up)), digits = 2),
#                                              "A_median_d5mC" = round(c(median(Down), median(Up), median(Down)), digits = 2),
#                                              "B_median_d5mC" = round(c(median(None), median(None), median(Up)), digits = 2),
#                                              "pval" = c(t_DN$p.value, t_UN$p.value, t_DU$p.value)))
#     }
#   }
#   ## Clean output
#   t_Output$Mark <- factor(t_Output$Mark, levels = Mark_classes)
#   t_Output$GroupA <- factor(t_Output$GroupA, levels = c("down", "none", "up"))
#   t_Output$GroupB <- factor(t_Output$GroupB, levels = c("down", "none", "up"))
#   t_Output <- t_Output |> arrange(GroupA, GroupB, Mark)
#   t_Output <- t_Output |> relocate(Mark, .after = GroupB)
#   t_Output$paj <- p.adjust(t_Output$pval, method = "bonferroni")
#   t_Output$sig <- t_Output$paj < 0.05
#   ## Round to four significant decimals
#   t_Output$pval <- signif(t_Output$pval, digits = 5)
#   t_Output$paj <- signif(t_Output$paj, digits = 5)
#
#   #cat("Calculating two-sided wilcox test (d5mC ",group2,"-⁠",group1,") per RNA direction (actv+actvM+prim merged)\n\n", sep = "", file = Name_log, append = TRUE)
#   #write.table(x = t_Output, file = Name_log, append = TRUE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
#   #cat("\n", sep = "", file = Name_log, append = TRUE)
#   ## Clean workspace
#   rm(tmp_aap, Down, None, Up, t_DN, t_UN, t_DU, t_Output, j)
# }
#
# rm(Enh_int_pivot, Enh_int_pivot_dedup)
# ## ----
