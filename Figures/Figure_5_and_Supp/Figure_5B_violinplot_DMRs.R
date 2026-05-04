## 5hmC in hMATAC peaks; separate by DMR intersection
## Figure 5B
library(ggplot2)
library(dplyr)
library(bedtoolsr)
library(stringr)
library(annotatr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/Measure_DNAme/DNAme_ROIs/WT_hMATAC_peakome_RmProm/")

## ----------------------------- ##
## ---- Annotate 5hmC peaks ---- ##
## ----------------------------- ##

## --------
## Make annotatr object from peakome
WT_hMATAC <- read_regions(con = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Data/peaks/hMATAC/WT_D4_hMATAC_mergedreps_2507.filtered.bed", genome = "mm10", format = "bed")
## Build annotations
Annotations <- build_annotations(genome = 'mm10', annotations = c("mm10_basicgenes", "mm10_cpgs"))

## Annotate WT_hMATAC peaks
WT_hMATAC <- annotate_regions(regions = WT_hMATAC, annotations = Annotations, ignore.strand = TRUE)
# Coerce to data.frame
WT_hMATAC <- data.frame(WT_hMATAC)
## Add coordinate
WT_hMATAC$coord <- str_c(str_c(WT_hMATAC$seqnames, WT_hMATAC$start, sep = ":"), WT_hMATAC$end, sep = "-")
## Keep only annot.type to remove duplicates (i.e. if a peak maps to multiple intron ids, it will only be counted once)
WT_hMATAC <- WT_hMATAC |> dplyr::select(!c(annot.seqnames,
                                           annot.start,
                                           annot.end,
                                           annot.width,
                                           annot.strand,
                                           annot.id,
                                           annot.tx_id,
                                           annot.gene_id,
                                           annot.symbol)) |> distinct()

## Prioritize annotation; keep only one row per peak
unique(WT_hMATAC$annot.type)

## Order (outside of promoters):
Order_priority <- c("mm10_cpg_islands",
                    "mm10_cpg_shores",
                    "mm10_genes_5UTRs",
                    "mm10_genes_3UTRs",
                    "mm10_genes_exons",
                    "mm10_genes_introns",
                    "mm10_genes_1to5kb",
                    "mm10_cpg_shelves",
                    "mm10_cpg_inter")

## Make function for prioritization; it also divides promoters into CGI and not-CGI.
prioritize_annot <- function(Sub_coord, Order_priority) {
  ## Input data.frame of multiple rows (Sub_coord), return data.frame of one row (Output_coord)
  ## If it is a promoter
  if("mm10_genes_promoters" %in% Sub_coord$annot.type) {
    ## If it is a CGI
    if("mm10_cpg_islands" %in% Sub_coord$annot.type) {
      new_annot <- "mm10_genes_promoters_cgi"
    } else {
      new_annot <- "mm10_genes_promoters_noncgi"
    }
    ## Return single row
    Output_coord <- Sub_coord |> filter(annot.type == "mm10_genes_promoters")
    Output_coord$annot.type <- new_annot
    return(Output_coord)
  } else {
    ## If it's not a promoter, follow regular prioritization
    Output_coord <- Sub_coord |> arrange(factor(annot.type, levels = Order_priority))
    Output_coord <- Output_coord[1,]
    return(Output_coord)
  }
}

## Run prioritization function over WT_hMATAC_annotated_simple
Coords <- unique(WT_hMATAC$coord)
WT_hMATAC_one_annot <- WT_hMATAC[0,]

## It takes ~5 min.
for(i in 1:length(Coords)){
  Sub_coord <- WT_hMATAC |> filter(coord == Coords[i])
  CpG_onerow <- prioritize_annot(Sub_coord, Order_priority)
  WT_hMATAC_one_annot <- rbind(WT_hMATAC_one_annot, CpG_onerow)
}
rm(Coords, i, Sub_coord, CpG_onerow)

## Factorize annotation type
WT_hMATAC_one_annot$annot.type <- factor(WT_hMATAC_one_annot$annot.type, levels = c("mm10_genes_promoters_cgi", "mm10_genes_promoters_noncgi", Order_priority))

## Clean
WT_hMATAC_one_annot <- WT_hMATAC_one_annot |> arrange(seqnames, start, end)
colnames(WT_hMATAC_one_annot) <- c("chr", "start", "end", "width", "strand", "annot.type", "coord")
## Correct start position (shifted by +1 by GRanges)
WT_hMATAC_one_annot$start <- WT_hMATAC_one_annot$start - 1
## Re-do coordinates
WT_hMATAC_one_annot$coord <- str_c(str_c(WT_hMATAC_one_annot$chr, WT_hMATAC_one_annot$start, sep = ":"), WT_hMATAC_one_annot$end, sep = "-")

# ## Export
# write.table(x = WT_hMATAC_one_annot, append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE,
#             file = "hMATAC/WT_D4_hMATAC_mergedreps_2507.filtered_annotated.bed")

rm(WT_hMATAC, Order_priority, prioritize_annot, Annotations)
## --------

## --------------------------------


## --------------------------------------- ##
## ---- Incorporate 5hmC measurements ---- ##
## --------------------------------------- ##

## If continuing from above, rename variable
hMATAC_annotation <- WT_hMATAC_one_annot
rm(WT_hMATAC_one_annot)

# ## If jumping from outside, read table
# hMATAC_annotation <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Data/peaks/hMATAC/WT_D4_hMATAC_mergedreps_2507.filtered_annotated.bed", header = TRUE, sep = '\t')

## -- Read 5hmC measurements
WT_hMATAC <- read.table(file = "meanDNAme_WT_D4_hMATAC_mergedreps__Coordinates_WT_D4_hMATAC_mergedreps_2507_filteredBroadPeak__2025_12_10.txt", header = TRUE, sep = '\t')
## Add coord
WT_hMATAC$coord <- str_c(str_c(WT_hMATAC$chr, WT_hMATAC$start, sep = ":"), WT_hMATAC$end, sep = "-")

## Merge annotation with measurements
WT_hMATAC <- merge(x = WT_hMATAC, y = hMATAC_annotation, all.x = TRUE, all.y = FALSE)
## Clean
WT_hMATAC <- WT_hMATAC |> dplyr::select(!c(width, strand))
WT_hMATAC <- WT_hMATAC |> relocate(coord, .after = last_col())
WT_hMATAC <- WT_hMATAC |> arrange(chr, start, end)
rm(hMATAC_annotation)

## --------------------------------

## -- Read DMRs
DMR_WGBS <- read.table(file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Results/MethyLasso/MethyLasso_WTvsTet12DKO_D4_WGBS_v2/WT_D4_vs_DKO_D4_dmrs_annot.tsv", header = TRUE, sep = '\t')

## -- Annotate 5hmC peaks according to whether they intersect with a DMR or not
WT_hMATAC_i_DMR <- bt.intersect(a = WT_hMATAC, b = DMR_WGBS, c = TRUE)
colnames(WT_hMATAC_i_DMR) <- c(colnames(WT_hMATAC), 'Num_ov')

## Add additional column on whether it overlaps DMR or not
WT_hMATAC_i_DMR <- WT_hMATAC_i_DMR |> mutate("DMR_intersect" = if_else(condition = Num_ov > 0,
                                                                       true = "In_DMR",
                                                                       false = "Out_DMR"))
WT_hMATAC_i_DMR <- WT_hMATAC_i_DMR[, colnames(WT_hMATAC_i_DMR) != "Num_ov"]

## Factorize DMR intersect
WT_hMATAC_i_DMR$DMR_intersect <- factor(WT_hMATAC_i_DMR$DMR_intersect, levels = c("Out_DMR", "In_DMR"))

## Add a column for grouping annotations
WT_hMATAC_i_DMR <- WT_hMATAC_i_DMR |> mutate("annot.simple" = case_when(annot.type == "mm10_genes_promoters_cgi" ~ "prom_cgi",
                                                                        annot.type == "mm10_genes_promoters_noncgi" ~ "prom_noncgi",
                                                                        .default = "Non_prom"))
WT_hMATAC_i_DMR$annot.simple <- factor(WT_hMATAC_i_DMR$annot.simple, levels = c("prom_cgi", "prom_noncgi", "Non_prom"))


## --------------------------------

## ------------------------- ##
## ---- Plots and Stats ---- ##
## ------------------------- ##

my_theme <- theme_bw() + theme(legend.position = "none",
                               panel.grid = element_blank(),
                               plot.title = element_text(hjust = 0.5, size = 9.5),
                               axis.title = element_text(size = 9),
                               axis.text.x = element_text(size = 7.5),
                               axis.title.x = element_blank())

## CGI promoter
p_promCGI <- ggplot(WT_hMATAC_i_DMR[WT_hMATAC_i_DMR$annot.simple == "prom_cgi",],
                    aes(x = DMR_intersect, y = meth.wmean, fill = DMR_intersect)) +
  geom_violin(quantile.linetype = "solid", bounds = c(0,100)) +
  scale_fill_manual(values = c("gray65","#BD327F")) +
  labs(y = "5hmC %", title = "CGI promoter") +
  coord_cartesian(ylim = c(0, 8)) +
  my_theme
p_promCGI

## Non-CGI promoter
p_promNonCGI <- ggplot(WT_hMATAC_i_DMR[WT_hMATAC_i_DMR$annot.simple == "prom_noncgi",],
                       aes(x = DMR_intersect, y = meth.wmean, fill = DMR_intersect)) +
  geom_violin(quantile.linetype = "solid", bounds = c(0,100)) +
  scale_fill_manual(values = c("gray65","#BD327F")) +
  labs(y = "5hmC %", title = "Non-CGI promoter") +
  coord_cartesian(ylim = c(0, 40)) +
  my_theme
p_promNonCGI

## Outside of promoters
p_NonProm <- ggplot(WT_hMATAC_i_DMR[WT_hMATAC_i_DMR$annot.simple == "Non_prom",],
                    aes(x = DMR_intersect, y = meth.wmean, fill = DMR_intersect)) +
  geom_violin(quantile.linetype = "solid", bounds = c(0,100)) +
  scale_fill_manual(values = c("gray65","#BD327F")) +
  labs(y = "5hmC %", title = "Non-promoter") +
  coord_cartesian(ylim = c(0, 40)) +
  my_theme
p_NonProm

## Combine
p_Master <- gridExtra::grid.arrange(p_promCGI, p_promNonCGI, p_NonProm, nrow = 1)
# ggsave(filename = "WT_D4_hMATAC_5hmC_in-out_DMRs_byCat.pdf", plot = p_Master, device = "pdf", width = 12, height = 6, units = "cm")

## --------------------------------

## -- Summary
Summary_hMATAC <- WT_hMATAC_i_DMR |> group_by(annot.simple, DMR_intersect) |> summarise("Count" = n(),
                                                                                        "Mean_5hmC" = round(mean(meth.wmean), digits = 4),
                                                                                        "Q25_5hmC" = round(quantile(meth.wmean, 0.25), digits = 4),
                                                                                        "Q50_5hmC" = round(quantile(meth.wmean, 0.50), digits = 4),
                                                                                        "Q75_5hmC" = round(quantile(meth.wmean, 0.75), digits = 4))
Summary_hMATAC

## -- Wilcox testing

## CGI promoter
Vals_OutDMR <- WT_hMATAC_i_DMR |> filter(annot.simple == "prom_cgi" & DMR_intersect == "Out_DMR") |> pull(meth.wmean)
Vals_InDMR <- WT_hMATAC_i_DMR |> filter(annot.simple == "prom_cgi" & DMR_intersect == "In_DMR") |> pull(meth.wmean)
T_Out_In <- wilcox.test(x = Vals_OutDMR, y = Vals_InDMR, alternative = "less")
T_results <- data.frame("annot.simple" = "prom_cgi",
                        "DMR_intersect" = c("Out_DMR", "In_DMR"),
                        "p_val" = c(NA, T_Out_In$p.value))

## Non-CGI promoter
Vals_OutDMR <- WT_hMATAC_i_DMR |> filter(annot.simple == "prom_noncgi" & DMR_intersect == "Out_DMR") |> pull(meth.wmean)
Vals_InDMR <- WT_hMATAC_i_DMR |> filter(annot.simple == "prom_noncgi" & DMR_intersect == "In_DMR") |> pull(meth.wmean)
T_Out_In <- wilcox.test(x = Vals_OutDMR, y = Vals_InDMR, alternative = "less")
T_results <- rbind(T_results, data.frame("annot.simple" = "prom_noncgi",
                                         "DMR_intersect" = c("Out_DMR", "In_DMR"),
                                         "p_val" = c(NA, T_Out_In$p.value)))

## Outside of promoters
Vals_OutDMR <- WT_hMATAC_i_DMR |> filter(annot.simple == "Non_prom" & DMR_intersect == "Out_DMR") |> pull(meth.wmean)
Vals_InDMR <- WT_hMATAC_i_DMR |> filter(annot.simple == "Non_prom" & DMR_intersect == "In_DMR") |> pull(meth.wmean)
T_Out_In <- wilcox.test(x = Vals_OutDMR, y = Vals_InDMR, alternative = "less")
T_results <- rbind(T_results, data.frame("annot.simple" = "Non_prom",
                                         "DMR_intersect" = c("Out_DMR", "In_DMR"),
                                         "p_val" = c(NA, T_Out_In$p.value)))

## Trim p-value digits
T_results$p_val <- scales::scientific(T_results$p_val, digits = 5)

## Append results to summary stats
Summary_hMATAC <- merge(Summary_hMATAC, T_results)
## Sort
Summary_hMATAC <- Summary_hMATAC |> arrange(annot.simple, DMR_intersect)

## Export
# write.table(x = Summary_hMATAC, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE,
#             file = "Violin_plot_hMATAC_in_DMRs/Info.txt")

