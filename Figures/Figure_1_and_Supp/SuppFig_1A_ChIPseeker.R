## ChIPseeker annotation of M-ATAC peaks
## Supplemental Figure 1A
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)
library(dplyr)
setwd("/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Data/peaks/MATAC_macs2_broad/")

Current_peak_file <- "WT_MATAC_macs2_broad_D0toD4_merged.bed"

Peaks <- read.table(file = Current_peak_file, header = FALSE, sep = '\t')
colnames(Peaks) <- c("chrom", "start", "end")
Peaks <- makeGRangesFromDataFrame(df = Peaks)

## Annotate without "Downstream (<=300)". --> I tried and it didn't really eliminate it. I guess I can overwrite.
Peaks <- annotatePeak(Peaks, tssRegion=c(-3000, 3000),
                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                      annoDb="org.Mm.eg.db",
                      addFlankGeneInfo=FALSE,
                      #genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic", "Downstream"),
                      verbose = FALSE)
#plotAnnoPie(Peaks)
plotAnnoBar(Peaks)

## ------------------------------
## Custom plotting:
Peaks_df <- as.data.frame(Peaks)

## Simplify annotation label for easier filtering
Peaks_df$annotation_cat <- gsub(pattern = " \\(ENSMUST.*?\\)", replacement = "", x = Peaks_df$annotation)
## Replace "Downstream" with "Distal intergenic"
Peaks_df$annotation_cat[grepl(pattern = "Downstream", x = Peaks_df$annotation_cat)] <- "Distal Intergenic"

table(Peaks_df$annotation_cat)

## Create summary
Annot_summary <- Peaks_df |> group_by(annotation_cat) |> summarise("total" = n(), "perc" = n()/nrow(Peaks_df) * 100)
## Add empty column as x axis
Annot_summary$x <- "x"
## Factorize categories
Annot_summary$annotation_cat <- factor(Annot_summary$annotation_cat, levels = c("Promoter (<=1kb)",
                                                                                "Promoter (1-2kb)",
                                                                                "Promoter (2-3kb)",
                                                                                "5' UTR",
                                                                                "3' UTR",
                                                                                "Exon",
                                                                                "Intron",
                                                                                "Distal Intergenic"))
## Invert levels for plotting
Annot_summary$annotation_cat <- factor(Annot_summary$annotation_cat, levels = rev(levels(Annot_summary$annotation_cat)))

## Plot
p <- ggplot(Annot_summary, aes(x = x, y = perc, fill = annotation_cat)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Promoter (<=1kb)" = "#A6CEE3",
                               "Promoter (1-2kb)" = "#1F78B4",
                               "Promoter (2-3kb)" = "#B2DF8A",
                               "5' UTR" = "#33A02C",
                               "3' UTR" = "#FB9A99",
                               "Exon" = "#E31A1C",
                               "Intron" = "#FF7F00",
                               "Distal Intergenic" = "#FFFF99")) +
  labs(y = "Percentage (%)", fill = "Category") +
  guides(fill = guide_legend(reverse = TRUE)) + coord_flip() +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
p
# ggsave(filename = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Paper_my_figs/Supp_Fig_1/Supp_Fig1A_ChIPseeker_Annot_v2.pdf",
#        plot = p, device = "pdf", width = 14, height = 7, units = "cm")


## Re-order levels and format table
Annot_summary$annotation_cat <- factor(Annot_summary$annotation_cat, levels = rev(levels(Annot_summary$annotation_cat)))
Annot_summary <- Annot_summary |> arrange(annotation_cat)
Annot_summary <- Annot_summary[,c("annotation_cat", "total", "perc")]
colnames(Annot_summary) <- c("Category", "N_peaks", "Perc_peaks")
Annot_summary$Perc_peaks <- round(Annot_summary$Perc_peaks, digits = 2)

## Export statistics
#write.table(x = Annot_summary, file = "/Volumes/MMF_IJM/Greenberg_Lab/M-ATAC_Project/Paper_my_figs/Supp_Fig_1/Info.txt", append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)

