## WT heatmaps: M-ATAC signal 
## Figure 1C and Supplemental Figure 1B,D

## -- Plot DNA accessibility and methylation with ordered clusters
cd /scratch/Marlet/

## ----------------------------------------------
## ------------- Dynamic (DA) peaks -------------

## Accessibility
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D0_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D2_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D4_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATACsignal_K5_sorted.gz \
--numberOfProcessors 12

plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATACsignal_K5_sorted.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_1/Panel_tC/WT_D0-D2-D4_MATACsignal_K5_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 1.0 --colorMap RdBu_r \
--yMin 0 --yMax 0.9 --yAxisLabel "Accessibility" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"

## Methylation
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D0_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D2_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D4_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_D0_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_EpiD2_WT_H3_doxPos_WGBS_AD_067_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_D4_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
--upstream 500 \
--downstream 500 \
--regionBodyLength 500 \
--binSize 50 \
--sortRegions keep \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted.gz \
--numberOfProcessors 12

plotHeatmap \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_1/Supp_Heatmap_5mC/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_Heatmap.pdf \
--sortRegions keep --interpolationMethod nearest \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 28 \
--heatmapWidth 4 \
--zMin 0 --zMax 100 --colorMap RdPu --missingDataColor silver \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0_ATAC" "WT-D2_ATAC" "WT-D4_ATAC" "WT-D0_WGBS" "WT-D2_WGBS" "WT-D4_WGBS" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"

plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_1/Supp_Heatmap_5mC/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_Profile.pdf \
--plotHeight 6 \
--plotWidth 5 \
--plotType lines \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 100 --yAxisLabel "5mC" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0_ATAC" "WT-D2_ATAC" "WT-D4_ATAC" "WT-D0_WGBS" "WT-D2_WGBS" "WT-D4_WGBS" \
--startLabel "5'" --endLabel "3'"


## ----------------------------------------
## ------------- Static peaks -------------

## Accessibility
sbatch computeMatrix scale-regions \
-R /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D0_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D2_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D4_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions ascend \
--sortUsingSamples 1 \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATACsignal_K5_static.gz \
--numberOfProcessors 12

plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATACsignal_K5_static.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_1/Panel_tC/SUPP_WT_D0-D2-D4_MATACsignal_static_regions.pdf \
--colorMap RdBu_r \
--sortRegions descend \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--xAxisLabel "NonProm static peak" \
--zMin 0 --zMax 1.0 --colorMap RdBu_r \
--yMin 0 --yMax 0.9 --yAxisLabel "Accessibility" \
--startLabel "5'" \
--endLabel "3'" \
--outFileSortedRegions /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3

## Methylation
sbatch computeMatrix scale-regions \
-R /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D0_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D2_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D4_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_D0_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_EpiD2_WT_H3_doxPos_WGBS_AD_067_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_D4_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
--upstream 500 \
--downstream 500 \
--regionBodyLength 500 \
--binSize 50 \
--sortRegions keep \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_static.gz \
--numberOfProcessors 12

plotHeatmap \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_static.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_1/Supp_Heatmap_5mC/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_static_Heatmap.pdf \
--sortRegions keep --interpolationMethod nearest \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 28 \
--heatmapWidth 4 \
--zMin 0 --zMax 100 --colorMap RdPu --missingDataColor silver \
--samplesLabel "WT-D0_ATAC" "WT-D2_ATAC" "WT-D4_ATAC" "WT-D0_WGBS" "WT-D2_WGBS" "WT-D4_WGBS" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"

plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_static.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_1/Supp_Heatmap_5mC/WT_D0-D2-D4_MATAC-WGBS_5mC_K5_static_Profile.pdf \
--plotHeight 6 \
--plotWidth 5 \
--plotType lines \
--colors "#266799" "#266799" "#266799" "#266799" "#266799" "#266799" \
--yMin 0 --yMax 100 --yAxisLabel "5mC" \
--samplesLabel "WT-D0_ATAC" "WT-D2_ATAC" "WT-D4_ATAC" "WT-D0_WGBS" "WT-D2_WGBS" "WT-D4_WGBS" \
--startLabel "5'" --endLabel "3'"
