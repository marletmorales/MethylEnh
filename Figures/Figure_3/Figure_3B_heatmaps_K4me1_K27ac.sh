## H3K4me1 and H3K27ac heatmaps
## Figure 2B, Supplemental Figure 3A-C, and Supplemental Figure 4A-C

cd /scratch/Marlet/

## --------------------- ##
## -- K4me1 - Dynamic -- ##
## --------------------- ##

## -- WT K4me1 -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/WT_D0-D2-D4_K4me1signal_DS_K5_ud500bp_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.8 \
--colorMap YlOrBr \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K4me1" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"


## -- DKO K4me1 -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D0_K4me1_MCnT_mergedreps_250819_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D4_K4me1_MCnT_mergedreps_250819_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/DKO_D0-D2-D4_K4me1signal_DS_K5_ud500bp_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.8 \
--colorMap YlOrBr \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K4me1" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "DKO-D0" "DKO-D2" "DKO-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"


## -- TKO K4me1 -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/TKO_D0-D2-D4_K4me1signal_DS_K5_ud500bp_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.8 \
--colorMap YlOrBr \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K4me1" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "TKO-D0" "TKO-D2" "TKO-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"


## -- Plot profiles
## WT K4me1
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/WT_D0-D2-D4_K4me1signal_DS_K5_ud500bp_sorted_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 0.9 --yAxisLabel "K4me1" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--startLabel "5'" --endLabel "3'"

## DKO K4me1
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/DKO_D0-D2-D4_K4me1signal_DS_K5_ud500bp_sorted_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 0.9 --yAxisLabel "K4me1" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "DKO-D0" "DKO-D2" "DKO-D4" \
--startLabel "5'" --endLabel "3'"

## TKO K4me1
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K4me1signal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/TKO_D0-D2-D4_K4me1signal_DS_K5_ud500bp_sorted_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 0.9 --yAxisLabel "K4me1" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "TKO-D0" "TKO-D2" "TKO-D4" \
--startLabel "5'" --endLabel "3'"


## -------------------- ##
## -- K4me1 - Static -- ##
## -------------------- ##

## -- WT K4me1 static -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K4me1signal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K4me1signal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/WT_D0-D2-D4_K4me1signal_DS_K5_ud500bp_static.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.8 \
--colorMap YlOrBr \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K4me1" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"


## -- DKO K4me1 static -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D0_K4me1_MCnT_mergedreps_250819_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D4_K4me1_MCnT_mergedreps_250819_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K4me1signal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K4me1signal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/DKO_D0-D2-D4_K4me1signal_DS_K5_ud500bp_static.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.8 \
--colorMap YlOrBr \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K4me1" \
--samplesLabel "DKO-D0" "DKO-D2" "DKO-D4" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"


## -- TKO K4me1 static -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K4me1signal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K4me1signal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/TKO_D0-D2-D4_K4me1signal_DS_K5_ud500bp_static.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.8 \
--colorMap YlOrBr \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K4me1" \
--samplesLabel "TKO-D0" "TKO-D2" "TKO-D4" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"


## -- K4me1 static profiles: Plot three genotypes together

## D0 (WT, DKO, TKO)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D0_K4me1_MCnT_mergedreps_250819_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D0_K4me1signal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## D2 (WT, DKO, TKO)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D2_K4me1signal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## D4 (WT, DKO, TKO)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D4_K4me1_MCnT_mergedreps_250819_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D4_K4me1signal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## -- Plot profiles
## D0 (WT, DKO, TKO)
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D0_K4me1signal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/WT-DKO-TKO_D0_K4me1signal_DS_K5_ud500bp_static_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--perGroup \
--colors "#266799" "#4C4C4C" "#F68532" \
--yMin 0 --yMax 0.9 --yAxisLabel "K4me1" \
--regionsLabel "static" \
--samplesLabel "WT-D0" "DKO-D0" "TKO-D0" \
--startLabel "5'" --endLabel "3'"

## D2 (WT, DKO, TKO)
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D2_K4me1signal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/WT-DKO-TKO_D2_K4me1signal_DS_K5_ud500bp_static_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--perGroup \
--colors "#266799" "#4C4C4C" "#F68532" \
--yMin 0 --yMax 0.9 --yAxisLabel "K4me1" \
--regionsLabel "static" \
--samplesLabel "WT-D2" "DKO-D2" "TKO-D2" \
--startLabel "5'" --endLabel "3'"

## D4 (WT, DKO, TKO)
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D4_K4me1signal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K4me1/WT-DKO-TKO_D4_K4me1signal_DS_K5_ud500bp_static_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--perGroup \
--colors "#266799" "#4C4C4C" "#F68532" \
--yMin 0 --yMax 0.9 --yAxisLabel "K4me1" \
--regionsLabel "static" \
--samplesLabel "WT-D4" "DKO-D4" "TKO-D4" \
--startLabel "5'" --endLabel "3'"


## -----------------------------------------------------
## -----------------------------------------------------

## --------------------- ##
## -- K27ac - Dynamic -- ##
## --------------------- ##

## -- WT K27ac -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/WT_D0-D2-D4_K27acsignal_DS_K5_ud500bp_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.6 \
--colorMap YlGn \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K27ac" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"


## -- DKO K27ac -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/DKO_D0-D2-D4_K27acsignal_DS_K5_ud500bp_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.6 \
--colorMap YlGn \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K27ac" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "DKO-D0" "DKO-D2" "DKO-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"


## -- TKO K27ac -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/TKO_D0-D2-D4_K27acsignal_DS_K5_ud500bp_sorted.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.6 \
--colorMap YlGn \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K27ac" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "TKO-D0" "TKO-D2" "TKO-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"


## -- Plot profiles
## WT K27ac
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/WT_D0-D2-D4_K27acsignal_DS_K5_ud500bp_sorted_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 1.1 --yAxisLabel "K27ac" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--startLabel "5'" --endLabel "3'"

## DKO K27ac
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/DKO_D0-D2-D4_K27acsignal_DS_K5_ud500bp_sorted_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 1.1 --yAxisLabel "K27ac" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "DKO-D0" "DKO-D2" "DKO-D4" \
--startLabel "5'" --endLabel "3'"

## TKO K27ac
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K27acsignal_DS_K5_sorted_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/TKO_D0-D2-D4_K27acsignal_DS_K5_ud500bp_sorted_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 1.1 --yAxisLabel "K27ac" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "TKO-D0" "TKO-D2" "TKO-D4" \
--startLabel "5'" --endLabel "3'"


## -------------------- ##
## -- K27ac - Static -- ##
## -------------------- ##

## -- WT K27ac static -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K27acsignal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT_D0-D2-D4_K27acsignal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/WT_D0-D2-D4_K27acsignal_DS_K5_ud500bp_static.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.6 \
--colorMap YlGn \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K27ac" \
--samplesLabel "WT-D0" "WT-D2" "WT-D4" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"


## -- DKO K27ac static -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K27acsignal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_K27acsignal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/DKO_D0-D2-D4_K27acsignal_DS_K5_ud500bp_static.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.6 \
--colorMap YlGn \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K27ac" \
--samplesLabel "DKO-D0" "DKO-D2" "DKO-D4" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"


## -- TKO K27ac static -- #
## Calculate matrix
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K27acsignal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## Plot keeping original sorting
plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/TKO_D0-D2-D4_K27acsignal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/TKO_D0-D2-D4_K27acsignal_DS_K5_ud500bp_static.pdf \
--sortRegions keep \
--zMin 0 --zMax 0.6 \
--colorMap YlGn \
--whatToShow "heatmap and colorbar" \
--yAxisLabel "K27ac" \
--samplesLabel "TKO-D0" "TKO-D2" "TKO-D4" \
--xAxisLabel "NonProm static peak" --startLabel "5'" --endLabel "3'"


## -- K27ac static profiles: Plot three genotypes together

## D0 (WT, DKO, TKO)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D0_K27acsignal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## D2 (WT, DKO, TKO)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D2_K27acsignal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## D4 (WT, DKO, TKO)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Static_peaks_sorted_Descend_WT_D0-D2-D4.bed3 \
--upstream 500 \
--downstream 500 \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/DKO_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/TKO_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D4_K27acsignal_DS_K5_static_ud500bp.gz \
--numberOfProcessors 12

## -- Plot profiles
## D0 (WT, DKO, TKO)
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D0_K27acsignal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/WT-DKO-TKO_D0_K27acsignal_DS_K5_ud500bp_static_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--perGroup \
--colors "#266799" "#4C4C4C" "#F68532" \
--yMin 0 --yMax 1.1 --yAxisLabel "K27ac" \
--regionsLabel "static" \
--samplesLabel "WT-D0" "DKO-D0" "TKO-D0" \
--startLabel "5'" --endLabel "3'"

## D2 (WT, DKO, TKO)
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D2_K27acsignal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/WT-DKO-TKO_D2_K27acsignal_DS_K5_ud500bp_static_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--perGroup \
--colors "#266799" "#4C4C4C" "#F68532" \
--yMin 0 --yMax 1.1 --yAxisLabel "K27ac" \
--regionsLabel "static" \
--samplesLabel "WT-D2" "DKO-D2" "TKO-D2" \
--startLabel "5'" --endLabel "3'"

## D4 (WT, DKO, TKO)
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-DKO-TKO_D4_K27acsignal_DS_K5_static_ud500bp.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/K27ac/WT-DKO-TKO_D4_K27acsignal_DS_K5_ud500bp_static_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--perGroup \
--colors "#266799" "#4C4C4C" "#F68532" \
--yMin 0 --yMax 1.1 --yAxisLabel "K27ac" \
--regionsLabel "static" \
--samplesLabel "WT-D4" "DKO-D4" "TKO-D4" \
--startLabel "5'" --endLabel "3'"

