## Figure 2A: Heatmaps D4

cd /scratch/Marlet/

## Accessibility - D4
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D4_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/TKO_D4_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/DKO_D4_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-TKO-DKO_D4_MATACsignal_K5_sorted.gz \
--numberOfProcessors 16

plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-TKO-DKO_D4_MATACsignal_K5_sorted.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/Plots_same_scale/WT-TKO-DKO_D4_MATACsignal_K5_sorted.pdf \
--sortRegions keep \
--whatToShow "heatmap and colorbar" \
--zMin 0 --zMax 1.0 --colorMap RdBu_r \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D4" "TKO-D4" "DKO-D4" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"

## ---------------------

## Accessibility - D2
sbatch -n 8 computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D2_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/TKO_D2_MATAC_mergedreps34_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/DKO_D2_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-TKO-DKO_D2_MATACsignal_K5_sorted.gz \
--numberOfProcessors 8

plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-TKO-DKO_D2_MATACsignal_K5_sorted.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/Plots_same_scale/WT-TKO-DKO_D2_MATACsignal_K5_sorted.pdf \
--sortRegions keep \
--whatToShow "heatmap and colorbar" \
--zMin 0 --zMax 1.0 --colorMap RdBu_r \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D2" "TKO-D2" "DKO-D2" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"

## ---------------------

## Accessibility - D0
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/WT_D0_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/TKO_D0_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/DKO_D0_MATAC_mergedreps_trimV6_mm10_rmDup_q10_b1_s0_CPM.bw \
--binSize 25 \
--sortRegions keep \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-TKO-DKO_D0_MATACsignal_K5_sorted.gz \
--numberOfProcessors 16

plotHeatmap -m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/WT-TKO-DKO_D0_MATACsignal_K5_sorted.gz \
-o /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/Plots/Plots_same_scale/WT-TKO-DKO_D0_MATACsignal_K5_sorted.pdf \
--sortRegions keep \
--whatToShow "heatmap and colorbar" \
--zMin 0 --zMax 1.0 --colorMap RdBu_r \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "WT-D0" "TKO-D0" "DKO-D0" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"

