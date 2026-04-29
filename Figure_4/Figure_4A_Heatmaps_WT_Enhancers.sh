## ---- Enhancers_v6 using DS libraries for visualization ---- ##
## See 2025_12_18_epipax_command_log.txt for code


## NEED TO UPDATE:
## WT enhancers ('actv', 'prim', 'inac') for D0, D2, D4

## ---------------------
## -- Calculate Matrix

## D0
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D0_actv.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D0_prim.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D0_inac.bed \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_MATAC_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D0_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--upstream 500 \
--downstream 500 \
--regionBodyLength 1000 \
--sortRegions descend \
--sortUsingSamples 1 \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D0_ud500bp_v2.gz \
--outFileSortedRegions /scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D0_sorted_actv-prim-inac_K4me1_ud500bp.bed \
--numberOfProcessors 16
# Submitted batch job 32788

## D2
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D2_actv.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D2_prim.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D2_inac.bed \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_MATAC_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D2_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--upstream 500 \
--downstream 500 \
--regionBodyLength 1000 \
--sortRegions descend \
--sortUsingSamples 1 \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D2_ud500bp_v2.gz \
--outFileSortedRegions /scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D2_sorted_actv-prim-inac_K4me1_ud500bp.bed \
--numberOfProcessors 16
# Submitted batch job 32789

## D4
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D4_actv.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D4_prim.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D4_inac.bed \
-S \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_K4me1_MCnT_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_MATAC_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
/data/marlet_data/epimethyltag/Downsample/bigWigs/WT_D4_MCnT-EMseq_K27ac_mergedreps_trimV6_bsmrk_mm10_mrkdup_DS_kpDup_q10_b1_s0_CPM.bw \
--upstream 500 \
--downstream 500 \
--regionBodyLength 1000 \
--sortRegions descend \
--sortUsingSamples 1 \
--missingDataAsZero \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D4_ud500bp_v2.gz \
--outFileSortedRegions /scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v5/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v5_D4_sorted_actv-prim-inac_K4me1_ud500bp.bed \
--numberOfProcessors 16
# Submitted batch job 32790


## -----------------
## -- Plot heatmap

## D0
plotHeatmap -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D0_ud500bp_v2.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Plots/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D0_ud500bp_v2.pdf \
--sortRegions keep \
--colorMap YlOrBr RdBu_r YlGn \
--zMin 0 0.15 0 --zMax 1.6 0.8 0.8 \
--whatToShow "heatmap and colorbar" \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D0_K4me1" "WT-D0_MATAC" "WT-D0_K27ac" \
--xAxisLabel "Enhancers D0" --startLabel "5'" --endLabel "3'"

## D2
plotHeatmap -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D2_ud500bp_v2.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Plots/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D2_ud500bp_v2.pdf \
--sortRegions keep \
--colorMap YlOrBr RdBu_r YlGn \
--zMin 0 0.15 0 --zMax 1.6 0.8 0.8 \
--whatToShow "heatmap and colorbar" \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D2_K4me1" "WT-D2_MATAC" "WT-D2_K27ac" \
--xAxisLabel "Enhancers D2" --startLabel "5'" --endLabel "3'"

## D4
plotHeatmap -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D4_ud500bp_v2.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Plots/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D4_ud500bp_v2.pdf \
--sortRegions keep \
--colorMap YlOrBr RdBu_r YlGn \
--zMin 0 0.15 0 --zMax 1.6 0.8 0.8 \
--whatToShow "heatmap and colorbar" \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D4_K4me1" "WT-D4_MATAC" "WT-D4_K27ac" \
--xAxisLabel "Enhancers D4" --startLabel "5'" --endLabel "3'"


## -----------------
## -- Plot profile

## D0
plotProfile \
-m /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D0_ud500bp_v2.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Plots/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D0_ud500bp_v2_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E5614E" "#F0AC4B" "#77BEF0" \
--yAxisLabel "Enhancers" \
--yMin 0 0 0 --yMax 1.8 0.7 2.0 \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D0_K4me1" "WT-D0_MATAC" "WT-D0_K27ac" \
--startLabel "5'" --endLabel "3'"

## D2
plotProfile \
-m /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D2_ud500bp_v2.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Plots/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D2_ud500bp_v2_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E5614E" "#F0AC4B" "#77BEF0" \
--yAxisLabel "Enhancers" \
--yMin 0 0 0 --yMax 1.8 0.7 2.0 \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D2_K4me1" "WT-D2_MATAC" "WT-D2_K27ac" \
--startLabel "5'" --endLabel "3'"

## D4
plotProfile \
-m /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Matrix/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D4_ud500bp_v2.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v5/Plots/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v5_D4_ud500bp_v2_Profile.pdf \
--plotType lines --plotHeight 6 --plotWidth 6 \
--colors "#E5614E" "#F0AC4B" "#77BEF0" \
--yAxisLabel "Enhancers" \
--yMin 0 0 0 --yMax 1.8 0.7 2.0 \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D4_K4me1" "WT-D4_MATAC" "WT-D4_K27ac" \
--startLabel "5'" --endLabel "3'"
