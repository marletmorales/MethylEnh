## WT enhancers: heatmaps for H3K4me1, ATAC, and H3K27ac
## Figure 4A and Supplemental Figure 5A

## ---------------------
## -- Calculate Matrix

## D0 (mESCs)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D0_actv.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D0_prim.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D0_inac.bed \
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
--outFileName /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D0_ud500bp.gz \
--outFileSortedRegions /scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D0_sorted_actv-prim-inac_K4me1_ud500bp.bed \
--numberOfProcessors 16

## D2 (EpiLCs)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D2_actv.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D2_prim.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D2_inac.bed \
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
--outFileName /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D2_ud500bp.gz \
--outFileSortedRegions /scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D2_sorted_actv-prim-inac_K4me1_ud500bp.bed \
--numberOfProcessors 16

## D4 (ExEpiLCs)
sbatch computeMatrix scale-regions \
-R \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D4_actv.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D4_prim.bed \
/scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D4_inac.bed \
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
--outFileName /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D4_ud500bp.gz \
--outFileSortedRegions /scratch/Marlet/peaks/macs2_broad/Enhancers_intersect_v6/WT_Enhancer_peakome/WT_Master_Enhancer_Peakome_v6_D4_sorted_actv-prim-inac_K4me1_ud500bp.bed \
--numberOfProcessors 16


## -----------------
## -- Plot heatmap

## D0 (mESCs)
plotHeatmap -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D0_ud500bp.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Plots/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D0_ud500bp.pdf \
--sortRegions keep \
--colorMap YlOrBr RdBu_r YlGn \
--zMin 0 0.15 0 --zMax 1.6 0.8 0.8 \
--whatToShow "heatmap and colorbar" \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D0_K4me1" "WT-D0_MATAC" "WT-D0_K27ac" \
--xAxisLabel "Enhancers D0" --startLabel "5'" --endLabel "3'"

## D2 (EpiLCs)
plotHeatmap -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D2_ud500bp.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Plots/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D2_ud500bp.pdf \
--sortRegions keep \
--colorMap YlOrBr RdBu_r YlGn \
--zMin 0 0.15 0 --zMax 1.6 0.8 0.8 \
--whatToShow "heatmap and colorbar" \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D2_K4me1" "WT-D2_MATAC" "WT-D2_K27ac" \
--xAxisLabel "Enhancers D2" --startLabel "5'" --endLabel "3'"

## D4 (ExEpiLCs)
plotHeatmap -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D4_ud500bp.gz \
-o /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Plots/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D4_ud500bp.pdf \
--sortRegions keep \
--colorMap YlOrBr RdBu_r YlGn \
--zMin 0 0.15 0 --zMax 1.6 0.8 0.8 \
--whatToShow "heatmap and colorbar" \
--regionsLabel "actv" "prim" "inac" \
--samplesLabel "WT-D4_K4me1" "WT-D4_MATAC" "WT-D4_K27ac" \
--xAxisLabel "Enhancers D4" --startLabel "5'" --endLabel "3'"


## -----------------
## -- Plot profile

# ## D0 (mESCs)
# plotProfile \
# -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D0_ud500bp.gz \
# -o /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Plots/WT_D0_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D0_ud500bp_Profile.pdf \
# --plotType lines --plotHeight 6 --plotWidth 6 \
# --colors "#E5614E" "#F0AC4B" "#77BEF0" \
# --yAxisLabel "Enhancers" \
# --yMin 0 0 0 --yMax 1.8 0.7 2.0 \
# --regionsLabel "actv" "prim" "inac" \
# --samplesLabel "WT-D0_K4me1" "WT-D0_MATAC" "WT-D0_K27ac" \
# --startLabel "5'" --endLabel "3'"
# 
# ## D2 (EpiLCs)
# plotProfile \
# -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D2_ud500bp.gz \
# -o /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Plots/WT_D2_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D2_ud500bp_Profile.pdf \
# --plotType lines --plotHeight 6 --plotWidth 6 \
# --colors "#E5614E" "#F0AC4B" "#77BEF0" \
# --yAxisLabel "Enhancers" \
# --yMin 0 0 0 --yMax 1.8 0.7 2.0 \
# --regionsLabel "actv" "prim" "inac" \
# --samplesLabel "WT-D2_K4me1" "WT-D2_MATAC" "WT-D2_K27ac" \
# --startLabel "5'" --endLabel "3'"
# 
# ## D4 (ExEpiLCs)
# plotProfile \
# -m /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Matrix/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D4_ud500bp.gz \
# -o /scratch/Marlet/computeMatrix/Enhancers_intersect_v6/Plots/WT_D4_K4me1-MATAC-K27ac_WT_Master_Enhancer_Peakome_v6_D4_ud500bp_Profile.pdf \
# --plotType lines --plotHeight 6 --plotWidth 6 \
# --colors "#E5614E" "#F0AC4B" "#77BEF0" \
# --yAxisLabel "Enhancers" \
# --yMin 0 0 0 --yMax 1.8 0.7 2.0 \
# --regionsLabel "actv" "prim" "inac" \
# --samplesLabel "WT-D4_K4me1" "WT-D4_MATAC" "WT-D4_K27ac" \
# --startLabel "5'" --endLabel "3'"
