## Supp Figure 2: DNA methylation heatmaps

## See 2025_07_25 and 2025_11_12
cd /scratch/Marlet/

## --------- ##
## -- DKO -- ##
## --------- ##


## -- Calculate DNA methylation matrix

sbatch -n 12 computeMatrix scale-regions \
-R \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_1_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_2_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_3_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_4_ordered.bed3 \
/scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/bed/Cluster_5_ordered.bed3 \
-S \
/data/priscillia_data/bigWigs/mergedreps/MATAC/DKO_D0_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/DKO_D2_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/mergedreps/MATAC/DKO_D4_MATAC_mergedreps_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_D0_Tet12DKO_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/priscillia_data/bigWigs/E14_D2_Tet12DKO_WGBS_CROP_trimV1_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_D4_Tet12DKO_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
--upstream 500 \
--downstream 500 \
--regionBodyLength 500 \
--binSize 50 \
--sortRegions keep \
--blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
--outFileName /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_v3.gz \
--numberOfProcessors 12
# Submitted batch job 34058
# 		The following chromosome names did not match between the bigwig files
# 		chromosome      length
# 		chr5_JH584296_random        199368
# 		chr4_JH584292_random         14945
# 		 chrUn_GL456367      42057
# 		 chrUn_GL456387      24685
# 		chr5_JH584297_random        205776
# 		 chrUn_GL456394      24323
# 		chr1_GL456210_random        169725
# 		 chrUn_GL456381      25871
# 		chrY_JH584303_random        158099
# 		chr5_JH584299_random        953012
# 		 chrUn_GL456379      72385
# 		 chrUn_GL456239      40056
# 		chr5_JH584298_random        184189
# 		 chrUn_GL456360      31704
# 		chr4_JH584293_random        207968
# 		chr5_GL456354_random        195993
# 		chr7_GL456219_random        175968
# 		 chrUn_GL456372      28664
# 		 chrUn_GL456382      23158
# 		chr4_GL456350_random        227966
# 		chrY_JH584301_random        259875
# 		chr4_JH584295_random          1976


## -- Plot heatmap
plotHeatmap \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_v3.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_2/Supp_DNAme/DKO_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_Heatmap_v3.pdf \
--sortRegions keep --interpolationMethod nearest \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 28 \
--heatmapWidth 4 \
--zMin 0 --zMax 100 --colorMap RdPu --missingDataColor silver \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "DKO-D0_ATAC" "DKO-D2_ATAC" "DKO-D4_ATAC" "DKO-D0_WGBS" "DKO-D2_WGBS" "DKO-D4_WGBS" \
--xAxisLabel "NonProm DA peak" --startLabel "5'" --endLabel "3'"

## -- Plot profile
plotProfile \
-m /scratch/Marlet/TCseq/WT_MATAC/Dynamic/K5/computeMatrix/DKO_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_v3.gz \
-o /scratch/Marlet/Manuscript_Plots/Figure_2/Supp_DNAme/DKO_D0-D2-D4_MATAC-WGBS_5mC_K5_sorted_Profile_v3.pdf \
--plotHeight 6 \
--plotWidth 5 \
--plotType lines \
--colors "#E4E319" "#90D643" "#2D908C" "#5A4791" "#571B67" \
--yMin 0 --yMax 100 --yAxisLabel "5mC" \
--regionsLabel "c1" "c2" "c3" "c4" "c5" \
--samplesLabel "DKO-D0_ATAC" "DKO-D2_ATAC" "DKO-D4_ATAC" "DKO-D0_WGBS" "DKO-D2_WGBS" "DKO-D4_WGBS" \
--startLabel "5'" --endLabel "3'"


