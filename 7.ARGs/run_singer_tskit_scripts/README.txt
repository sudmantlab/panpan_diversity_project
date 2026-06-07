
#'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
#'mPanTro3' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta",
#'mPonAbe1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta",
#'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"

#downloaded the chromosome sizes file
wget https://github.com/marbl/T2T-Browser/raw/main/T2Tgenomes/mPanTro3_v2.0_pri/mPanTro3_v2.0_pri.sizes
wget https://github.com/marbl/T2T-Browser/raw/main/T2Tgenomes/mPanPan1_v2.0_pri/mPanPan1_v2.0_pri.sizes
wget https://github.com/marbl/T2T-Browser/raw/main/T2Tgenomes/mPonAbe1_v2.0_pri/mPonAbe1_v2.0_pri.sizes
wget https://github.com/marbl/T2T-Browser/blob/main/T2Tgenomes/T2T-CHM13v2.0/chm13v2.0.sizes

wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/variants/mPanPan1_v2.0.accessibility_mask.bb

# Converted to BED format (start=0, end=chromosome_length)
awk '{print $1 "\t0\t" $2}' mPanTro3_v2.0_pri.sizes > chromosomes_full.bed
awk '{print $1 "\t0\t" $2}' mPanPan1_v2.0_pri.sizes > chromosomes_full.bed
awk '{print $1 "\t0\t" $2}' chm13v2.0.sizes > chromosomes_full.bed
awk '{print $1 "\t0\t" $2}' mPonAbe1_v2.0_pri.sizes > chromosomes_full.bed

#downloaded the telomeres and centromeres
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/pattern/mPanTro3_v2.0.GenomeFeature_v0.9.bb mPanTro3_genomefeatures.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/pattern/mPanPan1_v2.0.GenomeFeature_v1.0.bb mPanPan1_genomefeatures.bed
bigBedToBed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/CHM13/bbi/genomeFeature_v1.1.bb ht2t_genomefeatures.bed ht2t_genomefeatures.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/pattern/mPonAbe1_v2.0.GenomeFeature_v1.0.bb mPonAbe1_genomefeatures.bed

#filter to get coordinates for primary hap (hap1)
cat mPanTro3_genomefeatures.bed | grep 'hap1' > mPanTro3_hap1_genomefeatures.bed
awk '$4 == "Cen"' mPanTro3_hap1_genomefeatures.bed > centromeres.bed
awk '$4 == "Cen"' mPanPan1_genomefeatures.bed > centromeres.bed
awk '$4 == "CEN"' ht2t_genomefeatures.bed   > centromeres.bed
awk '$4 == "Cen"'  mPonAbe1_genomefeatures.bed   > centromeres.bed


# Sort chromosome regions and centromere annotations
sort -k1,1 -k2,2n chromosomes_full.bed > chromosomes_sorted.bed
sort -k1,1 -k2,2n centromeres.bed > centromeres_sorted.bed 

# manually removed MT from chromosomes_sorted.bed 
# manually added hap2 Y to centromeres_sorted.bed in mPanTro3
# manually removed unitigs chromosomes_sorted.bed from bonobo

# had to edit centromeres bed in bonobos to just keep the pat or mat cen that correspond to the primary assembly

# Split chromosomes into regions not overlapping centromeres
bedtools subtract \
    -a chromosomes_sorted.bed \
    -b centromeres_sorted.bed \
    > chromosomes_without_centromeres.bed

# Check output (example for chr1_hap1_hsa1):
grep "chr1_hap1_hsa1" chromosomes_without_centromeres.bed

# Expected output (if centromere was at 101,356,245-105,077,980):
# chr1_hap1_hsa1  0               101356245
# chr1_hap1_hsa1  105077980       231299448

bedtools makewindows -b chromosomes_without_centromeres.bed -w 5000000 > 5mb_windows.bed

tabix -p vcf pantros_mapped2mPanTro3.sorted.filtered.vcf.gz 
tabix -p vcf panpa_mapped2mPanPan1.sorted.filtered.vcf.gz 
tabix -p vcf hprc_hgsvc_mapped2ht2t.sorted.filtered.vcf.gz
tabix -p vcf allsamples_mapped2mPonAbe1.sorted.filtered.vcf.gz

# Split the VCF into chunks
mkdir -p vcf_chunks
while IFS=$'\t' read -r chr start end; do
  bcftools view \
    -r "${chr}:${start}-${end}" \
    -o "vcf_chunks/${chr}_${start}_${end}.vcf" \
    pantros_mapped2mPanTro3.sorted.filtered.vcf.gz
done < 5mb_windows.bed


while IFS=$'\t' read -r chr start end; do
  bcftools view \
    -r "${chr}:${start}-${end}" \
    -o "vcf_chunks/${chr}_${start}_${end}.vcf" \
    panpa_mapped2mPanPan1.sorted.filtered.vcf.gz
done < 5mb_windows.bed


##### HUMANS ##########################

python convert_vcf_chroms.py chrom_alias_ht2t.tsv hprc_hgsvc_mapped2ht2t.sorted.filtered.vcf.gz   hprc_hgsvc_mapped2ht2t.sorted.filtered_renamed.vcf.gz   
tabix -p vcf hprc_hgsvc_mapped2ht2t.sorted.filtered_renamed.vcf.gz 

while IFS=$'\t' read -r chr start end; do
  bcftools view \
    -r "${chr}:${start}-${end}" \
    -o "vcf_chunks/${chr}_${start}_${end}.vcf" \
    hprc_hgsvc_mapped2ht2t.sorted.filtered_renamed.vcf.gz
done < 5mb_windows.bed


python convert_vcf_chroms.py chrom_alias_ht2t.tsv human579_mapped2ht2t.sorted.filtered.vcf.gz    human579_mapped2ht2t.sorted.filtered_renamed.vcf.gz   
tabix -p vcf human579_mapped2ht2t.sorted.filtered_renamed.vcf.gz
mkdir vcf_chunks/
while IFS=$'\t' read -r chr start end; do
  bcftools view \
    -r "${chr}:${start}-${end}" \
    -o "vcf_chunks/${chr}_${start}_${end}.vcf" \
    human579_mapped2ht2t.sorted.filtered_renamed.vcf.gz
done < 5mb_windows.bed





#### ALL ##############
while IFS=$'\t' read -r chr start end; do
  bcftools view \
    -r "${chr}:${start}-${end}" \
    -o "vcf_chunks/${chr}_${start}_${end}.vcf" \
    allsamples_mapped2mPonAbe1.sorted.filtered.vcf.gz
done < 5mb_windows.bed


# snakemake -s singer_genomewide.py -n 24 --rerun-incomplete -n 

#download CenSats + SegDups + Short Read Accessibility + RepeatMasker - https://github.com/hloucks/CenSatData/tree/main/T2TPrimates

bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_CenSat_v2.0.bb mPanTro3_CenSat.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_CenSat_v2.0.bb mPanPan1_CenSat.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_CenSat_v2.0.bb mPonAbe1_CenSat.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0.SD_v1.0.bb mPanTro3_SegDups.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0.SD_v1.0.bb mPanPan1_SegDups.bed
bigBedToBed https://hgdownload.soe.ucsc.edu/gbdb/hs1/sedefSegDups/sedefSegDups.bb ht2t_SegDups.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0.SD_v1.0.bb mPonAbe1_SegDups.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/variants/mPanTro3_v2.0.accessibility_mask.bb mPanTro3_SR_mask.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0.RepeatMasker_v1.2.bb  mPanTro3_RepeatMasker.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/variants/mPanPan1_v2.0.accessibility_mask.bb mPanPan1_SR_mask.bed
bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0.RepeatMasker_v1.2.1.bb mPanPan1_RepeatMasker.bed
bigBedToBed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/hs1.combined_mask.bb ht2t_SR_mask.bed
bigBedToBed https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker/chm13v2.0_rmsk.bb ht2t_RepeatMasker.bed


############## BONOBOS #########################

### Subtract censats and funky stuff for mPanPan1 results
bedtools intersect -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results/genome-wide_metrics_annotated.csv) \
                   -b mPanPan1_CenSat.bed \
                   -v | \
cut -f4- | \
sed '1i chromosome,start,end,avg_tmrca,avg_pairwise_coalescence_time,genes' > results/genome-wide_metrics_annotated_noCensat.csv

### Keep only windows that overlap with the SR Accessibility mask
bedtools intersect \
    -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results/genome-wide_metrics_annotated_noCensat.csv) \
    -b <(awk '{print $1"\t"$2"\t"$3}' mPanPan1_SR_mask.bed) \
    -u \
    | cut -f4- \
    | sed '1i chromosome,start,end,avg_tmrca,avg_pairwise_coalescence_time,genes' > results/genome-wide_metrics_annotatedbyPop_SR_A_MASK_only.csv

############## CHIMPANZEEE #########################

#### SINGER OUTPUT ####

### Subtract censats for mPanTro3 results
bedtools intersect -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results/genome-wide_metrics_annotatedbyPop.csv) \
                   -b mPanTro3_CenSat.bed \
                   -v | \
cut -f4- | \
sed '1i chromosome,start,end,population,avg_tmrca,T_pooled,T_within,Tpooled_Twithin_ratio,genes' > results/genome-wide_metrics_annotatedbyPop_CenSat_removed.csv

### Keep only windows that overlap with the SR Accessibility mask
bedtools intersect -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results/genome-wide_metrics_annotatedbyPop_CenSat_removed.csv) \
                 -b mPanTro3_SR_mask.bed \
                 -u \
                 | cut -f4- \
                 | sed '1i chromosome,start,end,population,avg_tmrca,T_pooled,T_within,Tpooled_Twithin_ratio,genes' > results/genome-wide_metrics_annotatedbyPop_SR_A_MASK_only.csv



############## HUMANS #########################

#### SINGER OUTPUT ####

### censats are in genome features and funky stuff for humans results
### Keep only windows that overlap with the SR Accessibility mask
bedtools intersect -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results/genome-wide_metrics_annotatedbyPop.csv) \
                 -b ht2t_SR_mask.bed \
                 -u \
                 | cut -f4- \
                 | sed '1i chromosome,start,end,population,avg_tmrca,T_pooled,T_within,Tpooled_Twithin_ratio,genes' > results/genome-wide_metrics_annotatedbyPop_SR_A_MASK_only.csv




#### PIXY OUTPUT ####

bedtools intersect -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results/annotated_pi_ratios_final.csv) \
                   -b mPanTro3_CenSat.bed \
                   -v | \
cut -f4- | \
sed '1i chromosome,window_pos_1,window_pos_2,pop,avg_pi,no_sites,count_diffs,count_comparisons,count_missing,pi_pooled,pi_ratio,name' > results/annotated_pi_ratios_final_CenSat_removed.csv

### Keep only windows that overlap with the SR Accessibility mask
bedtools intersect -a <(awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$0}' results_per_pop/annotated_pi_ratios_final_CenSat_removed.csv) \
                 -b mPanTro3_SR_mask.bed \
                 -u \
                 | cut -f4- \
                 | sed '1i chromosome,window_pos_1,window_pos_2,pop,avg_pi,no_sites,count_diffs,count_comparisons,count_missing,pi_pooled,pi_ratio,name' > results_per_pop/annotated_pi_ratios_SR_A_MASK_only.csv



