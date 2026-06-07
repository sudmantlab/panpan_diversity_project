#minimap2 -cx asm5 --cs ../ncbi_dataset/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna  ../mPanTro3.pri.cur.20231031.fasta -t 24  > ncbi_october.paf &&  

### THANKSGIVING HALLOWEEN --> wrong order
#minimap2 -cx asm5 --cs -t 2 mPanTro3.pri.cur.20231122.fasta mPanTro3.pri.cur.20231031.fasta > thanksgiving_halloween.paf
#gsort -k6,6V -k8,8n thanksgiving_halloween.paf > thanksgiving_halloween.srt.paf
#k8 ./paftools.js call thanksgiving_halloween.srt.paf > thanksgiving_halloween.var.txt
#k8 ./paftools.js call -f  mPanTro3.pri.cur.20231122.fasta  -s thanksgiving_halloween  thanksgiving_halloween.srt.paf > thanksgiving_halloween.vcf
#k8 ./paftools.js stat thanksgiving_halloween.srt.paf
#paftools.js call ncbi_october.srt.paf  > ncbi_october.var.txt

### NCBI HALLOWEEN
#minimap2 -cx asm5 --cs -t 4 ncbi.fna mPanTro3.pri.cur.20231031.fasta > ncbi_halloween.paf
#gsort -k6,6V -k8,8n ncbi_halloween.paf > ncbi_halloween.srt.paf
#k8 ./paftools.js cal lncbi_halloween.srt.paf > ncbi_halloween.var.txt
#k8 ./paftools.js call -f  ncbi.fna  -s ncbi_halloween  ncbi_halloween.srt.paf > ncbi_halloween.vcf

### NCBI THANKSGIVING
#minimap2 -cx asm5 --cs -t 4 ncbi.fna mPanTro3.pri.cur.20231122.fasta > ncbi_thanksgiving.paf
#gsort -k6,6V -k8,8n ncbi_thanksgiving .paf > ncbi_thanksgiving.srt.paf
#k8 ./paftools.js call ncbi_thanksgiving.srt.paf > ncbi_thanksgiving.var.txt
#k8 ./paftools.js call -f  ncbi.fna  -s ncbi_thanksgiving ncbi_thanksgiving.srt.paf > ncbi_thanksgiving.vcf

### LIFTOFF ANNOTATIONS

#liftoff -g ncbi.renamed.gff.gz -o mPanTro3.pri.cur.20231031_liftoff_from_ncbi.gff.gz -chroms chroms_map.txt mPanTro3.pri.cur.20231031.fasta.gz ncbi.fna.gz 
#liftoff -g ncbi.renamed.gff.gz -o mPanTro3.pri.cur.20231122_liftoff_from_ncbi.gff.gz -chroms chroms_map.txt -p 1 mPanTro3.pri.cur.20231122.fasta.gz ncbi.fna.gz


### HALLOWEEN THANKSGIVING - Friday July 2025
minimap2 -cx asm5 --cs -t 8 mPanTro3.pri.cur.20231031.fasta mPanTro3.pri.cur.20231122.fasta > pafs_ref_debug/thanksgiving_halloween.paf
sort -k6,6V -k8,8n pafs_ref_debug/thanksgiving_halloween.paf > pafs_ref_debug/thanksgiving_halloween.srt.paf
./paf2chain --input thanksgiving_halloween.paf > thanksgiving_halloween.chain


### Bcftools liftover
export BCFTOOLS_PLUGINS=/global/scratch/users/joana_rocha/software/bcftools-1.20/plugins/
~/local/bin/bcftools +liftover pantros_mapped2mPanTro3_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz -Oz -o lifted_pantros_mapped2mPanTro3_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz -- \
    -s /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta \
    -f /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231122.fasta \
    -c /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/pafs_ref_debug/thanksgiving_halloween.chain

Lines   total/swapped/reference added/rejected: 47107694/0/0/0

export BCFTOOLS_PLUGINS=/global/scratch/users/joana_rocha/software/bcftools-1.20/plugins/
~/local/bin/bcftools +liftover pantros_mapped2mPanTro3_wholegenome.ALLSITES.sorted.filtered.vcf.gz -Oz -o lifted_pantros_mapped2mPanTro3_wholegenome.ALLSITES.sorted.filtered.vcf.gz -- \
   -s /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta \
   -f/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231122.fasta \
   -c /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/pafs_ref_debug/thanksgiving_halloween.chain

export BCFTOOLS_PLUGINS=/global/scratch/users/joana_rocha/software/bcftools-1.20/plugins/
~/local/bin/bcftools +liftover pantros_mapped2mPanTro3.sorted.filtered.vcf.gz  -Oz -o new_lifted_pantros_mapped2mPanTro3.sorted.filtered.vcf.gz  -- \
    -s /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta \
    -f/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231122.fasta \
    -c /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/pafs_ref_debug/thanksgiving_halloween.chain

Lines   total/swapped/reference added/rejected: 29638083/0/0/0
#your_input.vcf.gz: Your original VCF file, indexed with tabix.
#lifted_over.vcf.gz: The output file with coordinates in the new reference system.
#old_reference.fa: The source FASTA (the one your VCF is based on).
#new_reference.fa: The target FASTA (the one you are lifting over to).
#my_alignment.chain: The chain file you just created with paf2chain.

vcfs of all SVs
/global/scratch/users/nicolas931010/sv_detection/concat_vcf/ht2t/hprc.all_svs.vcf.gz
/global/scratch/users/nicolas931010/sv_detection/concat_vcf/mPanPan1/panpan-pp.all_svs.vcf.gz
/global/scratch/users/nicolas931010/sv_detection/concat_vcf/mPanTro3/panpan-pt.all_svs.vcf.gz

vcfs of SVs <100kb
/global/scratch/users/nicolas931010/sv_detection/truvari/ht2t/hprc.concat.vcf.gz
/global/scratch/users/nicolas931010/sv_detection/truvari/mPanPan1/panpan-pp.concat.vcf.gz
/global/scratch/users/nicolas931010/sv_detection/truvari/mPanTro3/panpan-pt.concat.vcf.gz

#### run vep SNVS from long reads made diploid phased (vcfs from PANPAN_Singer)

## human HPRC (long)
mkdir -p vep_snvs/ht2t
grep -v "#" /global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/ncbi_dataset/data/GCF_009914755.1/genomic.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > vep_snvs/ht2t/ht2t.gff.gz
tabix -p gff -C vep_snvs/ht2t/ht2t.gff.gz
vep -i hprc_mapped2ht2t.sorted.filtered.vcf.gz  \
    -o vep_snvs/ht2t/hprc_mapped2ht2t.phased_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/nicolas931010/sv_detection/reference/ht2t.fasta \
    --force_overwrite \
    --phased \
    --custom file=vep_snvs/ht2t/ht2t.gff.gz,short_name=ht2t,format=gff,type=overlap \
    &> vep_snvs/hprc_mapped2ht2t.phased_snvs.vep.log
   
## human HPRC + HGSVC3 (long)      
vep -i hprc_hgsvc_mapped2ht2t.sorted.filtered.vcf.gz  \
    -o vep_snvs/ht2t/hprc_hgsvc_mapped2ht2t.phased_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/nicolas931010/sv_detection/reference/ht2t.fasta \
    --force_overwrite \
    --phased \
    --custom file=vep_snvs/ht2t/ht2t.gff.gz,short_name=ht2t,format=gff,type=overlap \
    &> vep_snvs/hprc_hgsvc_mapped2ht2t.phased_snvs.vep.log  
   
### mPanTros
## long
mkdir -p vep_snvs/mPanTro3
grep -v "#" /global/scratch/users/joana_rocha/PANPAN/PANPAN_graph/pangene/mPanTro3_reference_proteins/mPanTro3_modified_chrom_alias.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > vep_snvs/mPanTro3/mPanTro3.gtf.gz
tabix -p gff -C vep_snvs/mPanTro3/mPanTro3.gtf.gz
vep -i lifted_pantros_mapped2mPanTro3.sorted.filtered.vcf.gz  \
    -o vep_snvs/mPanTro3/lifted_pantros_mapped2mPanTro3.phased_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231122.fasta  \
    --force_overwrite \
    --phased \
    --custom file=vep_snvs/mPanTro3/mPanTro3.gtf.gz,short_name=mPanTro3,format=gtf,type=overlap \
    &> vep_snvs/lifted_pantros_mapped2mPanTro3.phased_snvs.vep.log


#mkdir snpEff/data/mPanTro3
#copy gff and fa annotation and name it genes.gtf and sequence.fa
snpEff -Xmx8G build -gtf22 -v mPanTro3
snpEff -Xmx8G build -gtf22 -v mPanTro3 -noCheckCds -noCheckProtein
snpEff databases | grep mPanTro3
snpEff -v mPanTro3 -stats mPanTro3_report.html ../lifted_pantros_mapped2mPanTro3.sorted.filtered.vcf.gz > lifted_pantros_mapped2mPanTro3_snpEff_annotated.vcf

# short + long
mkdir -p vep_snvs/mPanTro3
ln -s /global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/vep_snvs/mPanTro3/mPanTro3.gtf.gz vep_snvs/mPanTro3/.
ln -s /global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/vep_snvs/mPanTro3/mPanTro3.gtf.gz.csi vep_snvs/mPanTro3/.
vep -i  lifted_pantros_mapped2mPanTro3_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz  \
    -o vep_snvs/mPanTro3/lifted_pantros_mapped2mPanTro3.shortlong_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231122.fasta  \
    --force_overwrite \
    --custom file=vep_snvs/mPanTro3/mPanTro3.gtf.gz,short_name=mPanTro3,format=gtf,type=overlap \
    &> vep_snvs/lifted_pantros_mapped2mPanTro3.shortlong_snvs.vep.log


### mPanPan1
#long
mkdir -p vep_snvs/mPanPan1
grep -v "#" /global/scratch/users/joana_rocha/PANPAN/PANPAN_graph/pangene/mPanPan1_reference_proteins/mPanPan1_modified_chrom_alias.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > vep_snvs/mPanPan1/mPanPan1.gtf.gz
tabix -p gff -C vep_snvs/mPanPan1/mPanPan1.gtf.gz
vep -i panpa_mapped2mPanPan1.sorted.filtered.vcf.gz  \
    -o vep_snvs/mPanPan1/panpa_mapped2mPanTro3.phased_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/nicolas931010/sv_detection/reference/mPanPan1.fasta  \
    --force_overwrite \
    --phased \
    --custom file=vep_snvs/mPanPan1/mPanPan1.gtf.gz,short_name=mPanPan1,format=gtf,type=overlap \
    &> vep_snvs/panpa_mapped2mPanPan1.phased_snvs.vep.log

# short + long
mkdir -p vep_snvs/mPanPan1
ln -s /global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/vep_snvs/mPanPan1/mPanPan1.gtf.gz vep_snvs/mPanPan1/.
ln -s /global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/vep_snvs/mPanPan1/mPanPan1.gtf.gz.csi vep_snvs/mPanPan1/.
vep -i /global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/BONOBOS/vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz  \
    -o vep_snvs/mPanPan1/panpa_mapped2mPanPan1.shortlong_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/nicolas931010/sv_detection/reference/mPanPan1.fasta  \
    --force_overwrite \
    --custom file=vep_snvs/mPanPan1/mPanPan1.gtf.gz,short_name=mPanPan1,format=gtf,type=overlap \
    &> vep_snvs/panpa_mapped2mPanPan1.shortlong_snvs.vep.log


#vep -i unphased_snvs.vcf.gz \
#    -o vep_output/unphased_snvs.vep.txt \
#    --fork 32 \
#    --fasta reference/ht2t.fasta \
#    --force_overwrite \
#    --custom file=vep/ht2t/ht2t.gff.gz,short_name=ht2t,format=gff,type=overlap \
#    &> vep_output/unphased_snvs.vep.log



#### run Haplossaurus SNVS from long reads made diploid phased (vcfs from PANPAN_Singer)