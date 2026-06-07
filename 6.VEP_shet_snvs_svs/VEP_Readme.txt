#### run vep SNVS from long reads made diploid phased

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
vep -i pantros_mapped2mPanTro3.sorted.filtered.vcf.gz  \
    -o vep_snvs/mPanTro3/pantros_mapped2mPanTro3.phased_snvs.vep.txt \
    --fork 32 \
    --fasta /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231122.fasta  \
    --force_overwrite \
    --phased \
    --custom file=vep_snvs/mPanTro3/mPanTro3.gtf.gz,short_name=mPanTro3,format=gtf,type=overlap \
    &> vep_snvs/lifted_pantros_mapped2mPanTro3.phased_snvs.vep.log


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
