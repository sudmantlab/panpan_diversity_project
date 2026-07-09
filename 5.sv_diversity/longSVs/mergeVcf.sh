#!/bin/bash

set -euo pipefail

### setup
SPP="mPanPan1"
# SPP="ht2t"
# SPP="mPanTro3"

inputDir="/global/scratch/users/joana_rocha/PANPAN"
outputDir="/global/scratch/users/scott_ferguson/panpan"
scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "$SPP" in
    "mPanPan1")
        refFasta="${inputDir}/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
        syriDir="${outputDir}/syri/real/mPanPan1/out"
        qryDir="${inputDir}/PANPAN_pairwise_aln/test_svbyeye/ragtag_mPanPan1_all"
        ;;
    "ht2t")
        refFasta="${inputDir}/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
        syriDir="${outputDir}/syri/real/ht2t/out"
        qryDir="${inputDir}/PANPAN_pairwise_aln/test_svbyeye/ragtag_ht2t_all"
        ;;
    "mPanTro3")
        refFasta="${inputDir}/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta"
        syriDir="${outputDir}/syri/real/mPanTro3/out"
        qryDir="${inputDir}/PANPAN_pairwise_aln/test_svbyeye/ragtag_mPanTro3_all"
        ;;
esac


### process
mkdir -p ${SPP}/vcf
cd ${SPP}

for x in ${qryDir}/*.fasta
do
    [[ -e "$x" ]] || continue
    sample=$(basename $x .fasta)
    qryFasta="${qryDir}/${sample}.fasta"
    syriOut="${syriDir}/${sample}_${SPP}-syri.out"

    [[ -f "$syriOut" ]] || continue
    
    python3 "${SCRIPT_DIR}/syriToVcf.py" $refFasta $qryFasta $syriOut vcf/${sample}.vcf
done

for x in vcf/*.vcf
do
    [[ -e "$x" ]] || continue
    bcftools sort -Oz $x > ${x}.gz
    bcftools index -t ${x}.gz
    rm $x
done

cohort_vcfs=(vcf/*.vcf.gz)
bcftools merge -Oz -m none ${cohort_vcfs[@]} > ${SPP}.vcf.gz
bcftools index -t ${SPP}.vcf.gz


echo "bash ${SCRIPT_DIR}/truvari.sh ${SPP}.vcf.gz" > commands.txt

echo "INFO: Launching Truvari processing environment..."
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate truvari
set -u

if [[ -f commands.txt ]]; then
    parallel --jobs $(nproc) < commands.txt
fi
