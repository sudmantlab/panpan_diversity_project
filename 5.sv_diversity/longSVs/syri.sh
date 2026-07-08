#!/bin/bash

### setup
species="human"
ref="/XXX/YYY/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
fastaDir="/XXX/YYY/fastas"
pafDir="/XXX/YYY/pafs"
slurmScript="/XXX/YYY/syri.slurm"

### check reference exists
if [[ ! -f $ref ]]
then
    echo "missing ref: $ref"
    exit 1
fi


mkdir -p $species
cd $species

for paf in $pafDir/*.paf
do 
    label=$(basename $paf .paf)
    qry="${fastaDir}/${label}.fasta"

    ### check query exists
    if [[ ! -f $qry ]]
    then
        echo "missing qry: $qry"
        exit 1
    fi

    export ref
    export qry
    export paf
    export label
    sbatch --export ref,qry,paf,label $slurmScript
done
