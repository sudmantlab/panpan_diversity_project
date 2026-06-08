import numpy as np
import pandas as pd

scaffold_ids= pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.bed', '\t', header=None, usecols=[0])[0].values


ref_fasta = '/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fna'

rule all:
    input:
        expand('sites/allchimps_{scaffold_id}_qc_pass.rf', scaffold_id=scaffold_ids),
        expand('sites/allchimps_{scaffold_id}_qc_fail.bz2', scaffold_id=scaffold_ids),
        expand('sites/allchimps_{scaffold_id}_qc_pass.vcf.bz2', scaffold_id=scaffold_ids),


rule run_snp_cleaner:
    input:
        'allchimps_bams_withStacy.txt'
    output:
        'sites/allchimps_{scaffold_id}_qc_pass.rf',
        'sites/allchimps_{scaffold_id}_qc_fail.bz2',
        'sites/allchimps_{scaffold_id}_qc_pass.vcf.bz2'
    shell: """
    bcftools mpileup --threads 20 -f {ref_fasta} -b {input} -r {wildcards.scaffold_id} -a SP,DP -q 30 -Q 20 | bcftools call --skip-variants indels -f GQ -c -| /global/home/users/joana_rocha/.local/ngsQC/snpCleaner/snpCleaner.pl -k 12 -u 1 -Q 30 -h 0 -H 1e-6 -b 1e-20 -S 1e-4 -f 1e-4 -e 1e-4 -v -B {output[0]} -p {output[1]} | bzip2 -c > {output[2]}
    """

# 40% of the individuals covered 1 
