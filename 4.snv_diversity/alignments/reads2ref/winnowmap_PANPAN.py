import numpy as np
import pandas as pd
import os 

#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/human_GRCh38.p14/hg38_HGSVC/hg38.no_alt.fa.gz" # human no ALT
#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fa"
#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta"
#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta"

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
print(sample_table)
samples = sample_table['Specimen'].unique()
print(samples)
sample_ids = sample_table['SampleID'].unique()
print(sample_ids)

def get_input_merged(wildcards):
    return expand(
      'bam_mapped/{sample_run}.flagfilt.sorted.rg.bam',
      sample_run = sample_table[sample_table['Specimen']==wildcards.sample]['SampleID'].values
      )

rule all:
    input:
        #expand('mPonAbe1/bam_final/{sample}.sorted.bam', sample=samples),
        expand('mPanTro3/bam_final_PBmixRevio1474_1_D01/{sample}.sorted.bam', sample=['AG18352_2']),

#rule minimap2_fastqs_ref:
#    input: '../PanPan_HiFiAdapterFilter/{sample_run}.ccs.filt.fastq.gz'
#    output: temp('bam_mapped/{sample_run}.flagfilt.sorted.rg.bam')
#    shell: """
#    pbmm2 align {ref_path} {input} {output} --sort --preset CCS --rg '@RG\tID:{wildcards.sample_run}\tSM:{wildcards.sample_run}'    
#    """

rule meryl_kmer:
    output: 
        directory('mPanTro3/merylDB'),
        'mPanTro3/repetitive_k15_mPanTro3.txt'
    shell: """
    meryl count k=15 output {output[0]} {ref_path} &&
    meryl print greater-than distinct=0.9998 {output[0]} > {output[1]} 
    """

rule winnowmap2_fastqs_ref:
    input: 
        'mPanTro3/repetitive_k15_mPanTro3.txt',
        '../PanPan_HiFiAdapterFilter/{sample_run}.ccs.filt.fastq.gz',
    output: 
        temp('bam_mapped/{sample_run}.flagfilt.sorted.rg.bam')
    shell: """ 
    winnowmap -W {input[0]} -x map-pb -a -Y -L --eqx  --cs {ref_path} {input[1]} -R '@RG\tID:{wildcards.sample_run}'  | samtools view -q 10  -hbT {ref_path} -o {output}
    """

rule merge_to_samples:
    input:  get_input_merged 
    output: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.bam'
    shell: 'samtools merge -f {output} {input}'

# -SM tag 
rule add_readgroup_SMtag:
    input: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.bam'
    output: temp('bam_merge/{sample}.flagfilt.sorted.rg.merged.SM.bam')
    shell: 'samtools addreplacerg -r SM:{wildcards.sample} -r ID:{wildcards.sample} -r PL:PACBIO -o {output} {input}'

rule sort_merged:
    input: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.SM.bam'
    output: 'mPanTro3/bam_final_PBmixRevio1474_1_D01/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """




    
