import numpy as np
import pandas as pd
import os 

ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta.gz", ### no difference to thanksgiving/final primary
#ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta"
sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC_PacBio_pepsamples.tsv', sep='\t') #has the correspondence between HPRC sample run IDs with the HPRC sample names
#sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC_PacBio_pepsamples_HG005.tsv', sep='\t') #has the correspondence between HPRC sample run IDs with the HPRC sample names


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
        #expand("../HPRC_raw/{sample_id}.ccs.filt.fastq.gz", sample_id=sample_ids),
        #expand("HPRC_raw/{sample_id}.filt.fastq.gz", sample_id=HG005_ids),
        #expand('ht2t/bam_final/{sample}.sorted.bam', sample=samples),
        #expand('mPonAbe1/bam_final/{sample}.sorted.bam', sample=['NA20129']), ### wait until mPonAbe after Dec1
        expand('ht2t/bam_final/{sample}.sorted.bam', sample=['NA20129']), ### wait until mPonAbe after Dec1
        #expand('mPanTro3/bam_final/{sample}.sorted.bam', sample=samples) 

rule HiFiAdapterFilt:
    input: "../HPRC_raw/{sample_id}.ccs.bam" 
    output: 
        "../HPRC_raw/{sample_id}.ccs.filt.fastq.gz",
    shell: """
        export PATH=$PATH:/global/scratch/users/joana_rocha/PANPAN/code/HiFiAdapterFilt
        export PATH=$PATH:/global/scratch/users/joana_rocha/PANPAN/code/HiFiAdapterFilt/DB
        cd ../HPRC_raw/ &&
        /global/scratch/users/joana_rocha/PANPAN/code/HiFiAdapterFilt/pbadapterfilt.sh -p {wildcards.sample_id} -t 20 
        cd -
    """
    
#rule HiFiAdapterFilt_fastq:
#    input: "HPRC_raw/{special_id}.fastq.gz" 
#    output: 
#        "HPRC_raw/{special_id}.filt.fastq.gz"
#    shell: """
#        export PATH=$PATH:/global/scratch/users/joana_rocha/PANPAN/code/HiFiAdapterFilt
#        export PATH=$PATH:/global/scratch/users/joana_rocha/PANPAN/code/HiFiAdapterFilt/DB
#        cd HPRC_raw/ &&
#        /global/scratch/users/joana_rocha/PANPAN/code/HiFiAdapterFilt/pbadapterfilt.sh -p {wildcards.special_id} -t 20 
#        cd -
#    """

#rule minimap2_fastqs_ref:
#    input: 'HPRC_raw/{sample_run}.ccs.filt.fastq.gz'
#    output: temp('bam_mapped/{sample_run}.flagfilt.bam')
#    shell: './minimap2 -ax map-hifi -t 20 {ref_path} {input} | samtools view -q 10  -bT {ref_path} -o {output}'

#rule minimap2_sort:
#    input: 'bam_mapped/{sample_run}.flagfilt.bam'
#    output: temp('bam_mapped/{sample_run}.flagfilt.sorted.bam')
#    shell: 'samtools sort -o {output} {input}'

#rule add_readgroup:
 #   input: 'bam_mapped/{sample_run}.flagfilt.sorted.bam'
  #  output: temp('bam_mapped/{sample_run}.flagfilt.sorted.rg.bam')
   # shell: 'samtools addreplacerg -r ID:{wildcards.sample_run} -r PL:PACBIO -o {output} {input}'

#rule pbmm2_fastqs_ref:
#    input: 'HPRC_raw/{sample_run}.ccs.filt.fastq.gz'
#    output: temp('bam_mapped/{sample_run}.flagfilt.sorted.rg.bam')
#    shell: """
#    pbmm2 align {ref_path} {input} {output} --sort --preset CCS --rg '@RG\tID:{wildcards.sample_run}'    
#    """

# making the k-mer table with meryl mery
rule meryl_kmer:
    output: 
        directory('merylDB'),
        #'mPonAbe1/repetitive_k15_mPonAbe1.txt'
        #'mPanTro3/repetitive_k15_mPanTro3.txt'
        'ht2t/repetitive_k15_CHM13.txt'
    shell: """
    meryl count k=15 output {output[0]} {ref_path} &&
    meryl print greater-than distinct=0.9998 {output[0]} > {output[1]} 
    """

rule winnowmap2_fastqs_ref_other:
    input: 
        #'mPonAbe1/repetitive_k15_mPonAbe1.txt',
        'ht2t/repetitive_k15_CHM13.txt',
        #'mPanTro3/repetitive_k15_mPanTro3.txt',
        '../HPRC_raw/{sample_run}.ccs.filt.fastq.gz',
        #'../HPRC_raw/{sample_run}.filt.fastq.gz',
    output: 
        'bam_mapped/{sample_run}.flagfilt.sorted.rg.bam'
    shell: """ 
    winnowmap -W {input[0]} -t 20 -x map-pb -a -Y -L --eqx  --cs {ref_path} {input[1]} -R '@RG\tID:{wildcards.sample_run}'  | samtools view -q 10  -hbT {ref_path} -o {output}
    """

rule merge_to_samples:
    input:  get_input_merged 
    output: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.bam'
    shell: 'samtools merge -f {output} {input}'

# -SM tag 
rule add_readgroup_SMtag:
    input: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.bam'
    output: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.SM.bam'
    shell: 'samtools addreplacerg -r SM:{wildcards.sample} -r ID:{wildcards.sample} -r PL:PACBIO -o {output} {input}'

rule sort_merged:
    input: 
        'bam_merge/{sample}.flagfilt.sorted.rg.merged.SM.bam'
    output:
        'ht2t/bam_final/{sample}.sorted.bam'
        #'mPanTro3/bam_final/{sample}.sorted.bam'  
	# 'mPonAbe1/bam_final/{sample}.sorted.bam' 
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """
