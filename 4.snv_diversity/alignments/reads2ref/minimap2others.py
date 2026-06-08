import numpy as np
import pandas as pd
import os 

ref_path = "/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fna"
#ref_path = "/global/scratch/users/joana_rocha/PANPAN/output/hifiasm-fasta/PR01227/joana_settings_shortcut/no_opts/PR01227.p_ctg.fa" 
sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', '\t')
print(sample_table)
samples = sample_table['Specimen'].unique()
print(samples)

def get_input_merged(wildcards):
    return expand(
      'bam_mapped/{sample_run}.flagfilt.sorted.rg.bam',
      sample_run = sample_table[sample_table['Specimen']==wildcards.sample]['SampleID'].values
      )

rule all:
    input:
        expand('bam_final_Clint/{sample}.sorted.bam', sample=samples)

rule minimap2_fastqs_ref:
    input: '{sample_run}.ccs.filt.fastq.gz'
    output: temp('bam_mapped/{sample_run}.flagfilt.bam')
    shell: './minimap2 -ax map-hifi -t 20 {ref_path} {input} | samtools view -q 10  -bT {ref_path} -o {output}'


rule minimap2_sort:
    input: 'bam_mapped/{sample_run}.flagfilt.bam'
    output: temp('bam_mapped/{sample_run}.flagfilt.sorted.bam')
    shell: 'samtools sort -o {output} {input}'

#no need to remove duplicates because there are no optical duplicates
rule add_readgroup:
    input: 'bam_mapped/{sample_run}.flagfilt.sorted.bam'
    output: temp('bam_mapped/{sample_run}.flagfilt.sorted.rg.bam')
    shell: 'samtools addreplacerg -r ID:{wildcards.sample_run} -r PL:PACBIO -o {output} {input}'

rule merge_to_samples:
    input:  get_input_merged 
    output: temp('bam_merge/{sample}.flagfilt.sorted.rg.merged.bam')
    shell: 'samtools merge -f {output} {input}'

rule sort_merged:
    input: 'bam_merge/{sample}.flagfilt.sorted.rg.merged.bam'
    output: 'bam_final_Clint/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """




    
