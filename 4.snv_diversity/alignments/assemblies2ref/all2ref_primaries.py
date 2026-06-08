import numpy as np
import pandas as pd
import os
from snakemake.io import expand

# Paths to reference genomes
ref_paths = {
    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    'mPanTro3': "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta",
    'mPonAbe1': "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta",
    'mPanPan1': "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
}

# Load sample tables
#sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/chimpsamples.tsv', sep='\t')
sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/bonobosamples.tsv', sep='\t')
#sample_table_hprc = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_names.txt', names=['names'])
#sample_table_hic = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples_hic.tsv', sep='\t')

# Extract unique samples
references = ['ht2t', 'mPanTro3', 'mPanPan1', 'mPonAbe1']
samples = sample_table['Specimen'].unique()
#samples = ['mPanTro3']
#human_samples = sample_table_hprc['names'].unique()


# Print sample haplotypes for debugging
print(samples)
#print(human_samples)

rule all:
    input:
        # Add both rules to the all rule to ensure they are both run
        #directory(expand('ragtag_{ref}_all/{sample}.ragtag_output', sample=samples, ref=['mPanTro3'])),
        directory(expand('ragtag_{ref}_all/{anything}.ragtag_output', anything=samples, ref=['mPanPan1'])),
        #expand('ragtag_{ref}_all/{anything}.ragtag_output', anything=human_samples, ref=['ht2t']),
        #expand("ragtag_{ref}_all/{sample}.{hap}.chr.fasta",  sample=samples, hap=haps, ref=['mPonAbe1']),
        expand("ragtag_{ref}_all/{anything}.chr.fasta",  anything=samples, ref=['mPanPan1']),
        #expand("ragtag_{ref}_all/{samplehap}.chr.fasta",  samplehap=human_samples,  ref=['ht2t']),
        #expand("ragtag_{ref}_all/all2ref/{samplehap}_{ref}.paf", samplehap=sample_haps, ref=['mPonAbe1']),
        expand("ragtag_{ref}_all/all2ref/{sample}_{ref}.paf", sample=samples, ref=['mPanPan1']),
        #expand("ragtag_{ref}_all/all2ref/{samplehap}_{ref}.paf", samplehap=human_samples, ref=['ht2t']),
        
def get_fasta_path(wildcards):
    """
    Look for verkko, HiC, or regular assembly files in order of preference
    """
    verkko_path = f"/global/scratch/users/joana_rocha/PANPAN/Verkko-fasta_shortcut/{wildcards.sample}.verkko.fasta"
    hic_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta-HiC_shortcut/{wildcards.sample}.p_ctg.hic.fa"
    hifi_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{wildcards.sample}.p_ctg.fa"
    
    if os.path.exists(verkko_path):
        return verkko_path
    elif os.path.exists(hic_path):
        return hic_path
    else:
        return hifi_path

rule panpan2ref:
    input:
        fasta=get_fasta_path
        #fasta='mPanTro3.pri.cur.20231031.fasta'
    output: 
        directory('ragtag_{ref}_all/{sample}.ragtag_output')
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    ragtag.py scaffold -C -u --mm2-params '-x asm5' --unimap-params '-x asm5 -t 20' -o {output} {params.ref} {input.fasta} &&
    cd {output} &&
    for filename in ragtag.*; do
        newname="$(echo "$filename" | sed 's/ragtag/{wildcards.sample}/')"
        echo "Renaming $filename to $newname"
        mv "$filename" "$newname"
    done
    """

rule filter_ragtag:
    input: 
        "ragtag_{ref}_all/{sample}.ragtag_output/{sample}.scaffold.fasta"
    output:
        "ragtag_{ref}_all/{sample}.chr.fasta"
    shell: """
    cat {input} | seqkit replace -p '_RagTag' -r ''| seqkit grep -w 0 -vrnip 'Chr0' | seqkit sort -N  > {output} 
    """

rule ref_alignments_wholegenome:
    input:
        fasta="ragtag_{ref}_all/{sample}.chr.fasta"
    output:
        paf="ragtag_{ref}_all/all2ref/{sample}_{ref}.paf"
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell:
        "minimap2 -x asm5 -c --eqx --cs -t 20 -secondary=no {params.ref} {input.fasta} > {output.paf}"