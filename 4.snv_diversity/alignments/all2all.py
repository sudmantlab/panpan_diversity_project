import numpy as np
import pandas as pd
import os
import itertools
ref_paths = {
    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    'ht2t_ucsc': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/T2T_CHM13/hs1.fa",
    'mPanTro3' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta", ### no difference to thanksgiving/final primary
    'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"

}

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
sample_table_hprc=pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_names.txt', names=['names']) ### removed Y chromosome
sample_table_hic=pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples_hic.tsv', sep='\t') 

references =  ['ht2t_ucsc', 'ht2t', 'mPanTro3', 'mPanPan1'] 
samples = sample_table['Specimen'].unique()
human_samples =sample_table_hprc['names'].unique()
hic_samples = sample_table_hic['specimen'].unique()

haps = ["hap1", "hap2"]
sample_haps = expand("{sample}.{hap}", sample=samples, hap=haps)
sample_haps_hic = expand("{sample}.{hap}", sample=hic_samples, hap=haps)

# Prioritize hic samples when there's overlap
preferred_samples = list(set(sample_haps_hic).union(set(sample_haps)))
print(preferred_samples)

all_samples= sample_haps + list(human_samples)
print(all_samples)


rule all:
    input:
        expand('chains/{ref}/{anything}.chain', anything=human_samples, ref=['ht2t_ucsc']),
        expand('chains/{ref}/{anything}.chain', anything=sample_haps, ref=['mPanPan1'])
        #expand("{sample_pair[0]}/{sample_pair[0]}_{sample_pair[1]}.paf", sample_pair=sample_pairs_haps),

def get_fasta_path(wildcards):
    verkko_path = f"/global/scratch/users/joana_rocha/PANPAN/Verkko-fasta_shortcut/{wildcards.sample.replace('.hic', '')}.verkko.fasta"
    hic_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta-HiC_shortcut/{wildcards.sample.replace('.hic', '')}.p_ctg.hic.fa"
    hifi_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{wildcards.sample.replace('.hic', '')}.p_ctg.fa"
    
    # Check if verkko file exists
    if os.path.exists(verkko_path):
        return verkko_path
    # If verkko file does not exist, check if hic file exists
    elif os.path.exists(hic_path):
        return hic_path
    # If neither verkko nor hic files exist, return the non-hic path
    else:
        return hifi_path



rule run_rustybam:
    input:
        '{ref}/{anything}.paf',
    output:
        '{ref}/{anything}.rb.saffire.paf', 
    shell: """
    rb trim-paf {input} | rb break-paf --max-size 5000  | rb orient | rb filter --paired-len 100000 | rb stats --paf > {output}
    """

#trims back alignments that align the same query sequence more than once` \
#breaks the alignment into smaller pieces on indels of 5000 bases or more` \
#orients each contig so that the majority of bases are forward aligned` \
#rb filter --paired-len 100000 `#filters for query sequences that have at least 100,000 bases aligned to a target across all alignments.` \
#rb stats --paf `#calculates statistics from the trimmed paf file` \

rule all2all_haplotye_pairs:
    input:
        "/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{samplehap1}.p_ctg.fa",
        "/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{samplehap2}.p_ctg.fa"
    output:
        "{samplehap1}/{samplehap1}_{samplehap2}.paf"
    shell: """
    minimap2 -x asm20 -t 20 -c --eqx --cs -D -P –dual=no  {input[0]} {input[1]} > {output}
   """
