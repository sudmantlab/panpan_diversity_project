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
#sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/bonobosamples.tsv', sep='\t')
#sample_table_hprc = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_names.txt', names=['names'])
#sample_table_hgsvc = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_names.txt', names=['names'])
#sample_table_hic = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples_hic.tsv', sep='\t')

# Extract unique samples
references = ['ht2t', 'mPanTro3', 'mPanPan1', 'mPonAbe1']
#samples = sample_table['Specimen'].unique()
#samples = ['mPanPan1', 'mPanTro3']
#samples = ['PR001227']
#human_samples = sample_table_hprc['names'].unique()
#human_samples = sample_table_hgsvc['names'].unique()
#human_samples = sample_table_hgsvc['names'].tolist()
#human_samples = ['HG02554.1', 'HG02554.2']
# Generate sample-haplotype combinations
#haps = ["hap1", "hap2"]
#sample_haps = expand("{sample}.{hap}", sample=samples, hap=haps)


# Print sample haplotypes for debugging

# print(sample_haps)
#print(human_samples)

rule all:
    input:
        # Add both rules to the all rule to ensure they are both run
       #directory(expand('ragtag_{ref}_all/{sample_hap}.ragtag_output', sample_hap=sample_haps, ref=['mPanTro3'])),
       # expand('ragtag_{ref}_all/{anything}.ragtag_output', anything=sample_haps, ref=['mPanPan1']),
       #expand("ragtag_{ref}_all/{anything}.ragtag_output", anything=human_samples, ref=['ht2t']),
        #expand("ragtag_{ref}_all/{sample}.{hap}.chr.fasta",  sample=samples, hap=haps, ref=['mPonAbe1']),
        #expand("ragtag_{ref}_all/{sample}.{hap}.chr.fasta",  sample=samples, hap=haps, ref=['mPanPan1']),
        #expand("ragtag_{ref}_all/{samplehap}.chr.fasta",  samplehap=human_samples,  ref=['ht2t']),
        #expand("ragtag_{ref}_all/all2ref/{samplehap}_{ref}.paf", samplehap=sample_haps, ref=['mPonAbe1']),
        #expand("ragtag_{ref}_all/all2ref/{samplehap}_{ref}.paf", samplehap=sample_haps, ref=['mPanTro3']),
        #expand("ragtag_{ref}_all/all2ref/{samplehap}_{ref}.paf", samplehap=sample_haps, ref=['mPanPan1']),
        #expand("ragtag_{ref}_all/all2ref_noeqx/{samplehap}_{ref}.paf",  samplehap=sample_haps, ref=['mPanPan1']),
        #expand("ragtag_{ref}_all/all2ref_noeqx/{samplehap}_{ref}.paf", samplehap=human_samples, ref=['ht2t']),
        
def get_fasta_path(wildcards):
    sample_hap = f"{wildcards.sample}.{wildcards.hap}"
    verkko_path = f"/global/scratch/users/joana_rocha/PANPAN/Verkko-fasta_shortcut/{sample_hap.replace('.hic', '')}.verkko.fasta"
    hic_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta-HiC_shortcut/{sample_hap.replace('.hic', '')}.p_ctg.hic.fa"
    hifi_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{sample_hap.replace('.hic', '')}.p_ctg.fa"
    
    if os.path.exists(verkko_path):
        return verkko_path
    elif os.path.exists(hic_path):
        return hic_path
    else:
        return hifi_path

rule panpan2ref:
    input:
        #fasta=get_fasta_path
        fasta='{sample}.{hap}.verkko.fasta'
    output: 
        directory('ragtag_{ref}_all/{sample}.{hap}.ragtag_output')
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    ragtag.py scaffold -C -u --mm2-params '-x asm5' --unimap-params '-x asm5 -t 20' -o {output} {params.ref} {input.fasta} &&
    cd {output} &&
    for filename in ragtag.*; do
        newname="$(echo "$filename" | sed 's/ragtag/{wildcards.sample}.{wildcards.hap}/')"
        echo "Renaming $filename to $newname"
        mv "$filename" "$newname"
    done
    """



rule hprc2ref:
    input: "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_assemblies/{human_sample}.fa"
    output: directory('ragtag_{ref}_all/{human_sample}.ragtag_output')
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    ragtag.py scaffold -C -u --mm2-params '-x asm5' --unimap-params '-x asm5 -t 20' -o {output} {params.ref} {input} &&
    cd ragtag_{wildcards.ref}_all/{wildcards.human_sample}.ragtag_output/ &&
    for filename in ragtag.* ; do mv "$filename" "$(echo "$filename" | sed 's/ragtag/{wildcards.human_sample}/')"; done;
    """

rule hgsvc2ref:
    input: "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_assemblies/{human_sample}.fasta.gz"
#    input: "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HG02554_scott/genome/{human_sample}.fasta"
    output: directory('ragtag_{ref}_all/{human_sample}.ragtag_output')#
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    ragtag.py scaffold -C -u --mm2-params '-x asm5' --unimap-params '-x asm5 -t 20' -o {output} {params.ref} {input} &&
    cd ragtag_{wildcards.ref}_all/{wildcards.human_sample}.ragtag_output/ &&
    for filename in ragtag.* ; do mv "$filename" "$(echo "$filename" | sed 's/ragtag/{wildcards.human_sample}/')"; done;
    """


rule decompress_and_index:
    input: 
        "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_assemblies/{human_sample}.fasta.gz"
       # "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HG02554_scott/genome/{human_sample}.fasta"
    output: 
        fasta=temp("temp/{human_sample}.fasta"),
        fai=temp("temp/{human_sample}.fasta.fai")  # Index file
    shell: 
        """
        # Create temp dir if needed
        mkdir -p temp
        
        # Decompress and validate
        cat {input} > {output.fasta}
        
        # Index the FASTA
        samtools faidx {output.fasta}
        
        # Verify both files exist
        if [[ ! -s {output.fasta} || ! -s {output.fai} ]]; then
            echo "ERROR: Decompression/indexing failed" >&2
            exit 1
        fi
        """

rule hgsvc2ref:
    input: "temp/{human_sample}.fasta"
    output: directory('ragtag_{ref}_all/{human_sample}.ragtag_output')

    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    ragtag.py scaffold -C -u --mm2-params '-x asm5' --unimap-params '-x asm5 -t 8' -o {output} {params.ref} {input} &&
    cd ragtag_{wildcards.ref}_all/{wildcards.human_sample}.ragtag_output/ &&
    for filename in ragtag.* ; do mv "$filename" "$(echo "$filename" | sed 's/ragtag/{wildcards.human_sample}/')"; done;
    """

#rule hgsvc2ref:
#    input:
#        fasta="temp/{human_sample}.fasta",
#        ref=ref_paths["ht2t"]
#    output:
#        "ragtag_{ref}_all/{human_sample}.ragtag_output/{human_sample}.scaffold.fasta"
#    params:
#        outdir=lambda wildcards: f"ragtag_{wildcards.ref}_all/{wildcards.human_sample}.ragtag_output"
#    log:
#        "logs/ragtag_{ref}_{human_sample}.log"
#    shell:
#        """
#        mkdir -p {params.outdir}
#        ragtag.py scaffold -C -u \
#            --mm2-params '-x asm5 -t 8' \
#            --unimap-params '-x asm5 -t 8' \
#            -o {params.outdir} {input.ref} {input.fasta} > {log} 2>&1
            
        # Verify output exists
#        if [[ ! -f "{params.outdir}/ragtag.scaffold.fasta" ]]; then
#            echo "Error: ragtag.scaffold.fasta not found" >> {log}
#            exit 1
#        fi
#        
#        # Rename the output file
#        mv "{params.outdir}/ragtag.scaffold.fasta" "{output}"
#        """

rule filter_ragtag:
    input: 
        "ragtag_{ref}_all/{samplehap}.ragtag_output/{samplehap}.scaffold.fasta",
    output:
        "ragtag_{ref}_all/{samplehap}.chr.fasta",
    shell: """
    cat {input} | seqkit replace -p '_RagTag' -r ''| seqkit grep -w 0 -vrnip 'Chr0' | seqkit sort -N  > {output} 
    """

rule ref_alignments_wholegenome:
    input:
        fasta="ragtag_{ref}_all/{samplehap}.chr.fasta"
    output:
        paf="ragtag_{ref}_all/all2ref/{samplehap}_{ref}.paf"
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell:
        "minimap2 -x asm5 -c --cs --eqx -t 20 -secondary=no {params.ref} {input.fasta} > {output.paf}"


rule alignments_ref_wholegenome:
    input:
        fasta="ragtag_{ref}_all/{samplehap}.chr.fasta"
    output:
        paf="ragtag_{ref}_all/all2ref_noeqx/{samplehap}_{ref}.paf"
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell:
        "minimap2 -cx asm5 -t 20 {input.fasta} {params.ref}  > {output.paf}"
