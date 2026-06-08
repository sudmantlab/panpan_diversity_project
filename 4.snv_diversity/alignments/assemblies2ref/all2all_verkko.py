import numpy as np
import pandas as pd
import os
import itertools
from itertools import combinations
import glob
import gzip
import csv
import sys
import logging

ref_paths = {
#    'hg38_noALT': "/global/scratch/users/joana_rocha/PANPAN/reference/human_GRCh38.p14/hg38_HGSVC/hg38.no_alt.fa.gz",
#    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    'mPanTro3' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta",
#    'mPonAbe1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta.gz",
#    'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta.gz"
}

references = ['mPanTro3']
sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/chimpsamples.tsv', sep='\t')
samples = sample_table['Specimen'].unique()
print(samples)
haps = ["hap1", "hap2"]
sample_haps = [f"{sample}.{hap}" for sample in samples for hap in haps]
unique_pairs = list(combinations(sample_haps, 2))

print(sample_haps)

chromo_list = pd.read_csv('chromosomes/mPanTro3.txt', header=None) 
chromo_names = chromo_list[0].tolist()

regions = ["chr3_hap1_hsa4:144659238-146113772"]

#def get_fasta(sample, hap):
#    pattern = f"{sample}.{hap}.*.fasta"
#    fasta_files = glob.glob(pattern)
#    if not fasta_files:
#        raise ValueError(f"No uncompressed FASTA files found for sample {sample} and haplotype {hap}.")
#    return fasta_files[0]

def get_fasta(sample, hap):
    pattern = f"{sample}.{hap}.*.fa*" #to include non-verkko outputs
    fasta_files = glob.glob(pattern)
    if not fasta_files:
        raise ValueError(f"No uncompressed FASTA files found for sample {sample} and haplotype {hap}.")
    return fasta_files[0]

def extract_chromosome_names(fasta_path):
    chrom_names = []
    opener = gzip.open if fasta_path.endswith('.gz') else open
    with opener(fasta_path, 'rt') as fasta:
        for line in fasta:
            if line.startswith('>'):
                chrom_name = line.split()[0][1:]
                chrom_names.append(chrom_name)
    return chrom_names

def get_fasta_paths(wildcards):
    # Extract sample and hap identifiers from the wildcards
    samplehap1, hap1 = wildcards.samplehap1.split('.')
    samplehap2, hap2 = wildcards.samplehap2.split('.')
    
    # Construct file paths
    hap1_fasta = f"ragtag_{wildcards.ref}/{samplehap1}.{hap1}.chr.fasta"
    hap2_fasta = f"ragtag_{wildcards.ref}/{samplehap2}.{hap2}.chr.fasta"
    return hap1_fasta, hap2_fasta

ruleorder:  split_by_chromosome > subset_byregion

rule all:
    input:
        directory(expand("chimps_ragtag_{ref}/{sample}.{hap}.ragtag_output", sample=samples, hap=haps, ref=references)),
        #expand("chimps_ragtag_{ref}/{sample}.{hap}.chr.fasta",  sample=samples, hap=haps, ref=references),
        #expand("ragtag_{ref}/all2ref/{samplehap}_{ref}.paf", samplehap=sample_haps, ref=references),
        #expand("ragtag_{ref}/all2ref_renamed/{samplehap}_{ref}.paf", samplehap=sample_haps, ref=references),
        #expand("ragtag_{ref}/{chrom}/{sample}.{hap}.{chrom}.fasta",  sample=samples, hap=haps, ref=references, chrom=chromo_names),
        #expand("ragtag_{ref}/{region}/{samplehap}.{region}.fasta", samplehap=sample_haps, ref=references, region=regions),
    #expand("ragtag_{ref}/all2all/{pair[0]}_{pair[1]}.paf", ref=references, pair=unique_pairs),
    # expand("ragtag_{ref}/all2all_renamed/{pair[0]}_{pair[1]}.paf", ref=references, pair=unique_pairs),
        #expand("ragtag_{ref}/{region}/all2ref_renamed/{samplehap}_{ref}.paf", samplehap=sample_haps, region=regions, ref=references),
        #expand("ragtag_{ref}/{region}/all2ref_renamed/{samplehap}_{ref}.paf", samplehap=sample_haps, region=regions, ref=references),
        #expand("ragtag_{ref}/{region}/all2ref_renamed/{samplehap}_{ref}.paf", samplehap=sample_haps, region=chromo_names, ref=references),
        #'ragtag_mPanTro3/chr3_hap1_hsa4:144659238-146113772/all2all_renamed/ava_concatenated_chr3_hap1_hsa4:144659238-146113772.paf',
        #'ragtag_mPanTro3/chr3_hap1_hsa4:145709505-145952261/all2all_renamed/ava_concatenated_chr3_hap1_hsa4:145709505-145952261.paf',
       
       

rule write_chromosome_names:
    input:
        lambda wildcards: ref_paths[wildcards.ref]
    output:
        "chromosomes/{ref}.txt"
    run:
        chrom_names = extract_chromosome_names(input[0])
        with open(output[0], 'w') as f:
            for name in chrom_names:
                f.write(f"{name}\n")

rule ragtag_ref:
    input:
        fasta=lambda wildcards: get_fasta(wildcards.sample, wildcards.hap)
    output: 
        directory('chimps_ragtag_{ref}/{sample}.{hap}.ragtag_output')
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    set -euo pipefail
    echo "Running ragtag.py scaffold for sample {wildcards.sample}.{wildcards.hap} with reference {params.ref}"
    ragtag.py scaffold -C -u --mm2-params '-x asm5' --unimap-params '-x asm5 -t 20' -o {output} {params.ref} {input}
    if [ $? -eq 0 ]; then
        echo "ragtag.py scaffold completed successfully for sample {wildcards.sample}.{wildcards.hap}"
        cd {output}
        for filename in ragtag.*; do 
            mv "$filename" "$(echo "$filename" | sed 's/ragtag/{wildcards.sample}.{wildcards.hap}/')"
        done
    else
        echo "ragtag.py scaffold failed for sample {wildcards.sample}.{wildcards.hap}"
        exit 1
    fi
    """
## -C ==  concatenate unplaced contigs and make 'chr0'

### filter the ragtagged genomes to only include chromosomes, and make sure that the name of the scaffolds match
rule filter_ragtag:
    input: 
        "chimps_ragtag_{ref}/{sample}.{hap}.ragtag_output/{sample}.{hap}.scaffold.fasta",
    output:
        "chimps_ragtag_{ref}/{sample}.{hap}.chr.fasta",
    shell: """
    cat {input} | seqkit replace -p '_RagTag' -r ''| seqkit grep -w 0 -vrnip 'Chr0' | seqkit sort -N  > {output} 
    """
#seqkit seq -m 5000000 will likeley remove sex chroms
#cat {input} | seqkit replace -p '_RagTag' -r ''| seqkit grep -w 0 -vrnip 'Chr0' | seqkit sort -N  > {output} 

rule ref_alignments_wholegenome:
    input:
        fasta="chimps_ragtag_{ref}/{samplehap}.chr.fasta"
    output:
        paf="chimps_ragtag_{ref}/all2ref/{samplehap}_{ref}.paf"
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell:
        "minimap2 -x asm5 -c --eqx --cs -t 20 -secondary=no {params.ref} {input.fasta} > {output.paf}"


rule process_wholegenome_all2ref_paf:
    input:
        paf="chimps_ragtag_{ref}/all2ref/{samplehap}_{ref}.paf"
    output:
        processed_paf="chimps_ragtag_{ref}/all2ref_renamed/{samplehap}_{ref}.paf"
    run:
        csv.field_size_limit(sys.maxsize)
        with open(input.paf, 'r') as infile, open(output.processed_paf, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            for row in reader:
                # Check if the row has at least 12 columns to avoid IndexError and ensure qname and tname are the same chromosome
                if len(row) >= 12:
                    row[0] = f"{wildcards.samplehap}_{row[0]}"
                    row[5] = f"{wildcards.ref}_{row[5]}"
                    writer.writerow(row)


rule split_by_chromosome:
    input:
        "ragtag_{ref}/{sample}.{hap}.chr.fasta"
    output:
        "ragtag_{ref}/{chrom}/{sample}.{hap}.{chrom}.fasta"
    shell: """
    seqkit grep -w 0 -rnip '{wildcards.chrom}' {input} > {output} 
    """

rule subset_byregion:
    input:
        fasta_file = "ragtag_{ref}/{sample}.{hap}.chr.fasta",
    output:
        subset_fasta = "ragtag_{ref}/{region}/{sample}.{hap}.{region}.fasta"
    params:
        target = "{region}"
    shell:
        """
        samtools faidx {input.fasta_file} {params.target} > {output.subset_fasta}
        """

rule cross_alignments_byregion:
    input:
        #fasta1=lambda wildcards: f"ragtag_{{ref}}/{{region}}/{wildcards.samplehap1.split('.')[0]}.{wildcards.samplehap1.split('.')[1]}.{{region}}.fasta",
        fasta1=lambda wildcards: f"ragtag_{{ref}}/{wildcards.samplehap1.split('.')[0]}.{wildcards.samplehap1.split('.')[1]}.chr.fasta",
        #fasta2=lambda wildcards: f"ragtag_{{ref}}/{{region}}/{wildcards.samplehap2.split('.')[0]}.{wildcards.samplehap2.split('.')[1]}.{{region}}.fasta"
        fasta2=lambda wildcards: f"ragtag_{{ref}}/{wildcards.samplehap2.split('.')[0]}.{wildcards.samplehap2.split('.')[1]}.chr.fasta"
    output:
        #paf="ragtag_{ref}/{region}/all2all/{samplehap1}_{samplehap2}.paf"
        paf="ragtag_{ref}/all2all/{samplehap1}_{samplehap2}.paf"
    shell:
        "minimap2 -x asm5 -c --eqx --cs -t 20 {input.fasta1} {input.fasta2} > {output.paf}"


#rule process_haplotype_pairs_paf:
    #input:
        #paf="ragtag_{ref}/{region}/all2all/{samplehap1}_{samplehap2}.paf"
     #   paf="ragtag_{ref}/all2all/{samplehap1}_{samplehap2}.paf"
    #output:
        #processed_paf="ragtag_{ref}/{region}/all2all_renamed/{samplehap1}_{samplehap2}.paf"
     #   processed_paf="ragtag_{ref}/all2all_renamed/{samplehap1}_{samplehap2}.paf"
    #run:
        #csv.field_size_limit(sys.maxsize)

        #with open(input.paf, 'r') as infile, open(output.processed_paf, 'w', newline='') as outfile:
            #reader = csv.reader(infile, delimiter='\t')
            #writer = csv.writer(outfile, delimiter='\t')
            
           # for row in reader:
                # Check if the row has at least 6 columns to avoid IndexError and ensure qname and tname are the same chromosome
                #if len(row) >= 6 and row[0] == row[5]:
                  #  row[0] = f"{wildcards.samplehap2}_{row[0]}"
                   # row[5] = f"{wildcards.samplehap1}_{row[5]}"
                   # writer.writerow(row)


rule process_haplotype_pairs_paf:
    input:
        paf="ragtag_{ref}/all2all/{samplehap1}_{samplehap2}.paf"
    output:
        processed_paf="ragtag_{ref}/all2all_renamed/{samplehap1}_{samplehap2}.paf"
    run:
        csv.field_size_limit(sys.maxsize)

        with open(input.paf, 'r') as infile, open(output.processed_paf, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')

            for row in reader:
                if len(row) >= 12:  # Check if the row has at least 12 columns
                    row[0] = f"{wildcards.samplehap2}_{row[0]}"
                    row[5] = f"{wildcards.samplehap1}_{row[5]}"
                    writer.writerow(row)
                else:
                    print(f"Skipping row: {row}")

 
rule ava_paf_byregion:
    input:
        paf_files=expand("ragtag_{ref}/{region}/all2all_renamed/{pair[0]}_{pair[1]}.paf", pair=unique_pairs, region=regions, ref=references)
    output:
        'ragtag_{ref}/{region}/all2all_renamed/ava_concatenated_{region}.paf'
    shell:
        'cat {input.paf_files} > {output}'


rule ref_alignments_byregion:
    input:
        fasta=lambda wildcards: f"ragtag_{{ref}}/{{region}}/{wildcards.samplehap.split('.')[0]}.{wildcards.samplehap.split('.')[1]}.{{region}}.fasta"
    output:
        paf="ragtag_{ref}/{region}/all2ref/{samplehap}_{ref}.paf"
    params:
        ref=lambda wildcards: ref_paths[wildcards.ref]
    shell:
        "minimap2 -x asm5 -c --eqx --cs -t 20 -secondary=no {params.ref} {input.fasta} > {output.paf}"


rule process_all2ref_paf:
    input:
        paf="ragtag_{ref}/{region}/all2ref/{samplehap}_{ref}.paf"
    output:
        processed_paf="ragtag_{ref}/{region}/all2ref_renamed/{samplehap}_{ref}.paf"
    run:
        csv.field_size_limit(sys.maxsize)
        with open(input.paf, 'r') as infile, open(output.processed_paf, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            for row in reader:
                # Check if the row has at least 12 columns to avoid IndexError and ensure qname and tname are the same chromosome
                if len(row) >= 12:
                    row[0] = f"{wildcards.samplehap}_{row[0]}"
                    row[5] = f"{wildcards.ref}_{row[5]}"
                    writer.writerow(row)
