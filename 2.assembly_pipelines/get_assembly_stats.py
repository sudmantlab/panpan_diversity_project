import numpy as np
import pandas as pd
import os 

ref_paths = {
    'hg38_noALT': "/global/scratch/users/joana_rocha/PANPAN/reference/human_GRCh38.p14/hg38_HGSVC/hg38.no_alt.fa.gz",
    'clint': "/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fa",
    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    'mPanTro3' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta", 
    'mPonAbe1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta",
    'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
}

ref_haps_paths = {
    'mPanTro3.hap1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.hap1.cur.20231031.fasta", 
    'mPanPan1.hap1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.mat.cur.20231122.fasta",
    'mPanTro3.hap2' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.hap2.cur.20231031.fasta", 
    'mPanPan1.hap2' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pat.cur.20231122.fasta"
}

references = ['ht2t', 'mPanTro3', 'mPanPan1']
HAPLOTYPES = ["hap1", "hap2"]
GENOME_SIZES = [3100000000]

# Load samples
sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
hifi_samples = [s for s in sample_table['Specimen'].unique() if not s.startswith('AG')]

sample_table_hic = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples_hic.tsv', sep='\t')
hic_samples = sample_table_hic['specimen'].unique()

verkko_samples = ['PR01227', 'PR01008', 'PR01010', 'PR00366']

rule all:
    input:
        expand("assembly_stats/hifiasm/{gname}.{hap}.{gsize}.csv",
               gname=hifi_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        expand("assembly_stats/hifiasm_hic/{gname}.{hap}.{gsize}.csv",
               gname=hic_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        expand("assembly_stats/verkko/{gname}.{hap}.{gsize}.csv",
               gname=verkko_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        expand('assembly_stats/references/{ref}.{gsize}.csv', 
               ref=references, gsize=GENOME_SIZES),
        expand('assembly_stats/reference_haps/{ref_key}.{gsize}.csv', 
               ref_key=ref_haps_paths.keys(), gsize=GENOME_SIZES),
        'assembly_stats/unified_assembly_stats.csv',
        'assembly_stats/allPan_references.tsv'

# Standard assembly stats rule for hifiasm/verkko
rule run_assembly_stats:
    input:
        lambda wildcards: (
            'Hifiasm-fasta_shortcut/{gname}.{hap}.p_ctg.fa'.format(**wildcards)
            if wildcards.type == 'hifiasm'
            else 'Hifiasm-fasta-HiC_shortcut/{gname}.{hap}.p_ctg.hic.fa'.format(**wildcards)
            if wildcards.type == 'hifiasm_hic'
            else 'Verkko-fasta_shortcut/{gname}.{hap}.verkko.fasta'.format(**wildcards)
        )
    output:
        'assembly_stats/{type}/{gname}.{hap}.{gsize}.csv'
    params:
        genomename = lambda wildcards: '--genomename {gname}_{hap}'.format(**wildcards)
    shell:
        'python ./assemblystats.py {input} {output} {params.genomename} --genomesize {wildcards.gsize}'

# Stats for full references (ht2t, etc)
rule run_reference_stats:
    input:
        lambda wildcards: ref_paths[wildcards.ref]
    output:
        'assembly_stats/references/{ref}.{gsize}.csv'
    params:
        genomename = lambda wildcards: '--genomename {ref}'.format(**wildcards)
    shell:
        'python ./assemblystats.py {input} {output} {params.genomename} --genomesize {wildcards.gsize}'

# Stats for reference haplotypes (mPanTro3.hap1, etc)
rule run_reference_hap_stats:
    input:
        lambda wildcards: ref_haps_paths[wildcards.ref_key]
    output:
        'assembly_stats/reference_haps/{ref_key}.{gsize}.csv'
    params:
        genomename = lambda wildcards: '--genomename {ref_key}'.format(**wildcards)
    shell:
        'python ./assemblystats.py {input} {output} {params.genomename} --genomesize {wildcards.gsize}'

rule concatenate_assembly_stats:
    input:
        hifiasm = expand("assembly_stats/hifiasm/{gname}.{hap}.{gsize}.csv",
               gname=hifi_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        hifiasm_hic = expand("assembly_stats/hifiasm_hic/{gname}.{hap}.{gsize}.csv", 
                             gname=hic_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        verkko = expand("assembly_stats/verkko/{gname}.{hap}.{gsize}.csv", 
                        gname=verkko_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        references = expand('assembly_stats/references/{ref}.{gsize}.csv', 
                            ref=references, gsize=GENOME_SIZES),
        reference_haps = expand('assembly_stats/reference_haps/{ref_key}.{gsize}.csv', 
                                ref_key=ref_haps_paths.keys(), gsize=GENOME_SIZES)
    output:
        unified_output='assembly_stats/unified_assembly_stats.csv'
    run:
        dataframes = []
        method_mapping = {
            'hifiasm': 'HiFi',
            'hifiasm_hic': 'HiFi+HiC',
            'verkko': 'HiFi+HiC+ONT',
            'references': 'Reference Genome',
            'reference_haps': 'Reference Haplotype'
        }

        for method_key, file_list in input.items():
            method_name = method_mapping[method_key]
            for file_path in file_list:
                df = pd.read_csv(file_path, sep='\t')
                filename = os.path.basename(file_path)
                
                # Logic to assign haplotype type column
                if 'hap1' in filename or 'mat' in filename:
                    df['type'] = 'hap1'
                elif 'hap2' in filename or 'pat' in filename:
                    df['type'] = 'hap2'
                else:
                    df['type'] = 'primary'
                
                df['method'] = method_name
                dataframes.append(df)

        concatenated_df = pd.concat(dataframes, ignore_index=True)
        concatenated_df.to_csv(output.unified_output, sep='\t', index=False)

rule combine_all_stats:
    input:
        unified_assemblies = 'assembly_stats/unified_assembly_stats.csv'
    output:
        'assembly_stats/allPan_references.tsv'
    run:
        df = pd.read_csv(input.unified_assemblies, sep='\t')
        df.to_csv(output[0], sep='\t', index=False)