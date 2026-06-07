import numpy as np
import pandas as pd
import os 

ref_paths = {
    'hg38_noALT': "/global/scratch/users/joana_rocha/PANPAN/reference/human_GRCh38.p14/hg38_HGSVC/hg38.no_alt.fa.gz",
    'clint': "/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fa",
    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    'mPanTro3' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta", ### no difference to thanksgiving/final primary
    'mPonAbe1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta",
    'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
}

ref_haps_paths = {
    'mPanTro3.hap1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.hap1.cur.20231031.fasta", ### no difference to thanksgiving/final primary
    'mPanPan1.hap1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.mat.cur.20231122.fasta",
    'mPanTro3.hap2' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.hap2.cur.20231031.fasta", ### no difference to thanksgiving/final primary
    'mPanPan1.hap2' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pat.cur.20231122.fasta"
}

references = ['ht2t', 'mPanTro3', 'mPanPan1']
reference_haps = list(ref_haps_paths.keys())

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
hifi_samples = sample_table['Specimen'].unique()
print(hifi_samples)

sample_table_hic = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples_hic.tsv', sep='\t')
hic_samples = sample_table_hic['specimen'].unique()
print(hic_samples)

verkko_samples=['PR01227', 'PR01008',  'PR01010', 'PR00366']
print(verkko_samples)

HAPLOTYPES = ["hap1", "hap2"]
GENOME_SIZES = [3100000000]

rule all:
    input:
        expand("assembly_stats/hifiasm/{gname}.{hap}.{gsize}.csv",
               gname=[s for s in hifi_samples if not s.startswith('AG')], hap=HAPLOTYPES, gsize=GENOME_SIZES),
        expand("assembly_stats/hifiasm_hic/{gname}.{hap}.{gsize}.csv",
               gname=hic_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        expand("assembly_stats/verkko/{gname}.{hap}.{gsize}.csv",
               gname=verkko_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        'assembly_stats/unified_assembly_stats.csv',
        expand('assembly_stats/references/{ref}.{gsize}.csv', ref=references, gsize=GENOME_SIZES),
        expand('assembly_stats/reference_haps/{ref}.{hap}.{gsize}.csv', ref=expand("{key}", key=ref_haps_paths.keys()), hap=HAPLOTYPES, gsize=GENOME_SIZES),
        'assembly_stats/allPan_references.tsv'

rule run_assembly_stats:
    input:
        lambda wildcards: (
            'Hifiasm-fasta_shortcut/{gname}.{hap}.p_ctg.fa'.format(**wildcards)
            if wildcards.type == 'hifiasm'
            else 'Hifiasm-fasta-HiC_shortcut/{gname}.{hap}.p_ctg.hic.fa'.format(**wildcards)
            if wildcards.type == 'hifiasm_hic'
            else 'Verkko-fasta_shortcut/{gname}.{hap}.verkko.fasta'.format(**wildcards)
            if wildcards.type == 'verkko'
            else "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/{ref}.{hap}.fasta.gz".format(ref=wildcards.ref, hap=wildcards.hap)
            if wildcards.type == 'reference_haps'
            else "invalid_input"
        )
    output:
        'assembly_stats/{type}/{gname}.{hap}.{gsize}.csv'
    params:
        genomename = lambda wildcards: '--genomename {gname}'.format(**wildcards)
    shell:
        '/global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats {input} {output} {params.genomename} --genomesize {wildcards.gsize}'

rule run_reference_stats:
    input:
        lambda wildcards: ref_paths[wildcards.ref]
    output:
        'assembly_stats/references/{ref}.{gsize}.csv'
    params:
        genomename = lambda wildcards: '--genomename {ref}'.format(**wildcards)
    shell:
        '/global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats {input} {output} {params.genomename} --genomesize {wildcards.gsize}'

rule run_reference_hap_stats:
    input:
        lambda wildcards: ref_haps_paths[wildcards.ref]
    output:
        'assembly_stats/reference_haps/{ref}.{hap}.{gsize}.csv'
    params:
        genomename = lambda wildcards: '--genomename {ref}'.format(**wildcards)
    shell:
        '/global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats {input} {output} {params.genomename} --genomesize {wildcards.gsize}'

rule concatenate_assembly_stats:
    input:
        hifiasm = expand("assembly_stats/hifiasm/{gname}.{hap}.{gsize}.csv",
               gname=[s for s in hifi_samples if not s.startswith('AG')], hap=HAPLOTYPES, gsize=GENOME_SIZES),
        hifiasm_hic = expand("assembly_stats/hifiasm_hic/{gname}.{hap}.{gsize}.csv", 
                             gname=hic_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        verkko = expand("assembly_stats/verkko/{gname}.{hap}.{gsize}.csv", 
                        gname=verkko_samples, hap=HAPLOTYPES, gsize=GENOME_SIZES),
        references = expand('assembly_stats/references/{ref}.{gsize}.csv', ref=references, gsize=GENOME_SIZES),
        reference_haps = expand('assembly_stats/reference_haps/{ref}.{hap}.{gsize}.csv', ref=expand("{key}", key=ref_haps_paths.keys()), hap=HAPLOTYPES, gsize=GENOME_SIZES)
    output:
        unified_output='assembly_stats/unified_assembly_stats.csv'
    run:
        import pandas as pd
        # Initialize an empty list to store dataframes
        dataframes = []

        # Define a mapping for method names
        method_mapping = {
            'hifiasm': 'HiFi',
            'hifiasm_hic': 'HiFi+HiC',
            'verkko': 'HiFi+HiC+ONT',
            'references': 'Reference Genome',
            'reference_haps': 'Reference Haplotype'
        }

        # Process each type
        for method_key, file_list in input.items():
            method_name = method_mapping[method_key]
            for file_path in file_list:
                df = pd.read_csv(file_path, sep='\t')
                df['type'] = file_path.split('/')[-1].split('.')[1]  # Extract hap type from filename
                df['method'] = method_name
                dataframes.append(df)

        # Concatenate all dataframes
        concatenated_df = pd.concat(dataframes, ignore_index=True)

        # Save the unified dataframe
        concatenated_df.to_csv(output.unified_output, sep='\t', index=False)

rule combine_all_stats:
    input:
        unified_assemblies = 'assembly_stats/unified_assembly_stats.csv'
    output:
        'assembly_stats/allPan_references.tsv'
    run:
        import pandas as pd

        # Read the unified assembly and reference statistics
        df_assemblies = pd.read_csv(input.unified_assemblies, sep='\t')

        # Save the combined DataFrame
        df_assemblies.to_csv(output[0], sep='\t', index=False)