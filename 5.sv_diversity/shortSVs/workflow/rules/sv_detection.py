import numpy as np
import pandas as pd
import os 

BASEDIR = '/global/scratch/users/nicolas931010/sv_detection'
BAMDIR = '/global/scratch/users/joana_rocha/PANPAN/minimap2others'
VCFDIR = '/global/scratch/users/joana_rocha/PANPAN/PANPAN_svs/svs'
SAMPLE_TABLE = pd.read_csv(BASEDIR + '/docs/sample_table_combined.tsv', sep='\t')
SAMPLES = SAMPLE_TABLE['sample_id'].unique()
REFS = [ 'hg38.no_alt', 'ht2t', 'clint', 'mPanTro3', 'mPanPan1', 'mPonAbe1' ]
DATASETS = [ 'hprc', 'panpan', 'panpan-pt', 'panpan-pp', 'combined' ]
ID_DATASET_DICT = SAMPLE_TABLE.set_index('sample_id')['dataset'].to_dict()

def get_bam_from_sample(wildcards):
    return f"{BAMDIR}/{ID_DATASET_DICT[wildcards.sample]}_winnowmaping/{wildcards.ref}/bam_final/{wildcards.sample}.sorted.bam"

def get_snfs_from_dataset(wildcards):
    sample_table = pd.read_csv(BASEDIR + f'/docs/sample_table_{wildcards.dataset}.tsv', sep='\t')
    snfs = sample_table.apply(lambda row: f'{BASEDIR}/sniffles/{wildcards.ref}/{row["sample_id"]}.snf', axis=1)
    return snfs.tolist()

def get_vcf_from_dataset(wildcards):
    sample_table = pd.read_csv(BASEDIR + f'/docs/sample_table_{wildcards.dataset}.tsv', sep='\t')
    snfs = sample_table.apply(lambda row: f'{VCFDIR}/{wildcards.ref}/{row["sample_id"]}/{row["sample_id"]}.sorted.vcf.gz', axis=1)
    return snfs.tolist()

def get_repeatmasker_database_from_dataset(wildcards):
    sample_table = pd.read_csv(BASEDIR + f'/docs/sample_table_{wildcards.dataset}.tsv', sep='\t')
    repeatmasker_database = sample_table.repeatmasker_database.unique().tolist()
    return repeatmasker_database

wildcard_constraints:
    sample = '|'.join([x for x in SAMPLES]),
    ref = '|'.join([x for x in REFS]),
    dataset = '|'.join([x for x in DATASETS]),

rule all:
    input:
        # expand(BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.pca.done', ref=[ 'ht2t' ], dataset='hprc'),
        # expand(BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.pca.done', ref=[ 'mPanTro3' ], dataset='panpan-pt'),
        # expand(BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.pca.done', ref=[ 'mPanPan1' ], dataset='panpan-pp'),
        # expand(BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.pca.done', ref=[ 'mPonAbe1' ], dataset='combined'),
        ## the following three lines are temporary
        expand(BASEDIR + '/{caller}/{ref}/{dataset}.filtered.done', ref=[ 'mPanPan1' ], dataset='panpan-pp', caller=[ 'sniffles', 'svim-asm']),
        expand(BASEDIR + '/{caller}/{ref}/{dataset}.filtered.done', ref=[ 'ht2t' ], dataset='hprc', caller=[ 'sniffles', 'svim-asm']),
        expand(BASEDIR + '/{caller}/{ref}/{dataset}.filtered.done', ref=[ 'mPanTro3' ], dataset='panpan-pt', caller=[ 'sniffles', 'svim-asm']),
        expand(BASEDIR + '/truvari/{ref}/{dataset}.pca_with_consensus_callset.done', ref=[ 'mPanPan1' ], dataset='panpan-pp'),
        expand(BASEDIR + '/truvari/{ref}/{dataset}.pca_with_consensus_callset.done', ref=[ 'ht2t' ], dataset='hprc'),
        expand(BASEDIR + '/truvari/{ref}/{dataset}.pca_with_consensus_callset.done', ref=[ 'mPanTro3' ], dataset='panpan-pt'),

include: "sniffles.py"
include: "truvari.py"
