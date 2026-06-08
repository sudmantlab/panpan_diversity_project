import numpy as np
import pandas as pd
import os 

sample_table = pd.read_csv('../metadata/pepsamples_old.tsv', '\t')
print(sample_table)
sample_table = sample_table.set_index("SampleID", drop=False)
print(sample_table)
#raise Exception()
samples = sample_table['Specimen'].unique()
print(samples)
sample_runs = sample_table['SampleID'].unique()
print(sample_runs)


rule all:
    input:
        'coverage_stats/deepTool_plots/coverage_plot',
        #expand('assembly/{sample}.p_ctg.fa', sample=samples),
        #expand('reads/{sample_run}.ccs.filt.fastq.gz', sample_run=sample_runs),
        #expand('bam_mapped/{sample}/{sample_run}.ccs.flagfilt.bam', zip, sample=sample_table.Specimen, sample_run=sample_table.SampleID),
        #expand('bam_final/{sample}.sorted.bam', sample=samples),
        #expand('coverage_stats/{sample}.default_hist.txt', sample=samples),
        #expand('coverage_stats/{sample}.default_hist_renamed.txt', sample=samples),
        #expand('coverage_stats/{sample}.bg_hist.txt', sample=samples),
        #expand('coverage_stats/{sample}.bg_hist_renamed.txt', sample=samples),
        #expand('coverage_stats/{sample}.bga_hist.txt', sample=samples),
        #expand('coverage_stats/{sample}.bga_hist_renamed.txt', sample=samples),
        
        
rule symlink_refs:
    input: '/global/scratch/users/joana_rocha/PANPAN/output/hifiasm-fasta/{sample}/joana_settings_shortcut/no_opts/{sample}.p_ctg.fa'
    output: 'assembly/{sample}.p_ctg.fa'
    shell: """
    ln -s {input} {output} 
    """


def get_input_symlink_reads(wildcards):
    row = sample_table.loc[wildcards.sample_run]
    return f'/global/scratch/users/joana_rocha/PANPAN/output/HiFi-adapterFiltered/{row.Specimen}/joana_settings_shortcut/{row.Facility}/{row.Lane}/{row.SampleID}.ccs.filt.fastq.gz'


rule symlink_reads:
    input: get_input_symlink_reads
    output: 'reads/{sample_run}.ccs.filt.fastq.gz'
    shell: """
    ln -s {input} {output} 
    """

rule minimap2_fastqs_ref:
    input:
        'assembly/{sample}.p_ctg.fa',
        'reads/{sample_run}.ccs.filt.fastq.gz',
    output: 'bam_mapped/{sample}/{sample_run}.ccs.flagfilt.bam'
    shell: 'minimap2 -ax map-hifi -t 20 {input[0]} {input[1]} | samtools view -q 10  -bT {input[0]}  -o {output}'

# sam flags potentially relevant to be filtered out (if any of those is present):
# 4     0x4     unmapped reads
# 256   0x100   alignment not primary
# 512   0x200   fail platform QC
# total: -F 772

rule minimap2_sort:
    input: 'bam_mapped/{sample}/{sample_run}.ccs.flagfilt.bam'
    output: temp('bam_mapped/{sample}/{sample_run}.flagfilt.sorted.bam')
    shell: 'samtools sort -o {output} {input}'

#no need to remove duplicates because there are no optical duplicates
rule add_readgroup:
    input: 'bam_mapped/{sample}/{sample_run}.flagfilt.sorted.bam'
    output: temp('bam_mapped/{sample}/{sample_run}.flagfilt.sorted.rg.bam')
    shell: 'samtools addreplacerg -r ID:{wildcards.sample_run} -r PL:PACBIO -o {output} {input}'


def get_input_merged(wildcards):
    return expand(
        'bam_mapped/{sample}/{sample_run}.flagfilt.sorted.rg.bam',
        sample = [wildcards.sample],
        sample_run = sample_table[sample_table['Specimen']==wildcards.sample]['SampleID'].values
        )

rule merge_to_samples:
    input:  get_input_merged 
    output: temp('bam_merged/{sample}.flagfilt.sorted.rg.merged.bam')
    shell: 'samtools merge -f {output} {input}'

rule sort_merged:
    input: 'bam_merged/{sample}.flagfilt.sorted.rg.merged.bam'
    output: 'bam_final/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """

rule estimate_coverage:
    input: 'bam_final/{sample}.sorted.bam'
    output: temp('coverage_stats/{sample}.default_hist.txt')
    shell: 'samtools view -b {input} | genomeCoverageBed -ibam - > {output}'

rule concat_coverage:
    input: 'coverage_stats/{sample}.default_hist.txt'
    output: 'coverage_stats/{sample}.default_hist_renamed.txt'
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.sample}"}}' | sed 's| |\t|g' >  {output}
    """
#cat *renamed.txt | sed 's| |\t|g'  > all.default_hist.txt 
rule estimate_coverage_d:
    input: 'bam_final/{sample}.sorted.bam'
    output: 'coverage_stats/{sample}.d_hist.txt'
    shell: 'samtools view -b {input} | genomeCoverageBed -d -ibam - > {output}'

rule concat_coverage_d:
    input: 'coverage_stats/{sample}.d_hist.txt'
    output: 'coverage_stats/{sample}.d_hist_renamed.txt'
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.sample}"}}' | sed 's| |\t|g' >  {output}
    """

rule estimate_coverage_bg:
    input: 'bam_final/{sample}.sorted.bam'
    output: 'coverage_stats/{sample}.bg_hist.txt'
    shell: 'samtools view -b {input} | genomeCoverageBed -bg -ibam - > {output}'

rule concat_coverage_bg:
    input: 'coverage_stats/{sample}.bg_hist.txt'
    output: 'coverage_stats/{sample}.bg_hist_renamed.txt'
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.sample}"}}' | sed 's| |\t|g' >  {output}
    """


rule estimate_coverage_bga:
    input: 'bam_final/{sample}.sorted.bam'
    output: 'coverage_stats/{sample}.bga_hist.txt'
    shell: 'samtools view -b {input} | genomeCoverageBed -bga -ibam - > {output}'

rule concat_coverage_bga:
    input: 'coverage_stats/{sample}.bga_hist.txt'
    output: 'coverage_stats/{sample}.bga_hist_renamed.txt'
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.sample}"}}' | sed 's| |\t|g' >  {output}
    """

rule plotcoverage:
    input: expand('bam_final/{sample}.sorted.bam', sample=samples)
    output: 'coverage_stats/deepTool_plots/coverage_plot'
    shell: """
    plotCoverage --bamfiles {input[0]} --plotFile {output[0]} -n 1000000 --plotTitle "Minimap to self" --outRawCounts coverage.tab --ignoreDuplicates --minMappingQuality 10 
    """
#plotCoverage -b bam_final/PR01227.sorted.bam bam_final/PR01223.sorted.bam bam_final/PR00512.sorted.bam bam_final/PR01225.sorted.bam bam_final/PR01228.sorted.bam bam_final/PR00834.sorted.bam bam_final/PR00838.sorted.bam bam_final/PR01100.sorted.bam bam_final/PR00366.sorted.bam bam_final/PR00445.sorted.bam bam_final/PR00249.sorted.bam  --plotFile coverage_stats/deepTool_plots/coverage_plot -n 1000000 --plotTitle "Minimap to self" --outRawCounts coverage.tab --ignoreDuplicates --minMappingQuality 10 