import os
import pandas as pd
import numpy as np
import itertools
import gzip

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
chimp_sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/chimpsamples.tsv', sep='\t')
bonobo_sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/bonobosamples.tsv', sep='\t')


sample_table_hprc=pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_names.txt', names=['names']) ### removed Y chromosome
sample_table_hprc['base_sample'] = sample_table_hprc['names'].apply(lambda x: x.split('.')[0])


references = ['mPanTro3']

samples = sample_table['Specimen'].unique()
chimpsamples = chimp_sample_table['Specimen'].unique()
bonobosamples = bonobo_sample_table['Specimen'].unique()

print(samples)
human_samples = sample_table_hprc['base_sample'].unique()
print(human_samples)


haps = ["hap1", "hap2"]
sample_haps = expand("{sample}.{hap}", sample=samples, hap=haps)

human_sample_haps =sample_table_hprc['names'].unique()

print(sample_haps)
print(human_sample_haps)

rule all:
    input:
        #expand("vcfs/ht2t/diploid/{sample}}_unphased.vcf.gz", sample=human_samples),
        #'vcfs/ht2t/diploid/hprc_mapped2ht2t_unphased.sorted.vcf.gz',
        #'plink/diploid/hprc_mapped2ht2t_unphased.eigenval', 
        #'plink/diploid/hprc_mapped2ht2t_unphased.eigenvec',
        expand("vcfs/mPanTro3/diploid/{sample}_unphased.vcf.gz", sample=chimpsamples),
        'vcfs/mPanTro3/diploid/pantros_mapped2mPanTro3.sorted_unphased.vcf.gz',
        'plink/diploid/pantros_mapped2mPanTro3_unphased.eigenval', #pantros (no bonobos)
        'plink/diploid/pantros_mapped2mPanTro3_unphased.eigenvec', #pantros (no bonobos)
  


rule process_panpan_vcf_diploid_unphased:
    input:
        "vcfs/{ref}/temp_{sample}.vcf.gz"
    output:
        temp("vcfs/{ref}/diploid/{sample}_unphased.vcf")
    run:
        input_vcf = input[0]
        output_vcf = output[0]
        hap1_col = wildcards.sample + ".hap1"
        hap2_col = wildcards.sample + ".hap2"
        dip_col = wildcards.sample

        # Read the input VCF file
        with gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r') as file:
            lines = file.readlines()

        # Extract header and data lines
        header_lines = [line for line in lines if line.startswith('#')]
        data_lines = [line for line in lines if not line.startswith('#')]

        # Create a dataframe from the data lines
        data = [line.strip().split('\t') for line in data_lines]
        columns = header_lines[-1].strip().split('\t')
        df = pd.DataFrame(data, columns=columns)

        # Function to apply the rules
        def apply_rules(row):
            hap1 = row[hap1_col]
            hap2 = row[hap2_col]
            if hap1 == './.' and hap2 == './.':
                return '0/0'
            elif hap1 == '1/1' and hap2 == './.':
                return '0/1'
            elif hap1 == './.' and hap2 == '1/1':
                return '0/1'
            elif hap1 == '1/1' and hap2 == '1/1':
                return '1/1'
            else:
                return './.'
        df[dip_col] = df.apply(apply_rules, axis=1)
        df = df.drop(columns=[hap1_col, hap2_col])
        final_header_lines = header_lines[:-1]
        final_header = header_lines[-1].strip().split('\t')
        final_header[-2:] = [dip_col]
        final_header_lines.append('\t'.join(final_header) + '\n')
        with open(output_vcf, 'w') as file:
            for line in final_header_lines:
                file.write(line)
            df.to_csv(file, sep='\t', index=False, header=False)

rule bgzip_and_index_vcf_unphased:
    input:
        "vcfs/{ref}/{diploid}/{sample}_unphased.vcf"
    output:
        "vcfs/{ref}/{diploid}/{sample}_unphased.vcf.gz",
        "vcfs/{ref}/{diploid}/{sample}_unphased.vcf.gz.tbi"
    shell: """
    bgzip -c {input} > {output[0]} &&
    bcftools index -t {output[0]}
    """

rule merge_individual_chimps_vcfs:
    input:
        vcf_files=expand('vcfs/{ref}/diploid/{sample}_unphased.vcf.gz', ref=['mPanTro3'], sample=chimpsamples),
    output:
        temp('vcfs/{ref}/diploid/pantros_mapped2{ref}_unphased.vcf.gz'),
        'vcfs/{ref}/diploid/pantros_mapped2{ref}.sorted_unphased.vcf.gz',
        'vcfs/{ref}/diploid/pantros_mapped2{ref}.sorted.filtered_unphased.vcf.gz',
    shell: """
    bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} &&
    bcftools index {output[1]}  &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]}  --threads 56 &&
    bcftools index {output[2]} 
    """

rule merge_individual_bonobos_vcfs:
    input:
        vcf_files=expand('vcfs/{ref}/diploid/{sample}_unphased.vcf.gz', ref=['mPanPan1'], sample=bonobosamples),
    output:
        temp('vcfs/{ref}/diploid/panpa_mapped2{ref}_unphased.vcf.gz'),
        'vcfs/{ref}/diploid/panpa_mapped2{ref}.sorted_unphased.vcf.gz',
        'vcfs/{ref}/diploid/panpa_mapped2{ref}.sorted.filtered_unphased.vcf.gz',
    shell: """
    bcftools merge  -m none --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} &&
    bcftools index {output[1]}  &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]}  --threads 56 &&
    bcftools index {output[2]} 
    """

rule merge_individual_HPRC_vcfs:
    input:
        vcf_files=expand('vcfs/{ref}/diploid/{sample}_unphased.vcf.gz', ref=['ht2t'], sample=human_samples), 
    output:
        temp('vcfs/{ref}/diploid/hprc_mapped2{ref}_unphased.vcf.gz'),
        'vcfs/{ref}/diploid/hprc_mapped2{ref}.sorted_unphased.vcf.gz',
        'vcfs/{ref}/diploid/hprc_mapped2{ref}.sorted.filtered_unphased.vcf.gz',
    shell: """
    bcftools merge  -m none --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
    bcftools sort {output[0]}  --output-type z -o {output[1]}  &&
    bcftools index {output[1]}  &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]} --threads 56  &&
    bcftools index {output[2]}
    """

rule pcas_intraspecies:
    input:
        'vcfs/{ref}/diploid/{dataset}_mapped2{ref}.sorted.filtered_unphased.vcf.gz'
    output: 
        'plink/diploid/{dataset}_mapped2{ref}_unphased.eigenvec',
        'plink/diploid/{dataset}_mapped2{ref}_unphased.eigenval'
    params: 'plink/diploid/{dataset}_mapped2{ref}_unphased'
    shell: """
    plink --vcf {input} --pca --maf 0.05 --double-id --allow-extra-chr --out {params}
    """