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
sample_table_hgsvc = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_names.txt', names=['names'])
sample_table_hgsvc['base_sample'] = sample_table_hgsvc['names'].apply(lambda x: x.split('.')[0])

samples = sample_table['Specimen'].unique()
chimpsamples = chimp_sample_table['Specimen'].unique()
bonobosamples = bonobo_sample_table['Specimen'].unique()

#input vcfs are those from minimap2snvs.py

human_samples_hprc = sample_table_hprc['base_sample'].unique()
human_samples_hprc_list = human_samples_hprc.tolist()
human_samples_hgsvc = sample_table_hgsvc['base_sample'].tolist()

all_human_samples = list(set(human_samples_hprc_list + human_samples_hgsvc))
#print(all_human_samples)

haps = ["hap1", "hap2"]
sample_haps = expand("{sample}.{hap}", sample=samples, hap=haps)
human_sample_haps =sample_table_hprc['names'].unique()

all_samples = list(all_human_samples) + list(chimpsamples) + list(bonobosamples)

#print(sample_haps)
# print(human_sample_haps)
print(all_samples)

rule all:
    input:
        expand("vcfs/ht2t/diploid/hgsvc/{sample}.vcf.gz", sample=human_samples_hgsvc),
        #expand("vcfs/mPonAbe1/diploid/{sample}.vcf.gz", sample=human_samples_hprc),
        #expand("vcfs/mPonAbe1/diploid/{sample}.vcf.gz", sample=human_samples_hgsvc),
        #'vcfs/ht2t/diploid/hprc_hgsvc_mapped2ht2t.sorted.vcf.gz',
        #expand("vcfs/mPanTro3/diploid/{sample}.vcf.gz", sample=chimpsamples),
        #expand("vcfs/mPonAbe1/diploid/{sample}.vcf.gz", sample=chimpsamples),
         #expand("vcfs/mPonAbe1/diploid/{sample}.vcf.gz", sample=chimpsamples),
        #'vcfs/mPanTro3/diploid/pantros_mapped2mPanTro3.sorted.vcf.gz',
        #'plink_diploid/pantros_mapped2mPanTro3.eigenval', #pantros (no bonobos)
        #'plink_diploid/pantros_mapped2mPanTro3.eigenvec', #pantros (no bonobos)
        #expand("vcfs/mPanPan1/diploid/{sample}.vcf.gz", sample=bonobosamples),
        #expand("vcfs/mPonAbe1/diploid/{sample}.vcf.gz", sample=bonobosamples),
        #expand("vcfs/mPonAbe1/diploid/{sample}.vcf.gz", sample=bonobosamples),
        #'vcfs/mPonAbe1/diploid/allsamples_mapped2mPonAbe1.sorted.vcf.gz',
        #'vcfs/mPanPan1/diploid/panpa_mapped2mPanPan1.sorted.vcf.gz',
        #'vcfs/mPonAbe1/diploid/allsamples_mapped2mPonAbe1.sorted.filtered.vcf.gz',
        #'plink_diploid/allsamples_mapped2mPonAbe1.eigenval', 
        #'plink_diploid/allsamples_mapped2mPonAbe1.eigenvec', 


# ==============================================================================
# STEP 1: Merge HGSVC haploids into a temporary 2-column file
# ==============================================================================
rule merge_hgsvc_temp:
    input:
        hap1="vcfs/{ref}/{sample}.1.vcf.gz",
        hap2="vcfs/{ref}/{sample}.2.vcf.gz"
    output:
        temp("vcfs/{ref}/diploid/hgsvc/{sample}.merged.vcf.gz") # Renamed to avoid loop
    wildcard_constraints:
        sample="[^/.]+"
    shell:
        """
        bcftools merge -m none -Oz -o {output} {input.hap1} {input.hap2} --threads 56
        """

# ==============================================================================
# STEP 2: Convert the 2-column HGSVC into a 1-column phased diploid
# ==============================================================================
rule process_hgsvc_diploid_phased:
    input:
        "vcfs/{ref}/diploid/hgsvc/{sample}.merged.vcf.gz"
    output:
        temp("vcfs/{ref}/diploid/hgsvc/{sample}.vcf")
    wildcard_constraints:
        sample="[^/.]+"
    run:
        import gzip
        
        input_vcf = input[0]
        output_vcf = output[0]
        
        # Hardcoding the .1 and .2 suffixes based on your input files
        hap1_col = wildcards.sample + ".1"
        hap2_col = wildcards.sample + ".2"
        dip_col = wildcards.sample

        with gzip.open(input_vcf, 'rt') as fin, open(output_vcf, 'w') as fout:
            idx1, idx2 = None, None
            
            for line in fin:
                if line.startswith('##'):
                    fout.write(line)
                elif line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    idx1 = header.index(hap1_col)
                    idx2 = header.index(hap2_col)
                    
                    new_header = [h for i, h in enumerate(header) if i not in (idx1, idx2)]
                    new_header.append(dip_col)
                    fout.write('\t'.join(new_header) + '\n')
                else:
                    cols = line.strip().split('\t')
                    hap1 = cols[idx1]
                    hap2 = cols[idx2]

                    if hap1 == './.' and hap2 == './.':
                        dip_gt = '0|0'
                    elif hap1 == '1/1' and hap2 == './.':
                        dip_gt = '1|0'
                    elif hap1 == './.' and hap2 == '1/1':
                        dip_gt = '0|1'
                    elif hap1 == '1/1' and hap2 == '1/1':
                        dip_gt = '1|1'
                    else:
                        dip_gt = './.'

                    new_cols = [c for i, c in enumerate(cols) if i not in (idx1, idx2)]
                    new_cols.append(dip_gt)
                    fout.write('\t'.join(new_cols) + '\n')

# ==============================================================================
# STEP 3: Zip and Index the final phased HGSVC files
# ==============================================================================
rule bgzip_and_index_hgsvc:
    input:
        "vcfs/{ref}/diploid/hgsvc/{sample}.vcf"
    output:
        vcf="vcfs/{ref}/diploid/hgsvc/{sample}.vcf.gz",
        tbi="vcfs/{ref}/diploid/hgsvc/{sample}.vcf.gz.tbi"
    wildcard_constraints:
        sample="[^/.]+"
    shell: 
        """
        bgzip -c {input} > {output.vcf} &&
        bcftools index -t {output.vcf}
        """
        
rule merge_vcf:
    input:
        hap1="vcfs/{ref}/{sample}.hap1.vcf.gz",
        hap2="vcfs/{ref}/{sample}.hap2.vcf.gz"
    output:
        "vcfs/{ref}/temp_{sample}.vcf.gz"
    shell:
        "bcftools merge -m none -Oz -o {output} {input.hap1} {input.hap2} --threads 56"



#rule process_panpan_vcf_diploid_phased:
#    input:
#        "vcfs/mPonAbe1/temp_{sample}.vcf.gz"
#    output:
#        "vcfs/mPonAbe1/diploid/{sample}.vcf"
#    run:
#        input_vcf = input[0]
#        output_vcf = output[0]
#        hap1_col = wildcards.sample + ".hap1"
#        hap2_col = wildcards.sample + ".hap2"
#        dip_col = wildcards.sample

        # Read the input VCF file
#        with gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r') as file:
#            lines = file.readlines()#

#        # Extract header and data lines
#        header_lines = [line for line in lines if line.startswith('#')]
#        data_lines = [line for line in lines if not line.startswith('#')]

        # Create a dataframe from the data lines
#        data = [line.strip().split('\t') for line in data_lines]
#        columns = header_lines[-1].strip().split('\t')
#        df = pd.DataFrame(data, columns=columns)

        # Function to apply the rules
#        def apply_rules(row):
#            hap1 = row[hap1_col]
#            hap2 = row[hap2_col]
#            if hap1 == './.' and hap2 == './.':
#                return '0|0'
#            elif hap1 == '1/1' and hap2 == './.':
#                return '1|0'
#            elif hap1 == './.' and hap2 == '1/1':
#                return '0|1'
#            elif hap1 == '1/1' and hap2 == '1/1':
#                return '1|1'
#            else:
#                return './.'
#        df[dip_col] = df.apply(apply_rules, axis=1)
#        df = df.drop(columns=[hap1_col, hap2_col])
#        final_header_lines = header_lines[:-1]
#        final_header = header_lines[-1].strip().split('\t')
#        final_header[-2:] = [dip_col]
#        final_header_lines.append('\t'.join(final_header) + '\n')
#        with open(output_vcf, 'w') as file:
#            for line in final_header_lines:
#                file.write(line)
#            df.to_csv(file, sep='\t', index=False, header=False)


rule bgzip_and_index_vcf:
    input:
        "vcfs/{ref}/{diploid}/{sample}.vcf"
    output:
        "vcfs/{ref}/{diploid}/{hgsvc}/{sample}.vcf.gz",
        "vcfs/{ref}/{diploid}/{hgsvc}/{sample}.vcf.gz.tbi"
    shell: """
    bgzip -c {input} > {output[0]} &&
    bcftools index -t {output[0]}
    """

rule merge_individual_chimps_vcfs:
    input:
        vcf_files=expand('vcfs/{ref}/diploid/{sample}.vcf.gz', ref=['mPanTro3'], sample=chimpsamples),
    output:
        temp('vcfs/{ref}/diploid/pantros_mapped2{ref}.vcf.gz'),
        'vcfs/{ref}/diploid/pantros_mapped2{ref}.sorted.vcf.gz',
        'vcfs/{ref}/diploid/pantros_mapped2{ref}.sorted.filtered.vcf.gz',
    shell: """
    bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} &&
    bcftools index {output[1]}  &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]}  --threads 56 &&
    bcftools index {output[2]} 
    """

rule merge_individual_bonobos_vcfs:
    input:
        vcf_files=expand('vcfs/{ref}/diploid/{sample}.vcf.gz', ref=['mPanPan1'], sample=bonobosamples),
    output:
        temp('vcfs/{ref}/diploid/panpa_mapped2{ref}.vcf.gz'),
        'vcfs/{ref}/diploid/panpa_mapped2{ref}.sorted.vcf.gz',
        'vcfs/{ref}/diploid/panpa_mapped2{ref}.sorted.filtered.vcf.gz',
    shell: """
    bcftools merge  -m none --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} &&
    bcftools index {output[1]}  &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]}  --threads 56 &&
    bcftools index {output[2]} 
    """

rule merge_individual_humans_vcfs:
    input:
        vcf_files=expand('vcfs/{ref}/diploid/{sample}.vcf.gz', ref=['ht2t'], sample=all_human_samples),
    output:
        temp('vcfs/{ref}/diploid/hprc_hgsvc_mapped2{ref}.vcf.gz'),
        'vcfs/{ref}/diploid/hprc_hgsvc_mapped2{ref}.sorted.vcf.gz',
        'vcfs/{ref}/diploid/hprc_hgsvc_mapped2{ref}.sorted.filtered.vcf.gz',
    shell: """
    bcftools merge  -m none --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
    bcftools sort {output[0]}  --output-type z -o {output[1]}  &&
    bcftools index {output[1]}  &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]} --threads 56  &&
    bcftools index {output[2]}
    """


#rule merge_individual_samples_vcfs:
#    input:
#        vcf_files=expand('vcfs/{ref}/diploid/{sample}.vcf.gz', ref=['mPonAbe1'], sample=all_samples),
#    output:
#        temp('vcfs/{ref}/diploid/allsamples_mapped2{ref}.vcf.gz'),
#        'vcfs/{ref}/diploid/allsamples_mapped2{ref}.sorted.vcf.gz',
#        'vcfs/{ref}/diploid/allsamples_mapped2{ref}.sorted.filtered.vcf.gz',
#    shell: """
#    bcftools merge  -m none --output-type z -o {output[0]} {input.vcf_files} --threads 56 &&
#    bcftools sort {output[0]}  --output-type z -o {output[1]}  &&
#    bcftools index {output[1]}  &&
#    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]} --threads 56  &&
#    bcftools index {output[2]}
#    """

rule pcas_intraspecies_phased:
    input:
        'vcfs/{ref}/diploid/{dataset}_mapped2{ref}.sorted.filtered.vcf.gz'
    output: 
        'plink_diploid/{dataset}_mapped2{ref}.eigenvec',
        'plink_diploid/{dataset}_mapped2{ref}.eigenval'
    params: 'plink_diploid/{dataset}_mapped2{ref}'
    shell: """
    plink --vcf {input} --pca --double-id --allow-extra-chr --out {params}
    """