import os
import pandas as pd
import numpy as np
import itertools
import gzip

NAMES_FILE = "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/human579_assemblies/names.txt"
sample_table_hgsvc = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_names.txt', names=['names'])
sample_table_hgsvc['base_sample'] = sample_table_hgsvc['names'].apply(lambda x: x.split('.')[0])

# Prevent duplicate samples from crashing bcftools downstream
human_samples_hgsvc = sample_table_hgsvc['base_sample'].unique().tolist()

# ==============================================================================
# 1. Read the raw haploid assemblies
# ==============================================================================
human579_assemblies = []
if os.path.exists(NAMES_FILE):
    with open(NAMES_FILE) as f:
        human579_assemblies = [l.strip() for l in f if l.strip()]
else:
    print(f"Warning: {NAMES_FILE} not found.")

# ==============================================================================
# 2. Group the haploids by their true base sample name (ignoring leading numbers)
# ==============================================================================
human579_dict = {}
for asm in human579_assemblies:
    if '.' in asm:
        base_with_num, hap_info = asm.split('.', 1) 
        
        if '_' in base_with_num:
            prefix, true_base = base_with_num.split('_', 1) 
            if not prefix.isdigit():
                true_base = base_with_num
        else:
            true_base = base_with_num 

        if true_base not in human579_dict:
            human579_dict[true_base] = []
        
        human579_dict[true_base].append(asm)

# ==============================================================================
# 3. Filter for valid diploids (must have exactly 2 haplotypes)
# ==============================================================================
human579_base_samples = []
for base, exact_asms in human579_dict.items():
    if len(exact_asms) == 2:
        human579_base_samples.append(base) 
    else:
        print(f"⚠️ Warning: Skipping {base} because it has {len(exact_asms)} haplotypes instead of 2: {exact_asms}")

print(f"Successfully paired {len(human579_base_samples)} diploid human579 individuals.")

# ==============================================================================
# 4. Remove redundancy and build the Master File Lists
# ==============================================================================
human579_set = set(human579_base_samples)
unique_hgsvc_samples = [sample for sample in human_samples_hgsvc if sample not in human579_set]
print(f"Merging {len(human579_base_samples)} Human579 samples and {len(unique_hgsvc_samples)} unique HGSVC samples.")

# Generate paths pointing to the NEW diploid_phased folder for Human 579
human579_vcf_paths = expand("vcfs/ht2t/diploid_phased/{sample}.vcf.gz", sample=human579_base_samples)

# Generate paths pointing to the existing hgsvc folder for unique HGSVC
hgsvc_vcf_paths = expand("vcfs/ht2t/diploid/hgsvc/{sample}.vcf.gz", sample=unique_hgsvc_samples)

# Combine them into one master list
master_vcf_list = human579_vcf_paths + hgsvc_vcf_paths
master_tbi_list = [vcf + ".tbi" for vcf in master_vcf_list] 

# ==============================================================================
# DAG Target Rule
# ==============================================================================
rule all:
    input:
        #expand("vcfs/ht2t/diploid_phased/{sample}.vcf", sample=human579_base_samples),
        'vcfs/ht2t/diploid_phased/human579_mapped2ht2t.sorted.filtered.vcf.gz',
        'vcfs/ht2t/diploid_phased/human579_hgsvc_master_mapped2ht2t.sorted.filtered.vcf.gz'
        #'plink_diploid/human579_hgsvc_master_mapped2ht2t.eigenval', 
        #'plink_diploid/human579_hgsvc_master_mapped2ht2t.eigenvec'

# ==============================================================================
# Input function to dynamically grab the correct two haploids for each base sample
# ==============================================================================
def get_human579_haploids(wildcards):
    base = wildcards.sample
    exact_asms = human579_dict[base] 
    
    hap1_vcf = f"vcfs/{wildcards.ref}/{exact_asms[0]}.vcf.gz"
    hap2_vcf = f"vcfs/{wildcards.ref}/{exact_asms[1]}.vcf.gz"
    
    return [hap1_vcf, hap2_vcf]

# ==============================================================================
# STEP A: Merge two haploid files into one unphased, 2-column diploid VCF
# ==============================================================================
rule merge_human579_haploids:
    input:
        get_human579_haploids
    output:
        temp("vcfs/{ref}/diploid/{sample}.vcf.gz")
    wildcard_constraints:
        sample="[^/.]+"  # <-- Added the dot here!
    shell:
        """
        bcftools merge -m none -Oz -o {output} {input[0]} {input[1]} --threads 56
        """

# ==============================================================================
# STEP B: Convert the 2-column haploid VCF into a single 1-column phased diploid VCF
# ==============================================================================
rule process_human579_diploid_phased:
    input:
        "vcfs/{ref}/diploid/{sample}.vcf.gz" 
    output:
        temp("vcfs/{ref}/diploid_phased/{sample}.vcf") 
    wildcard_constraints:
        sample="[^/.]+"  # <-- Added the dot here!
    run:
        input_vcf = input[0]
        output_vcf = output[0]
        dip_col = wildcards.sample
        
        hap1_col = human579_dict[dip_col][0] 
        hap2_col = human579_dict[dip_col][1]

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
# STEP C: Zip and Index the new phased diploid files
# ==============================================================================
rule bgzip_and_index_diploid:
    input:
        "vcfs/{ref}/diploid_phased/{sample}.vcf"
    output:
        vcf = temp("vcfs/{ref}/diploid_phased/{sample}.vcf.gz"),
        tbi = temp("vcfs/{ref}/diploid_phased/{sample}.vcf.gz.tbi")
    wildcard_constraints:
        sample="[^/.]+"  # <-- Added the dot here!
    shell: 
        """
        bgzip -c {input} > {output.vcf} &&
        bcftools index -t {output.vcf}
        """


# ==============================================================================
# Merge ONLY Human579 individuals into a cohort VCF
# ==============================================================================
rule merge_only_human579_vcfs:
    input:
        vcf_files = human579_vcf_paths,
        tbi_files = [vcf + ".tbi" for vcf in human579_vcf_paths]
    output:
        vcf = temp('vcfs/{ref}/diploid_phased/human579_mapped2{ref}.vcf.gz'),
        sorted_vcf = temp('vcfs/{ref}/diploid_phased/human579_mapped2{ref}.sorted.vcf.gz'),
        filtered_vcf = 'vcfs/{ref}/diploid_phased/human579_mapped2{ref}.sorted.filtered.vcf.gz'
    shell: 
        """
        bcftools merge -m none -Oz -o {output.vcf} {input.vcf_files} --threads 56 &&
        bcftools sort -m 50G -T /global/scratch/users/joana_rocha/tmp {output.vcf} --output-type z -o {output.sorted_vcf} &&
        bcftools index {output.sorted_vcf} &&
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output.sorted_vcf} -o {output.filtered_vcf} --threads 56 &&
        bcftools index {output.filtered_vcf}
        """

# ==============================================================================
# STEP D: Merge ALL individuals (Human579 phased + HGSVC unique) into a master cohort
# ==============================================================================
rule merge_master_cohort:
    input:
        vcf_files = master_vcf_list,
        tbi_files = master_tbi_list 
    output:
        vcf = temp("vcfs/ht2t/diploid_phased/human579_hgsvc_master_mapped2ht2t.vcf.gz"),
        sorted_vcf = temp("vcfs/ht2t/diploid_phased/human579_hgsvc_master_mapped2ht2t.sorted.vcf.gz"),
        filtered_vcf = "vcfs/ht2t/diploid_phased/human579_hgsvc_master_mapped2ht2t.sorted.filtered.vcf.gz"
    shell:
        """
        bcftools merge -m none -Oz -o {output.vcf} {input.vcf_files} --threads 56 &&
        bcftools sort -m 50G -T /global/scratch/users/joana_rocha/tmp {output.vcf} --output-type z -o {output.sorted_vcf} &&
        bcftools index {output.sorted_vcf} &&
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output.sorted_vcf} -o {output.filtered_vcf} --threads 56 &&
        bcftools index {output.filtered_vcf}
        """

# ==============================================================================
# Run PCA on the master diploid cohort (Optional)
# ==============================================================================
rule pcas_master_diploid:
    input:
        'vcfs/{ref}/diploid/human579_hgsvc_master_merged.vcf.gz'
    output: 
        'plink_diploid/human579_hgsvc_master_mapped2{ref}.eigenvec',
        'plink_diploid/human579_hgsvc_master_mapped2{ref}.eigenval'
    params: 
        prefix='plink_diploid/human579_hgsvc_master_mapped2{ref}'
    shell: 
        """
        plink --vcf {input} --pca --double-id --allow-extra-chr --out {params.prefix}
        """