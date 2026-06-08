import numpy as np
import pandas as pd
import os
import itertools

AGC_FILE = "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/human579_plus/human_579_plus_HG02554.1.2.agc" 
NAMES_FILE = "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/human579_assemblies/names.txt" # Generated via: agc listset human_579_plus_HG02554.1.2.agc > names.txt

# Read the list of assemblies directly from the file
human579_assemblies = []
if os.path.exists(NAMES_FILE):
    with open(NAMES_FILE) as f:
        human579_assemblies = [l.strip() for l in f if l.strip()]
else:
    print(f"Warning: {NAMES_FILE} not found. Human579 specific rules may typically fail if this is missing.")


ref_paths = {
    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    'mPanTro3' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta", ### no difference to thanksgiving/final primary
    'mPonAbe1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta",
    'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
}

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
chimp_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/chimpsamples.tsv', sep='\t')
bonobo_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/bonobosamples.tsv', sep='\t')
sample_table_hprc=pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_names.txt', names=['names']) ### removed Y chromosome
sample_table_hgsvc = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_names.txt', names=['names'])

references = ['ht2t', 'mPanTro3', 'mPanPan1', 'mPonAbe1']

samples = sample_table['Specimen'].unique()

chimpsamples = chimp_table['Specimen'].unique()
print(chimpsamples)

bonobosamples = bonobo_table['Specimen'].unique()
print(bonobosamples)

haps = ["hap1", "hap2"]
sample_haps = expand("{sample}.{hap}", sample=samples, hap=haps)

chimp_sample_haps = expand("{sample}.{hap}", sample=chimpsamples, hap=haps)
bonobo_sample_haps = expand("{sample}.{hap}", sample=bonobosamples, hap=haps)

human_samples_hprc =sample_table_hprc['names'].unique()
human_samples_hgsvc = sample_table_hgsvc['names'].tolist()

# Combine all human samples: HGSVC + HPRC + Human579 (from AGC) with Dantu Hybrid reassembled 
#all_human_samples = list(human_samples_hgsvc) + list(human_samples_hprc) + list(human579_assemblies)
#human_base_sample_names = set(sample.split('.')[0] for sample in all_human_samples)
#print(all_human_samples)
#print(human579_assemblies)
print(human_samples_hgsvc)
#ruleorder: agc_map > hgsvc2ref

rule all:
    input:
        #expand('vcfs/{anyref}/{sample}.{hap}.vcf.gz', sample=bonobosamples, hap=haps, anyref=['mPanPan1']),
        #expand('vcfs/{anyref}/{sample}.{hap}.vcf.gz', sample=chimpsamples, hap=haps, anyref=['mPanTro3']),     
        # This now includes HPRC, HGSVC, AND Human579 individual VCFs
        #expand('vcfs/{anyref}/{anything}.srt.paf', anything=human579_assemblies, anyref=['ht2t']),
        #expand('vcfs/{anyref}/{anything}.vcf.gz', anything=human579_assemblies, anyref=['ht2t']),
        #expand('vcfs/{anyref}/{anything}.vcf.gz', anything=human_samples_hgsvc, anyref=['ht2t']),
        #expand('vcfs/{anyref}/{sample}.{hap}.vcf.gz', sample=bonobosamples, hap=haps, anyref=['mPonAbe1']),
        #expand('vcfs/{anyref}/{sample}.{hap}.vcf.gz', sample=chimpsamples, hap=haps, anyref=['mPonAbe1']),
        #expand('vcfs/{anyref}/{anything}.vcf.gz', anything=human_samples_hprc, anyref=['mPonAbe1']),
        #expand('vcfs/{anyref}/{anything}.vcf.gz', anything=human_samples_hgsvc, anyref=['mPonAbe1']),
        #'vcfs/mPanTro3/pantros_mapped2mPanTro3.sorted.filtered.vcf.gz',
        #'vcfs/mPanPan1/panpa_mapped2mPanPan1.sorted.filtered.vcf.gz',
        'plink/pantros_mapped2mPanTro3.eigenvec',
        'plink/pantros_mapped2mPanTro3.eigenval',
        #'plink/panpa_mapped2mmPanPan1.eigenvec',
        #'plink/panpa_mapped2mmPanPan1.eigenval'

def get_fasta_path(wildcards):
    sample = wildcards.sample
    hap = wildcards.hap
    
    #verkko_path = f"/global/scratch/users/joana_rocha/PANPAN/Verkko-fasta_shortcut/{sample}.{hap}.verkko.fasta"
    hic_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta-HiC_shortcut/{sample}.{hap}.p_ctg.hic.fa"
    hifi_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{sample}.{hap}.p_ctg.fa"
      
    # Check if verkko file exists
    #if os.path.exists(verkko_path):
    #    return verkko_path
    # If verkko file does not exist, check if hic file exists
    if os.path.exists(hic_path):
        return hic_path
    # If neither verkko nor hic files exist, return the non-hic path
    else:
        return hifi_path


rule agc_map:
    output: 
        paf = 'pafs/{anyref}/{assembly_name}.paf',
        srt_paf = 'pafs/{anyref}/{assembly_name}.srt.paf',
    params:
        ref_path = lambda wildcards: ref_paths[wildcards.anyref],
        agc_file = AGC_FILE,
        # The wildcard {assembly_name} IS the internal AGC name (e.g. 120001_CN1.pat)
        asm_name = lambda wildcards: wildcards.assembly_name
    threads: 24
    wildcard_constraints:
        # Only allow this rule if the assembly name matches our list
        assembly_name = "|".join(human579_assemblies) if human579_assemblies else "NOMATCH"
    run:
        shell("""
        agc getset {params.agc_file} {params.asm_name} | \
        minimap2 -cx asm5 -t 20 --cs {params.ref_path} - -t {threads} > {output.paf} && \
        sort -k6,6 -k8,8n {output.paf} > {output.srt_paf}
        """)


#rule pan2ref:
#    input:
#        get_fasta_path,
#    output: 
#        'pafs/{anyref}/{sample}.{hap}.paf',
#        'pafs/{anyref}/{sample}.{hap}.srt.paf',
#    params:
#        ref_path = lambda wildcards: ref_paths[wildcards.anyref],
#    shell: """
#    minimap2 -cx asm20 --cs {params.ref_path} {input} -t 24  > {output[0]} &&  
#    sort -k6,6 -k8,8n {output[0]} > {output[1]} 
#    """

#rule hprc2ref:
#    input: "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_assemblies/{human_sample_hprc}.fa"
#    output: 
#        "pafs/{anyref}/{human_sample_hprc}.paf",
#        "pafs/{anyref}/{human_sample_hprc}.srt.paf",
#    params:
#        ref_path = lambda wildcards: ref_paths[wildcards.anyref],
#    shell: """
#    minimap2 -cx asm20 --cs {params.ref_path} {input} -t 24  > {output[0]} &&  
#    sort -k6,6 -k8,8n {output[0]} > {output[1]}
#    """

rule hgsvc2ref:
    input:
        "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HGSVC_assemblies/{human_sample_hgsvc}.fasta.gz"
    output: 
        "pafs/{anyref}/{human_sample_hgsvc}.paf",
        "pafs/{anyref}/{human_sample_hgsvc}.srt.paf",
    params:
        ref_path = lambda wildcards: ref_paths[wildcards.anyref],
    shell: """
    minimap2 -cx asm20 --cs {params.ref_path} {input} -t 24  > {output[0]} &&  
    sort -k6,6 -k8,8n {output[0]} > {output[1]}
    """


rule paf2snvs:
    input: 
        'pafs/{anyref}/{anything}.srt.paf',
    output:
        'vcfs/{anyref}/{anything}.vcf',
        'logs/{anyref}/{anything}.log',
    params:
        ref_path = lambda wildcards: ref_paths[wildcards.anyref],
    shell: """
    paftools.js call -f {params.ref_path} -s {wildcards.anything} {input} - > {output[0]} 2> {output[1]}
    """

rule bgzip_vcfs:
    input:
        'vcfs/{anyref}/{anything}.vcf',
    output:
        vcf_gz='vcfs/{anyref}/{anything}.vcf.gz',
        vcf_gz_tbi='vcfs/{anyref}/{anything}.vcf.gz.tbi',
    shell: """
    bgzip -c {input} > {output.vcf_gz} &&
    tabix -p vcf {output.vcf_gz}
    """

rule merge_individual_chimps_vcfs: #bcftools query -l panpan_mapped2clint_filtered.vcf.gz > panpan_names_list.txt
    input:
        vcf_files=expand('vcfs/mPanTro3/{sample}.vcf.gz', sample=chimp_sample_haps),
    output:
        'vcfs/mPanTro3/pantros_mapped2mPanTro3.vcf.gz',
        'vcfs/mPanTro3/pantros_mapped2mPanTro3.sorted.vcf.gz',
        'vcfs/mPanTro3/pantros_mapped2mPanTro3.sorted.filtered.vcf.gz',
    shell: """
    bcftools merge --merge snps --missing-to-ref  --output-type z -o {output[0]} {input.vcf_files} &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} 
    bcftools index {output[1]} &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]}  &&
    bcftools index {output[2]} 
    """
#bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} &&

rule merge_individual_bonobos_vcfs: #bcftools query -l panpan_mapped2clint_filtered.vcf.gz > panpan_names_list.txt
    input:
        vcf_files=expand('vcfs/mPanPan1/{sample}.vcf.gz', sample=bonobo_sample_haps),
    output:
        'vcfs/mPanPan1/panpa_mapped2mPanPan1.vcf.gz',
        'vcfs/mPanPan1/panpa_mapped2mPanPan1.sorted.vcf.gz',
        'vcfs/mPanPan1/panpa_mapped2mPanPan1.sorted.filtered.vcf.gz',
    shell: """
    bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} 
    bcftools index {output[1]} &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[1]}  -o {output[2]}  &&
    bcftools index {output[2]} 
    """

rule merge_individual_HPRC_vcfs:
    input:
        vcf_files=expand('vcfs/ht2t/{sample}.vcf.gz', sample=human_samples_hprc),
    output:
        'vcfs/ht2t/hprc_mapped2ht2t.vcf.gz',
        'vcfs/ht2t/hprc_mapped2ht2t_filtered.vcf.gz',
    shell: """
    bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} &&
    bcftools index {output[0]} &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[0]}  -o {output[1]}  &&
    bcftools index {output[1]} 
    """

# ==============================================================================
# NEW MERGE RULE FOR HUMAN579
# ==============================================================================
rule merge_human579:
    input:
        # We use the raw assembly names (e.g. 120001_CN1.pat) as the identifier
        vcf_files = expand('vcfs/{anyref}/{assembly_name}.vcf.gz', anyref="{anyref}", assembly_name=human579_assemblies),
    output:
        'vcfs/{anyref}/human579_merged.vcf.gz',
        'vcfs/{anyref}/human579_merged_filtered.vcf.gz',
    shell: """
    bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} &&
    bcftools index {output[0]} &&
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[0]} -o {output[1]} &&
    bcftools index {output[1]} 
    """
# ==============================================================================


#rule merge_individual_HPRC_chimps_bonobos_vcfs:
#    input:
#        vcf_files=expand('vcfs/mPonAbe1/{sample}.vcf.gz', sample=all_samples),
#    output:
#        'vcfs/mPonAbe1/panpanhprc_mapped2mPonAbe1.vcf.gz',
#        'vcfs/mPonAbe1/panpanhprc_mapped2mPonAbe1_filtered.vcf.gz',
#    shell: """
#    bcftools merge --merge snps --missing-to-ref --output-type z -o {output[0]} {input.vcf_files} &&
#    bcftools index {output[0]} &&
#    bcftools view --types snps --min-alleles 2 --max-alleles 2 --output-type z {output[0]}  -o {output[1]}  &&
#    bcftools index {output[1]} 
#    """

rule pcas_chimps:
    input:
        'vcfs/mPanTro3/pantros_mapped2mPanTro3.sorted.filtered.vcf.gz'
    output: 
        'plink/pantros_mapped2mPanTro3.eigenvec',
        'plink/pantros_mapped2mPanTro3.eigenval'
    params: 'plink/pantros_mapped2mPanTro3'
    shell: """
    plink --vcf {input} --pca --maf 0.05 --double-id --allow-extra-chr --out {params}
    """

rule pcas_bonobos:
    input:
        'vcfs/mPanPan1/panpa_mapped2mPanPan1.sorted.filtered.vcf.gz'
    output: 
        'plink/panpa_mapped2mmPanPan1.eigenvec',
        'plink/panpa_mapped2mmPanPan1.eigenval'
    params: 'plink/panpa_mapped2mPanPan1'
    shell: """
    plink --vcf {input} --pca --maf 0.05 --double-id --allow-extra-chr --out {params}
    """


rule pcas_hprc:
    input: 
        'vcfs/ht2t/hprc_mapped2ht2t_filtered.vcf.gz'
    output: 
        'plink/hprc_mapped2ht2t.eigenvec',
        'plink/hprc_mapped2ht2t.eigenval',
    params: 'plink/hprc_mapped2ht2t'
    shell: """
    plink --vcf {input} --pca --maf 0.05 --double-id --allow-extra-chr --out {params}
    """

rule pcas_human579:
    input: 
        'vcfs/{anyref}/human579_merged_filtered.vcf.gz'
    output: 
        'plink/human579_mapped2{anyref}.eigenvec',
        'plink/human579_mapped2{anyref}.eigenval',
    params: 
        prefix = 'plink/human579_mapped2{anyref}'
    shell: """
    plink --vcf {input} --pca --maf 0.05 --double-id --allow-extra-chr --out {params.prefix}
    """

#rule vcf_regions: ### only interested in panpan hprc mapped to hg38 and take out the Y sample
#    input:
#        'vcfs/hg38/panpanhprc_mapped2hg38_filtered.vcf.gz'
#    output:
#        'vcfs/hg38/panpanhprc_mapped2hg38_filtered_{region}.vcf.gz'
#    params:
#        region_coordinate = lambda wildcards: region_coordinates[wildcards.region],
#    shell:  'bcftools view -m2 -M2 -v snps {input} -r {params.region_coordinate}  -Ov --output-type z -o {output}' 


#rule vcf2fasta:
#    input: 'vcfs/hg38/panpanhprc_mapped2hg38_filtered.vcf.gz'
#    output: 'trees/panpanhprc_mapped2hg38_filtered.min4.fasta'
#    params: 'panpanhprc_mapped2hg38_filtered'
#    shell: """
#    ./vcf2phylip.py -i {input} --phylip-disable --fasta  --output-folder trees  --output-prefix {params}
#    """

#rule vcf2fasta_region:
#    input: 'vcfs/hg38/panpanhprc_mapped2hg38_filtered_{region}.vcf.gz'
#    output: 'trees/panpanhprc_mapped2hg38_filtered_{region}.min4.fasta'
#    params: 'panpanhprc_mapped2hg38_filtered_{region}'
#    shell: """
#    ./vcf2phylip.py -i {input} --phylip-disable --fasta  --output-folder trees  --output-prefix {params}
#    """

#rule tree_region:
#    input: 'trees/panpanhprc_mapped2hg38_filtered_{region}.min4.fasta'
#    output: 'trees/panpanhprc_mapped2hg38_filtered_{region}_iqtree.log'
#    params: 'trees/panpanhprc_mapped2hg38_filtered_{region}'
#    shell: """
#    iqtree --prefix {params} -s {input} -m TEST -bb 1000 -nt 20  >  {output}
#    """

#rule tree_region_timecalibrated:
#    input: 
#        'trees/panpanhprc_mapped2hg38_filtered_{region}.min4.fasta',
#        'trees/panpanhprc_mapped2hg38_filtered_{region}.min4.treefile'
#    output: 'trees/panpanhprc_mapped2hg38_filtered_{region}_iqtree_timed.log'
#    params: 'trees/panpanhprc_mapped2hg38_filtered_{region}_timed'
#    shell: """
#    iqtree -s {input[0]} -te {input[1]} --prefix {params} -nt 20  --date-root -6000000 > {output}
#    """