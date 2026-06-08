import numpy as np
import pandas as pd
import os
import itertools
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

references = ['ht2t', 'mPanTro3', 'mPanPan1', 'mPonAbe1']

samples = sample_table['Specimen'].unique()
haps = ["hap1", "hap2"]


chimpsamples = chimp_table['Specimen'].unique()
print(chimpsamples)

bonobosamples = bonobo_table['Specimen'].unique()
print(bonobosamples)

human_samples =sample_table_hprc['names'].unique()
human_base_sample_names = set(sample.split('.')[0] for sample in human_samples)
print(human_base_sample_names)

all_samples= list(samples) + list(human_base_sample_names)
#all_samples= sample_haps + list(human_samples)
print(all_samples)

rule all:
    input:
        #panapn and hprc to mPonAbe1
        #expand('asm/{ref}/{anything}.srt.rg.bam', ref=['mPonAbe1'], anything=samples),  
        #expand('asm/{ref}/{anything}.srt.rg.bam', ref=['mPonAbe1'], anything=human_samples),
        expand('svs/{ref}/{anything}/{anything}.sorted.vcf.gz', ref=['mPonAbe1'], anything=samples),
        expand('svs/{ref}/{anything}/{anything}.sorted.vcf.gz', ref=['mPonAbe1'], anything=human_base_sample_names),
        #'svs/mPanPan1/panpan_hprc_vcf_files_raw_calls.txt',
        #'svs/svasm_truvari/mPanPan1/panpan_hprc.truvari.merged.vcf.gz',
        # hprc to ht2t
        #expand('asm/{ref}/{anything}.srt.rg.bam', ref=['ht2t'], anything=human_samples), 
        #expand('svs/{ref}/{anything}/{anything}.sorted.vcf', ref=['ht2t'], anything=human_base_sample_names),
        #'svs/ht2t/hprc_vcf_files_raw_calls.txt',
        #'svs/svasm_truvari/ht2t/hprc.truvari.merged.vcf.gz',
        #'svs/survivor/ht2t/hprc.merged.vcf',
        #pt to mPanTro3,  
        #expand('asm/{ref}/{anything}.srt.rg.bam', ref=['mPanTro3'], anything=chimpsamples),
        expand('svs/{ref}/{anything}/{anything}.sorted.vcf.gz', ref=['mPanTro3'], anything=chimpsamples),
        #'svs/mPanTro3/panpan-pt_vcf_files_raw_calls.txt',
        #'svs/svasm_truvari/mPanTro3/panpan-pt.truvari.merged.vcf.gz',
        # ppa to mPanPan1
        #expand('asm/{ref}/{anything}.srt.rg.bam', ref=['mPanPan1'], anything=bonobosamples),
        expand('svs/{ref}/{anything}/{anything}.sorted.vcf.gz', ref=['mPanPan1'], anything=bonobosamples),
        #'svs/mPanPan1/panpan-ppa_vcf_files_raw_calls.txt',
        #'svs/svasm_truvari/mPanPan1/panpan-ppa.truvari.merged.vcf.gz',
        #hprc to hg38
        #expand('asm/{ref}/{anything}.srt.rg.bam', ref=['hg38_noALT'], anything=human_samples), 
        #expand('svs/{ref}/{anything}/{anything}.sorted.vcf', ref=['hg38_noALT'], anything=human_base_sample_names),
        #expand('svs/{ref}/{anything}/{anything}.sorted.vcf.gz', ref=['hg38_noALT'], anything=human_base_sample_names),
        #'svs/hg38_noALT/hprc_vcf_files_raw_calls.txt',
        #'svs/svasm_truvari/hg38_noALT/hprc.truvari.merged.vcf.gz',
        #'svs/svasm_survivor/hg38_noALT/hprc.merged.vcf',
        #'pav/run.success',

   
def get_fasta_path(wildcards):
    verkko_path = f"/global/scratch/users/joana_rocha/PANPAN/Verkko-fasta_shortcut/{wildcards.sample.replace('.hic', '')}.verkko.fasta"
    hic_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta-HiC_shortcut/{wildcards.sample.replace('.hic', '')}.p_ctg.hic.fa"
    hifi_path = f"/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{wildcards.sample.replace('.hic', '')}.p_ctg.fa"
    
    # Check if verkko file exists
    #if os.path.exists(verkko_path):
    #    return verkko_path
    # If verkko file does not exist, check if hic file exists
    if os.path.exists(hic_path):
        return hic_path
    # If neither verkko nor hic files exist, return the non-hic path
    else:
        return hifi_path  
        
rule panpan2ref:
    input:
        get_fasta_path
    output: 
        temp('asm/{ref}/{sample}.sam'),
        temp('asm/{ref}/{sample}.srt.bam'),
        'asm/{ref}/{sample}.srt.rg.bam',
    params:
        ref = lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    minimap2 -a -x asm5 --cs -r2k -t 20 {params.ref} {input} > {output[0]} &&  
    samtools sort -m4G -@4 -o {output[1]} {output[0]} &&
    samtools addreplacerg -r ID:{wildcards.sample} -r PL:PACBIO -o {output[2]} {output[1]} 
    samtools index {output[2]}
    """


rule hprc2ref:
    input:
        "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_assemblies/{human_sample}.fa",
    output: 
        temp('asm/{ref}/{human_sample}.sam'),
        temp('asm/{ref}/{human_sample}.srt.bam'),
        'asm/{ref}/{human_sample}.srt.rg.bam',
    params:
        ref = lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    minimap2 -a -x asm5 --cs -r2k -t 20 {params.ref} {input[0]} > {output[0]} &&  
    samtools sort -m4G -@4 -o {output[1]} {output[0]} &&
    samtools addreplacerg -r ID:{wildcards.human_sample} -r PL:PACBIO -o {output[2]} {output[1]} 
    samtools index {output[2]}
    """

#### run svim asm 
#rule asm2svs:
#    input: 
#        'asm/{ref}/{anything}.srt.rg.bam',
#    output: 
#       temp('svs/{ref}/{anything}/{anything}.vcf'),
#    params:
#        ref = lambda wildcards: ref_paths[wildcards.ref],
#        anything = lambda wildcards: wildcards.anything,
#    shell: """
#    svim-asm haploid 'svs/{wildcards.ref}/{wildcards.anything}' {input} {params.ref} --min_sv_size 40 --query_names --interspersed_duplications_as_insertions &&
#    mv 'svs/{wildcards.ref}/{wildcards.anything}/variants.vcf' 'svs/{wildcards.ref}/{wildcards.anything}/{wildcards.anything}.vcf' &&
#    sed -i 's|Sample|{wildcards.anything}|g' 'svs/{wildcards.ref}/{wildcards.anything}/{wildcards.anything}.vcf' 
#    """

def construct_sample_pairs(sample_names, patterns):
    return {base_name: [pattern.format(ref='{ref}', base_name=base_name) for pattern in patterns] 
            for base_name in sample_names}

human_sample_pairs = construct_sample_pairs(human_base_sample_names, [
    'asm/{ref}/{base_name}.1.srt.rg.bam',
    'asm/{ref}/{base_name}.2.srt.rg.bam'
])

chimp_sample_pairs = construct_sample_pairs(chimpsamples, [
    'asm/{ref}/{base_name}.hap1.srt.rg.bam',
    'asm/{ref}/{base_name}.hap2.srt.rg.bam'
])

bonobo_sample_pairs = construct_sample_pairs(bonobosamples, [
    'asm/{ref}/{base_name}.hap1.srt.rg.bam',
    'asm/{ref}/{base_name}.hap2.srt.rg.bam'
])

sample_pairs = {**human_sample_pairs, **chimp_sample_pairs, **bonobo_sample_pairs}

rule asm2svs_diploid:
    input: 
        hap1=lambda wildcards: sample_pairs[wildcards.anything][0],
        hap2=lambda wildcards: sample_pairs[wildcards.anything][1],
    output: 
        temp('svs/{ref}/{anything}/{anything}.vcf'),
    params:
        ref = lambda wildcards: ref_paths[wildcards.ref],
        anything = lambda wildcards: wildcards.anything,
    shell: """
    svim-asm diploid 'svs/{wildcards.ref}/{wildcards.anything}' {input.hap1} {input.hap2} {params.ref} --min_sv_size 40 --query_names --interspersed_duplications_as_insertions &&
    mv 'svs/{wildcards.ref}/{wildcards.anything}/variants.vcf' 'svs/{wildcards.ref}/{wildcards.anything}/{wildcards.anything}.vcf' &&
    sed -i 's|Sample|{wildcards.anything}|g' 'svs/{wildcards.ref}/{wildcards.anything}/{wildcards.anything}.vcf'
    """

#SVs were called using svim-asm diploid with --query_names --interspersed_duplications_as_insertions --min_sv_size 40. 
#The resulting VCF files were sorted and indexed using BCFtools.

rule bcftools_sort: # getting bzip for truvari, and unzipped for svasm_survivor as none of them handles the other way around
    input:
        'svs/{ref}/{anything}/{anything}.vcf',
    output:
        'svs/{ref}/{anything}/{anything}.sorted.vcf', 
        'svs/{ref}/{anything}/{anything}.sorted.vcf.gz', 
    shell: """
    bcftools sort {input} -Ov -o {output[0]} && 
    bcftools sort {input} -Ob -o {output[1]} && 
    bcftools index -t {output[1]}
    """ 

#rule list_vcfs_panpanhprc2ref:
#    output: 'svs/svasm_survivor/{ref}/panpanhprc_vcf_files_raw_calls.txt'
#    params:
#        all_samples = all_samples  # Use the previously defined
#    run:
#        shell("rm -f {output}")  # Remove existing output file if exists
#        for sample in params.all_samples:
#            shell("ls svs/{wildcards.ref}/{sample}/*.sorted.vcf >> {output}")

#rule list_vcfs_panpan2ref:
#    output: 'svs/svasm_survivor/{ref}/panpan_vcf_files_raw_calls.txt'
#    params:
#      #  all_samples = preferred_samples  # Use the previously defined
#        all_samples = preferred_chimp_base_sample_names  # Use the previously defined
#    run:
#        shell("rm -f {output}")  # Remove existing output file if exists
#        for sample in params.all_samples:
#            shell("ls svs/{wildcards.ref}/{sample}/*.sorted.vcf >> {output}")

rule list_vcfs_pt2ref:
    output: 'svs/{ref}/panpan-pt_vcf_files_raw_calls.txt'
    params:
      #  all_samples = preferred_samples  # Use the previously defined
        all_samples = chimpsamples # Use the previously defined
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.all_samples:
            shell("ls svs/{wildcards.ref}/{sample}/*.sorted.vcf >> {output}")

rule list_vcfs_ppa2ref:
    output: 'svs/{ref}/panpan-ppa_vcf_files_raw_calls.txt'
    params:
      #  all_samples = preferred_samples  # Use the previously defined
        all_samples = bonobosamples
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.all_samples:
            shell("ls svs/{wildcards.ref}/{sample}/*.sorted.vcf >> {output}")


rule list_vcfs_hprc2ref:
    output: 'svs/{ref}/hprc_vcf_files_raw_calls.txt'
    params:
        #all_samples = human_samples  # Use the previously defined 
        all_samples = human_base_sample_names
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.all_samples:
            shell("ls svs/{wildcards.ref}/{sample}/*.sorted.vcf >> {output}")


#rule merge_vcfs_svasm_survivor:
#    input: 'svs/{ref}/{anything}_vcf_files_raw_calls.txt'
#    output: 'svs/svasm_survivor/{ref}/{anything}.merged.vcf' 
#    shell: """
#    ./SURVIVOR merge {input} 1000 1 1 0 0 30 {output}
#    """
#Max distance between breakpoints: 1000
#Minimum number of supporting caller: 1
#Take the type into account (1==yes, else no): 1
#Take the strands of SVs into account (1==yes, else no): 0
#Estimate distance based on the size of SV (1==yes, else no): 0
#Minimum size of SVs to be taken into account: 30

rule merge_vcfs_truvari:
    input: 
        #vcf_files=expand('svs/{ref}/{samples}/{samples}.sorted.vcf.gz', ref=['ht2t'], samples=human_base_sample_names)
        #vcf_files=expand('svs/{ref}/{samples}/{samples}.sorted.vcf.gz', ref=['mPanTro3'], samples=chimpsamples)
        #vcf_files=expand('svs/{ref}/{samples}/{samples}.sorted.vcf.gz', ref=['mPanPan1'], samples=bonobosamples)
        vcf_files=expand('svs/{ref}/{samples}/{samples}.sorted.vcf.gz', ref=['mPonAbe1'], samples=all_samples)
    output: 
        'svs/svasm_truvari/{ref}/{anything}.bcftools.merged.vcf.gz',
        'svs/svasm_truvari/{ref}/{anything}.truvari.merged.vcf.gz', 
        'svs/svasm_truvari/{ref}/{anything}.truvari.collapsed.vcf.gz'
    params:
        ref = lambda wildcards: ref_paths[wildcards.ref]
    shell: """
    bcftools merge -m none {input.vcf_files} | bgzip > {output[0]} &&
    bcftools index -t {output[0]} &&
    truvari collapse -i {output[0]} -o {output[1]} -c {output[2]}  -f {params.ref}
    """

#### run pav
#rule pav:
#    input:
#        config = "pav/config.json",
#        asm_table = "pav/assemblies.tsv",
#    output:
#        # placeholder for run finish right now
#        flag = "pav/run.success"
#    threads: 40
#    params:
#        mount_dir = "/global/scratch/users/joana_rocha",
#        out_dir = "pav",
#        cache_dir = '/global/scratch/users/joana_rocha/software/singularity'
#    shell:
#        """
#        cd {params.out_dir}
#        APPTAINER_CACHEDIR={params.cache_dir}
#
#        singularity run --writable-tmpfs \
#        --bind {params.mount_dir}:/{params.mount_dir} \
#        docker://becklab/pav:latest \
#        --rerun-incomplete \
#        --cores {threads}
#
#        cd -
#        touch {output.flag}
#        """


#rules to merge svim-asm+pbsv+sniffles sv-called multisample vcfs: 
#rule list_vcfs_panpanhprc_svasm_survivor_pbsv_sniffles:
#    output: 'svs/svasm_survivor/panpanhprc_{anything}_survivor_pbsv_sniffles_raw_calls.txt'
#    shell: """
#    ls svs/svasm_survivor/panpanhprc_{wildcards.anything}_merged* >> {output}
#    """

#rule merge_vcfs_panpanhprc_survivor_pbsv_sniffles:
#    input: 'svs/survivor/panpanhprc_{anything}_survivor_pbsv_sniffles_raw_calls.txt'
#    output: 'svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers.vcf'
#    shell: """
#    SURVIVOR merge {input} 1000 1 1 1 0 50 {output}
#    """

#rule make_VennDiagram_input:
#    input: 'svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers.vcf'
#    output: 'svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers_overlapp.txt'
#    shell: """
#    perl perl_script.pl {input} {output}
#    """

#rule venn_diagram:
#    input: 
#        t = "svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers_overlapp.txt"
#    output: 
#        "svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers_overlapp.tiff"
#    params:
#        r_script = "venn_diagram.R"  
#    shell: """
#    Rscript {params.r_script} {input.t} {output}
#    """


### recode vcfs for pca plots
rule recode_genotypes:
    input:
        "svs/survivor/{ref}/{anything}_merged.vcf"
    output:
        "svs/survivor/{ref}/{anything}_merged_modified.vcf"
    shell: "python recode_vcf.py {input} {output}"


rule pcas:
    input:
        'svs/survivor/{ref}/{anything}_merged_modified.vcf'
    output: 
        'plink/{ref}/{anything}_svs.eigenvec',
        'plink/{ref}/{anything}_svs.eigenval'
    params: 'plink/{ref}/{anything}_svs'
    shell: """
    plink --vcf {input} --pca  --double-id --allow-extra-chr --out {params}
    """
