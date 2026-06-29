rule sniffles_single_sample:
    input:
        bam = get_bam_from_sample,
        ref_path = BASEDIR + '/reference/{ref}.fasta',
        trf_path = BASEDIR + '/reference/{ref}.trf.bed', 
    output:
        vcf = BASEDIR + '/sniffles/{ref}/{sample}.vcf.gz',
        snf = BASEDIR + '/sniffles/{ref}/{sample}.snf',
        tbi = BASEDIR + '/sniffles/{ref}/{sample}.vcf.gz.tbi',
        done = touch(BASEDIR + '/sniffles/{ref}/{sample}.done'),
    conda:
        BASEDIR + '/workflow/envs/sniffles.yaml'
    threads:
        28
    log: BASEDIR + '/sniffles/{ref}/{sample}.log'
    params:
        outdir = BASEDIR + '/sniffles/{ref}',
    shell:
        '''
        mkdir -p {params.outdir}
        sniffles --input {input.bam} --vcf {output.vcf} --snf {output.snf} --ref {input.ref_path} --tandem-repeats {input.trf_path} --threads {threads}  &> {log}
        '''

rule sniffles_multi_sample:
    input:
        snf = get_snfs_from_dataset,
        ref_path = BASEDIR + '/reference/{ref}.fasta',
        trf_path = BASEDIR + '/reference/{ref}.trf.bed', 
    output:
        vcf = BASEDIR + '/sniffles/{ref}/{dataset}.vcf',
        done = touch(BASEDIR + '/sniffles/{ref}/{dataset}.done'),
    conda:
        BASEDIR + '/workflow/envs/sniffles.yaml'
    threads:
        28
    log: BASEDIR + '/sniffles/{ref}/{dataset}.log'
    shell:
        'sniffles --input {input.snf} --vcf {output.vcf} --ref {input.ref_path} --tandem-repeats {input.trf_path} --threads {threads} &> {log}'

rule filter_sniffles:
    input:
        vcf = BASEDIR + '/sniffles/{ref}/{dataset}.vcf',
    output: 
        indel = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.indel.vcf.gz',
        dupinv = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.dupinv.vcf.gz',
        done = touch(BASEDIR + '/sniffles/{ref}/{dataset}.filtered.done'),
    threads: 16
    params:
        indir = BASEDIR + '/sniffles/{ref}',
        rscript = BASEDIR + '/workflow/scripts/filter_sniffles.R'
    log: BASEDIR + '/sniffles/{ref}/{dataset}.filtered.log'
    shell:
        '''
        module load r
        Rscript {params.rscript} {params.indir} {wildcards.dataset} {threads} &> {log}
        '''
        
rule sort_filtered_sniffles:
    input:
        indel = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.indel.vcf.gz',
        dupinv = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.dupinv.vcf.gz',
    output:
        indel = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.indel.sorted.vcf.gz',
        dupinv = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.dupinv.sorted.vcf.gz',
        bcftools_concat = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.sorted.vcf.gz',
        done = touch(BASEDIR + '/sniffles/{ref}/{dataset}.filtered.sorted.done'),
    conda: 'bcftools'
    threads: 28
    shell:
        '''
        gunzip -f -c {input.indel} | bcftools sort | bgzip > {output.indel}
        bcftools index {output.indel}
        gunzip -f -c {input.dupinv} | bcftools sort | bgzip > {output.dupinv}
        bcftools index {output.dupinv}
        bcftools concat -a {output.indel} {output.dupinv} | bgzip > {output.bcftools_concat}
        bcftools index {output.bcftools_concat}
        '''
