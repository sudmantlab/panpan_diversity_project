rule truvari_merge_sample:
    input:
        vcf = get_vcf_from_dataset,
        ref_path = BASEDIR + '/reference/{ref}.fasta',
    output:
        bcftools_merged = BASEDIR + '/svim-asm/{ref}/{dataset}.bcftools.merged.vcf.gz',
        truvari_merged = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.vcf.gz',
        truvari_collapsed = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.collapsed.vcf.gz',
        done = touch(BASEDIR + '/svim-asm/{ref}/{dataset}.merge.done'),
    threads:
        6
    conda:
        "truvari-4.2.0"
    log: BASEDIR + '/svim-asm/{ref}/{dataset}.merge.log'
    shell:
        '''
        bcftools merge -m none {input.vcf} | bgzip > {output.bcftools_merged} 2>> {log}
        bcftools index -t {output.bcftools_merged} &>> {log}
        truvari collapse -i {output.bcftools_merged} -o {output.truvari_merged} -c {output.truvari_collapsed} -f {input.ref_path} --chain -r 1000 --sizemin 50 -p 0.9 -P 0.9 &>> {log}
        '''

rule missing_to_ref:
    input:
        truvari_merged = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.vcf.gz',
    output:
        missing2ref = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.vcf.gz',
        done = touch(BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.done'),
    threads:
        6
    conda:
        "bcftools"
    log: BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.log'
    shell:
        '''
        bcftools +missing2ref  {input.truvari_merged} | bgzip >  {output.missing2ref} 2>> {log}
        '''

rule plink_pca:
    input:
        missing2ref = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.vcf.gz',
    output:
        eigenval = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.eigenval',
        eigenvec = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.eigenvec',
        done = touch(BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.pca.done'),
    params:
        outbase = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref',
    threads:
        6
    conda:
        "plink"
    log: BASEDIR + '/svim-asm/{ref}/{dataset}.missing2ref.pca.log'
    shell:
        '''
        plink --vcf {input.missing2ref} --double-id --allow-extra-chr --pca --out {params.outbase} &> {log}
        '''
        
rule filter_svimasm:
    input:
        vcf = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.vcf.gz',
    output: 
        indel = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.indel.vcf.gz',
        dupinv = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.dupinv.vcf.gz',
        done = touch(BASEDIR + '/svim-asm/{ref}/{dataset}.filtered.done'),
    threads: 16
    params:
        indir = BASEDIR + '/svim-asm/{ref}',
        rscript = BASEDIR + '/workflow/scripts/filter_svimasm.R'
    log: BASEDIR + '/svim-asm/{ref}/{dataset}.filtered.log'
    shell:
        '''
        module load r
        Rscript {params.rscript} {params.indir} {wildcards.dataset} {threads} &> {log}
        '''

rule sort_filtered_svimasm:
    input:
        indel = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.indel.vcf.gz',
        dupinv = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.dupinv.vcf.gz',
    output:
        indel = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.indel.sorted.vcf.gz',
        dupinv = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.dupinv.sorted.vcf.gz',
        bcftools_concat = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.sorted.vcf.gz',
        done = touch(BASEDIR + '/svim-asm/{ref}/{dataset}.filtered.sorted.done'),
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

rule truvari_merge_caller:
    input:
        sniffles_indel = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.indel.sorted.vcf.gz',
        sniffles_dupinv = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.dupinv.sorted.vcf.gz',
        svimasm_indel = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.indel.sorted.vcf.gz',
        svimasm_dupinv = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.dupinv.sorted.vcf.gz',
        ref_path = BASEDIR + '/reference/{ref}.fasta',
    output:
        bcftools_merged_indel = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.bcftools.merged.indel.vcf.gz',
        truvari_merged_indel = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.merged.indel.vcf.gz',
        truvari_collapsed_indel = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.collapsed.indel.vcf.gz',
        truvari_merged_indel_sorted = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.merged.indel.sorted.vcf.gz',
        bcftools_merged_dupinv = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.bcftools.merged.dupinv.vcf.gz',
        truvari_merged_dupinv = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.merged.dupinv.vcf.gz',
        truvari_collapsed_dupinv = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.collapsed.dupinv.vcf.gz',
        truvari_merged_dupinv_sorted = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.merged.dupinv.sorted.vcf.gz',
        bcftools_concat = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.vcf.gz',
        done = touch(BASEDIR + '/truvari/{ref}/{dataset}.merge.done'),
    threads:
        6
    conda:
        "truvari-4.2.0"
    log: BASEDIR + '/truvari/{ref}/{dataset}.merge.log'
    shell:
        '''
        ##indel
        bcftools merge -m none --force-samples {input.sniffles_indel} {input.svimasm_indel} | bgzip > {output.bcftools_merged_indel} 2> {log}
        bcftools index -t {output.bcftools_merged_indel} 2>> {log}
        truvari collapse -i {output.bcftools_merged_indel} -o {output.truvari_merged_indel} -c {output.truvari_collapsed_indel} -f {input.ref_path} --chain -r 1000 --sizemin 50 -p 0.9 -P 0.9 &>> {log}
        bcftools sort {output.truvari_merged_indel} | bgzip > {output.truvari_merged_indel_sorted}  2>> {log}
        bcftools index {output.truvari_merged_indel_sorted}  &>> {log}
        ## dupinv
        bcftools merge -m none --force-samples {input.sniffles_dupinv} {input.svimasm_dupinv} | bgzip > {output.bcftools_merged_dupinv} 2>> {log}
        bcftools index -t {output.bcftools_merged_dupinv} 2> {log}
        truvari collapse -i {output.bcftools_merged_dupinv} -o {output.truvari_merged_dupinv} -c {output.truvari_collapsed_dupinv} -f {input.ref_path} --chain -r 1000 --sizemin 50 -p 0 -P 0.9 &>> {log}
        bcftools sort {output.truvari_merged_dupinv} | bgzip > {output.truvari_merged_dupinv_sorted} 2>> {log}
        bcftools index {output.truvari_merged_dupinv_sorted} &>> {log}
        ## concat
        bcftools concat -a -o {output.bcftools_concat} {output.truvari_merged_indel_sorted} {output.truvari_merged_dupinv_sorted} &>> {log}
        '''

rule generate_indel_sequences:
    input:
        vcf = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.vcf.gz',
        tandem_repeat = BASEDIR + '/reference/{ref}.trf.bed',
    output: 
        fasta = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.indel.fasta',
        vcf_summary_wide = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.summary.tsv',
        sv_id_table = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.id.tsv',
        done = touch(BASEDIR + '/truvari/{ref}/{dataset}.indel_sequences.done'),
    threads: 28
    params:
        indir = BASEDIR + '/truvari/{ref}',
        rscript = BASEDIR + '/workflow/scripts/generate_indel_sequences.R'
    log: BASEDIR + '/truvari/{ref}/{dataset}.indel_sequences.log'
    shell:
        '''
        module load r
        Rscript {params.rscript} {params.indir} {wildcards.dataset} {input.tandem_repeat} {threads} &> {log}
        '''

rule run_repeatmasker:
    input:
        fasta = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.indel.fasta',
    output: 
        gff = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.indel.fasta.out.gff',
        out = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.indel.fasta.out',
        done = touch(BASEDIR + '/truvari/{ref}/{dataset}.repeatmasker.done'),
    threads: 40
    log: BASEDIR + '/truvari/{ref}/{dataset}.repeatmasker.log'
    params:
        indir = BASEDIR + '/truvari/{ref}',
        repeatmasker_database = get_repeatmasker_database_from_dataset,
    conda:
        "repeatmasker"
    shell:
        '''
        cd {params.indir}
        RepeatMasker -pa {threads} -engine rmblast -nocut -gff -species {params.repeatmasker_database} {input.fasta} &> {log}
        '''

rule analyze_repeatmasker:
    input:
        out = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.indel.fasta.out',
        gff = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.indel.fasta.out.gff',
        vcf_summary_wide = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.summary.tsv',
        sv_id_table = BASEDIR + '/truvari/{ref}/{dataset}.sniffles_svim-asm.truvari.sorted.merged.id.tsv',
    output: 
        svimasm_consensus_id = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.consensus.id',
        sniffles_consensus_id = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.consensus.id',
        svimasm_only_id = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.only.id',
        sniffles_only_id = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.only.id',
        done = touch(BASEDIR + '/truvari/{ref}/{dataset}.analyze_repeatmasker.done'),
    threads: 6
    params:
        indir = BASEDIR + '/truvari/{ref}',
        rscript = BASEDIR + '/workflow/scripts/analyze_repeatmasker.R'
    log: BASEDIR + '/truvari/{ref}/{dataset}.analyze_repeatmasker.log'
    shell:
        '''
        module load r
        Rscript {params.rscript} {params.indir} {wildcards.dataset} &> {log}
        '''

rule generate_consensus_callset:
    input:
        svimasm_consensus_id = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.consensus.id',
        sniffles_consensus_id = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.consensus.id',
        svimasm_only_id = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.only.id',
        sniffles_only_id = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.only.id',
        svimasm_vcf = BASEDIR + '/svim-asm/{ref}/{dataset}.truvari.merged.missing2ref.filtered.sorted.vcf.gz',
        sniffles_vcf = BASEDIR + '/sniffles/{ref}/{dataset}.filtered.sorted.vcf.gz',
    output:
        svimasm_consensus_vcf = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.consensus.vcf.gz',
        sniffles_consensus_vcf = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.consensus.vcf.gz',
        svimasm_only_vcf = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.only.vcf.gz',
        sniffles_only_vcf = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.only.vcf.gz',
        concat_vcf = BASEDIR + '/truvari/{ref}/{dataset}.concat.vcf.gz',
        done = touch(BASEDIR + '/truvari/{ref}/{dataset}.consensus_callset.done'),
    conda: 'bcftools'
    threads: 28
    shell:
        '''
        bcftools view -i'ID=@{input.svimasm_consensus_id}' {input.svimasm_vcf} | awk '{{gsub("-svimasm", ""); print}}' | bgzip > {output.svimasm_consensus_vcf}
        bcftools view -i'ID=@{input.sniffles_consensus_id}' {input.sniffles_vcf} | awk '{{gsub("-sniffles", ""); print}}' | bgzip > {output.sniffles_consensus_vcf}
        bcftools view -i'ID=@{input.svimasm_only_id}' {input.svimasm_vcf} | awk '{{gsub("-svimasm", ""); print}}' | bgzip > {output.svimasm_only_vcf}
        bcftools view -i'ID=@{input.sniffles_only_id}' {input.sniffles_vcf} | awk '{{gsub("-sniffles", ""); print}}' | bgzip > {output.sniffles_only_vcf}
        bcftools index {output.svimasm_consensus_vcf}
        bcftools index {output.sniffles_consensus_vcf}
        bcftools index {output.svimasm_only_vcf}
        bcftools index {output.sniffles_only_vcf}
        bcftools concat -a {output.svimasm_consensus_vcf} {output.svimasm_only_vcf} {output.sniffles_only_vcf} | bgzip > {output.concat_vcf}
        bcftools index {output.concat_vcf}
        '''

rule pca_with_consensus_callset:
    input:
        svimasm_consensus_vcf = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.consensus.vcf.gz',
        sniffles_consensus_vcf = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.consensus.vcf.gz',
        svimasm_only_vcf = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.only.vcf.gz',
        sniffles_only_vcf = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.only.vcf.gz',
        concat_vcf = BASEDIR + '/truvari/{ref}/{dataset}.concat.vcf.gz',
    output:
        done = touch(BASEDIR + '/truvari/{ref}/{dataset}.pca_with_consensus_callset.done'),
    params:
        svimasm_consensus = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.consensus',
        sniffles_consensus = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.consensus',
        svimasm_only = BASEDIR + '/truvari/{ref}/{dataset}.svimasm.only',
        sniffles_only = BASEDIR + '/truvari/{ref}/{dataset}.sniffles.only',
        concat = BASEDIR + '/truvari/{ref}/{dataset}.concat',
    conda: 'plink'
    threads: 6
    shell:
        '''
        plink --vcf {input.svimasm_consensus_vcf} --double-id --allow-extra-chr --pca --out {params.svimasm_consensus}
        plink --vcf {input.sniffles_consensus_vcf} --double-id --allow-extra-chr --pca --out {params.sniffles_consensus}
        plink --vcf {input.svimasm_only_vcf} --double-id --allow-extra-chr --pca --out {params.svimasm_only}
        plink --vcf {input.sniffles_only_vcf} --double-id --allow-extra-chr --pca --out {params.sniffles_only}
        plink --vcf {input.concat_vcf} --double-id --allow-extra-chr --pca --out {params.concat}
        '''
