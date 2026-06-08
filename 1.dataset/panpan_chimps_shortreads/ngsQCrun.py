import numpy as np
# to get function np.unique (removes repetitions in a list)
# for example np.unique(['sample1', 'sample2', 'sample2', 'sample3']) -> ['sample1', 'sample2', 'sample3']

 # The following workflow applies to Paired-end Illumina NovaSeq data.
 # Please ensure you have the following working directories for script to work as it is. Otherwise just change it to incorporate your own folder name
 # ### fastq_raw -> folder containing all your fastq.gz folders
 # ### fastq_clean-> directory that will contain your trimmed folders
 # ### bam_fil -> where mapped, filtered and sorted files will be
 # ### bam_rmdup -> duplicate removed bam files
 # ### bam_indel_remap -> GATK output after running RealignTargetCreater and IndelRealigner
 # ### bam_map -> mapdamage outputs

READS = os.listdir('fastq_raw')
SAMPLE_RUNS = np.unique([('_').join(f.split('_')[:4]) for f in READS])
SAMPLES = np.unique([sr.split('_')[0] for sr in SAMPLE_RUNS])
print(SAMPLE_RUNS)
print(SAMPLES)
#ref_path = '/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fna'
#ref_path = '/global/scratch/users/joana_rocha/PANPAN/output/hifiasm-fasta/PR01227/joana_settings_shortcut/no_opts/PR01227.p_ctg.fa'
ref_path = '/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta'

#please make sure you ran the following comands on ref_path:

#bwa index /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta
#samtools faidx /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta
#awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta.fai > /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.bed
#java -jar /global/scratch/users/joana_rocha/software/picard-2.27.2/picard.jar CreateSequenceDictionary R=/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fna O=/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.dict


rule all:
    input:
        expand('fastq_clean/{sample_run}.trim_1P.fastq.gz', sample_run=SAMPLE_RUNS), 
        expand('fastq_clean/{sample_run}.trim_2P.fastq.gz', sample_run=SAMPLE_RUNS), 
        expand('bam_filt/{sample_run}.flagfilt.sorted.bam', sample_run = SAMPLE_RUNS),
        expand('bam_rmdup/{sample_run}.flagfilt.rmdup.sorted.rg.bam', sample_run=SAMPLE_RUNS),
        expand('bam_final_mapped2mPanTro3/{sample}.sorted.bam', sample=SAMPLES)
       
# 1. TRIM READS

#rule cut_adapt:
#    input:
#        'fastq_raw/{sample_run}_1.fq.gz',
#        'fastq_raw/{sample_run}_2.fq.gz'
#    output:
#        temp('fastq_clean/{sample_run}_trim_1.fastq.gz'),
#        temp('fastq_clean/{sample_run}_trim_2.fastq.gz')
#    shell: 'cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -n 2 -e 0.1 -O 1 -m 30 -q 20,20 --max-n 0.5 --nextseq-trim=20 --pair-filter any -o {output[0]} -p {output[1]} {input[0]} {input[1]}'

# use --nextseq-trim only if you have NovaSeq data


rule trim_clean:
    input:
        'fastq_raw/{sample_run}_1.fq.gz',
        'fastq_raw/{sample_run}_2.fq.gz'
    output:
        'fastq_clean/{sample_run}.trim_1P.fastq.gz',
        temp('fastq_clean/{sample_run}.trim_1U.fastq.gz'),
        'fastq_clean/{sample_run}.trim_2P.fastq.gz',
        temp('fastq_clean/{sample_run}.trim_2U.fastq.gz')
    shell:'java -jar /global/scratch/users/joana_rocha/software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -quiet -threads 8 -phred33 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} SLIDINGWINDOW:4:20 MINLEN:30'

# 2. MAP READS WHILE FILTERING
rule map_filt:
    input:
        'fastq_clean/{sample_run}.trim_1P.fastq.gz',
        'fastq_clean/{sample_run}.trim_2P.fastq.gz'
    output: temp('bam_filt/{sample_run}.flagfilt.bam')
    shell: 'bwa mem -t 12 {ref_path} {input[0]} {input[1]} | samtools view -@ 12 -q 15 -bT {ref_path} -F 780 -o {output}'

rule map_filt_sort:
    input: 'bam_filt/{sample_run}.flagfilt.bam'
    output: 'bam_filt/{sample_run}.flagfilt.sorted.bam'
    shell: 'samtools sort -o {output} {input}'

# 3. MARK DUPLICATES
rule picard_rmdup:
    input: 'bam_filt/{sample_run}.flagfilt.sorted.bam'
    output: temp('bam_rmdup/{sample_run}.flagfilt.rmdup.bam')
    shell: """
    mkdir bam_rmdup/tmp_0_{wildcards.sample_run};
    java -Xmx31G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=bam_rmdup/tmp_0_{wildcards.sample_run} -jar /global/scratch/users/joana_rocha/software/picard-2.27.2/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I={input}  O={output} M=bam_rmdup/{wildcards.sample_run}.rmdup.metrics
    """


rule rmdup_sort:
    input: 'bam_rmdup/{sample_run}.flagfilt.rmdup.bam'
    output: temp('bam_rmdup/{sample_run}.flagfilt.rmdup.sorted.bam')
    shell: 'samtools sort -o {output} {input}'


# 4. ADD READGROUPS AND MERGE LIBRARIES FROM DIFFERENT RUNS INTO SAMPLES
rule add_readgroup:
    input: 'bam_rmdup/{sample}_{library}_{flowcell}_{lane}.flagfilt.rmdup.sorted.bam'
    output: 'bam_rmdup/{sample}_{library}_{flowcell}_{lane}.flagfilt.rmdup.sorted.rg.bam'
    shell:
        'samtools addreplacerg -r ID:{wildcards.library}_{wildcards.flowcell}_{wildcards.lane} -r SM:{wildcards.sample} -r PU:{wildcards.lane} -r PL:ILLUMINA -o {output} {input}'


rule merge_to_samples:
    input: expand('bam_rmdup/{sample_run}.flagfilt.rmdup.sorted.rg.bam', sample_run=SAMPLE_RUNS)
    output: expand('bam_merge/{sample}.flagfilt.rmdup.sorted.rg.merged.bam', sample=SAMPLES)
    run:
        for sample in SAMPLES:
            input = 'bam_rmdup/{}_*.rg.bam'.format(sample)
            output = 'bam_merge/{}.flagfilt.rmdup.sorted.rg.merged.bam'.format(sample)
            cmd = 'samtools merge -f {} {}'.format(output, input)
            print(cmd)
            shell(cmd)


# 5. INDEL REAMMAPPING
# java -jar /home/joana/software/ngsQC/picard.jar CreateSequenceDictionary R=/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fastaa O=/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.dict  is necessary for gatk to work

rule indel_remap:
    input: 'bam_merge/{sample}.flagfilt.rmdup.sorted.rg.merged.bam'
    output:
        'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.intervals',
        'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam'
    shell:"""
    samtools index {input} &&
    mkdir bam_indelremap/tmp_0_{wildcards.sample} &&
    /global/scratch/users/joana_rocha/software/jre1.8.0_211/bin/java -Xmx5G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_indelremap/tmp_0_{wildcards.sample} -jar /global/scratch/users/joana_rocha/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {ref_path} -nt 16 -I {input} -o {output[0]} &&
    /global/scratch/users/joana_rocha/software/jre1.8.0_211/bin/java -Xmx4g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_indelremap/tmp_0_{wildcards.sample} -jar /global/scratch/users/joana_rocha/software//GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T IndelRealigner -I {input} -R {ref_path} -targetIntervals {output[0]} -o {output[1]}
    """


# 6. Create final files for snpCleaner
#the input either cames from bam_indelremap or bam_damage, depending on damage patterns
rule bam_final:
    input: 'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam'
    output: 'bam_final_mapped2mPanTro3/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """


rule plotcoverage:
    input: 
        bams = expand('bam_final_mapped2mPanTro3/{sample}.sorted.bam', sample=SAMPLES)
    output: 
        plot = 'coverage_stats/bam_final_mapped2mPanTro3_coverage_plot.pdf'
    params:
        bam_list = lambda wildcards, input: ' '.join(input.bams)
    shell: """
    plotCoverage --bamfiles {params.bam_list} --plotFile {output.plot} -n 1000000 --plotTitle "Chimpanzees (short-reads) mapped to mPanTro3" --outRawCounts coverage.tab --ignoreDuplicates --plotFileFormat pdf --plotHeight 20 --plotWidth 25 --minMappingQuality 15 -p 24 
    """