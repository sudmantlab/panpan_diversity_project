import numpy as np
import pandas as pd
import os
import itertools
from snakemake.io import glob_wildcards


ref_paths = {
    'mPonAbe1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPonAbe1/mPonAbe1.pri.cur.20231205.fasta",
    'mPanPan1' : "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"

}

region_table = pd.read_csv('mPanPan1.txt', sep='\t', names= ["chrom"])
chromosomes=region_table['chrom'].unique()
print(region_table)
print(chromosomes)
references = ['mPanPan1']

### chosen vcf for QC
VCF_INPUT = 'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz'
OUT_PREFIX = 'vcfs_from_reads/summary_stats/panpaniscus_mapped2mPanPan1_wholegenome.BIALLELIC_SNPS.sorted.filtered'

rule all:
    input:
        # VCF generation and processing
        #expand("vcfs_from_reads/per_chrom/panpaniscus_mapped2mPanPan1_{chromosome}.vcf.gz", chromosome=chromosomes),   
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.sorted.vcf.gz',
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.ALLSITES.sorted.filtered.vcf.gz',
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz',
        #'vcfs_from_reads/coverage_stats/coverage_bams.txt',
        #'vcfs_from_reads/coverage_stats/coverage_histograms.png',
        # Detailed VCF stats (using bcftools query)
        f"{OUT_PREFIX}.depth.gz",
        f"{OUT_PREFIX}.quality.gz",
        f"{OUT_PREFIX}.missing.gz",
        f"{OUT_PREFIX}.sample_missing.gz",
        f"{OUT_PREFIX}.sample_het.gz",
        #OUT_PREFIX + "/gq_stats.txt",
        #OUT_PREFIX + "/gq_per_sample.txt",
        #OUT_PREFIX + "/gq_per_site.txt",
        # from bcftools stats
        f"{OUT_PREFIX}.bcftools.stats",
        #f"{OUT_PREFIX}/plots/bcftools/quality.png",
        #f"{OUT_PREFIX}/plots/bcftools/depth.png",
        #f"{OUT_PREFIX}/plots/bcftools/indel_distribution.png",
        #f"{OUT_PREFIX}/plots/bcftools/substitution_types.png",
        #f"{OUT_PREFIX}/plots/bcftools/allele_frequency.png",
        #'plink/panpaniscus_nohybrids_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenvec',
        #'plink/panpaniscus_nohybrids_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenval',
        'plink/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenvec',
        'plink/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenval'

rule calculate_coverage_per_bams:
    input:
        bam_file_paths = "panpaniscus_long_short_mapped2mPanPan1_bam_file_paths.txt",
        vcf_file = VCF_INPUT
    output:
        temp("coverage_stats/{sample}.coverage.txt")
    params:
        outdir = "coverage_stats"
    shell:
        """
        set -e  # Exit on error
        set -x  # Print commands for debugging
        
        # Create output directory
        mkdir -p {params.outdir}
        
        # Get the BAM file path for this sample using a more flexible pattern
        bam_file=$(grep -E "/{wildcards.sample}[._]|/{wildcards.sample}\.sorted\.bam" {input.bam_file_paths} | head -n1 || true)
        
        echo "Looking for sample: {wildcards.sample}"
        echo "Found BAM file: $bam_file"
        
        if [ -z "$bam_file" ]; then
            echo "Error: No BAM file found for sample {wildcards.sample}"
            echo "Contents of BAM file paths:"
            cat {input.bam_file_paths}
            exit 1
        fi
        
        # Check if BAM file exists
        if [ ! -f "$bam_file" ]; then
            echo "Error: BAM file $bam_file does not exist"
            exit 1
        fi
        
        # Create genome file from BAM header
        echo "Creating genome file from BAM header..."
        samtools view -H "$bam_file" | \
            grep -P '^@SQ' | \
            cut -f 2,3 | \
            sed 's/SN://;s/LN://' | \
            awk '{{print $1"\t1\t"$2}}' > {params.outdir}/genome.{wildcards.sample}.tmp
        
        # Check if genome file was created and has content
        if [ ! -s {params.outdir}/genome.{wildcards.sample}.tmp ]; then
            echo "Error: genome file is empty or was not created"
            exit 1
        fi
        
        echo "Calculating coverage..."
        # Calculate coverage
        bedtools genomecov -ibam "$bam_file" -g {params.outdir}/genome.{wildcards.sample}.tmp | \
            awk -v sample="{wildcards.sample}" '{{print $0"\t"sample}}' > {output}
        
        # Check if output was created
        if [ ! -s {output} ]; then
            echo "Error: Output file is empty"
            exit 1
        fi
        
        # Clean up
        rm -f {params.outdir}/genome.{wildcards.sample}.tmp
        """

rule concatenate_coverage:
    input:
        coverage_files = lambda wildcards: expand(
            "coverage_stats/{sample}.coverage.txt",
            sample=get_sample_names(VCF_INPUT)
        )
    output:
        "vcfs_from_reads/coverage_stats/coverage_bams.txt"
    shell:
        """
        cat {input.coverage_files} > {output}
        """

rule plot_coverage_histograms:
    input:
        "vcfs_from_reads/coverage_stats/coverage_bams.txt"
    output:
        "vcfs_from_reads/coverage_stats/coverage_histograms.png"
    script:
        "coverage_stats/plot_coverage_histograms.R"



#### all sites per chromosome
rule run_bcftools_bonobos:
    input: 
        bam_files= "panpaniscus_long_short_mapped2mPanPan1_bam_file_paths.txt",
    output: 
        "vcfs_from_reads/per_chrom/panpaniscus_mapped2mPanPan1_{chromosome}.vcf.gz"
    shell: """
    bcftools mpileup --threads 20 -q 30 -Q 20 -a AD,DP,SP -Ou -f /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta -b {input.bam_files} -r {wildcards.chromosome}  | bcftools call -f GQ,GP -mO z -o {output}
    """
#-a - Annotate the vcf - here we add allelic depth (AD), genotype depth (DP) and strand bias (SP).
#-f - format fields for the vcf - here they are genotype quality (GQ) and genotype probability (GP).
#-m- use bcftools multiallelic caller
# -v for variant sites only has not been use to keep all sites
# -q 30 -Q 20  are base and mapping quality filters 

rule index_vcfs:
    input:
        "vcfs_from_reads/per_chrom/panpaniscus_mapped2mPanPan1_{chromosome}.vcf.gz",
    output:
        "vcfs_from_reads/per_chrom/panpaniscus_mapped2mPanPan1_{chromosome}.vcf.gz.tbi"
    shell:
        """
        bcftools index -t {input}
        """
#### all sites whole genome
rule concat_bcfs:
    input:
        vcf_files=expand("vcfs_from_reads/per_chrom/panpaniscus_mapped2mPanPan1_{chromosome}.vcf.gz", chromosome=chromosomes),
        index_files=expand("vcfs_from_reads/per_chrom/panpaniscus_mapped2mPanPan1_{chromosome}.vcf.gz.tbi", chromosome=chromosomes)
    output:
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.vcf.gz',
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.sorted.vcf.gz'
    params:
        tmpdir = "/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/bcftools_tmp/"
    shell: """
    mkdir -p {params.tmpdir} &&
    bcftools concat {input.vcf_files} -a -Oz --threads 20 -o {output[0]} &&
    bcftools sort {output[0]}  --output-type z -o {output[1]} -T {params.tmpdir} &&
    bcftools index -t {output[1]} --threads 24 
    """

rule filter_vcf:
    input:
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.sorted.vcf.gz'
    output:
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.ALLSITES.sorted.filtered.vcf.gz'
    params:
        tmpdir = "/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/bcftools_tmp/"
    shell: """
    mkdir -p {params.tmpdir} &&
    vcftools --gzvcf {input} \
        --remove-indels \
        --max-alleles 2 \
        --max-missing 0.9 \
        --min-meanDP 10 \
        --max-meanDP 55 \
        --minDP 5 \
        --maxDP 200 \
        --recode \
        --recode-INFO-all \
        --stdout | gzip -c > {output}
    """

#depth params rationale set as: 
#min-meanDP 200/18 and max-meanDP 1000/18
#(18 is total long+short-read chimps mapped to mPanPan1)


rule filter_vcf_biallelicsnps:
    input:
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.sorted.vcf.gz'
    output:
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz'
    params:
        tmpdir = "/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/bcftools_tmp/"
    shell: """
    mkdir -p {params.tmpdir} &&
    vcftools --gzvcf {input} \
        --remove-indels \
        --min-alleles 2 \
        --max-alleles 2 \
        --max-missing 0.9 \
        --min-meanDP 10 \
        --max-meanDP 55 \
        --minDP 5 \
        --maxDP 200 \
        --minQ 30 \
        --minGQ 20 \
        --recode \
        --recode-INFO-all \
        --stdout | gzip -c > {output} 
    """

# -minGQ 20 
#--remove-indels - remove all indels (SNPs only)
#--min-alleles 2 --max-alleles 2  (Biallelic SNPS only)
#--max-missing - set minimum non-missing data. A little counterintuitive - 0 is totally missing, 1 is none missing. Here 0.9 means we will tolerate 10% missing data.
#--minQ - this is just the minimum quality score required for a site to pass our filtering threshold. Here we set it to 30.
#--min-meanDP - the minimum mean depth for a site. 
#--max-meanDP - the maximum mean depth for a site. 
#--minDP - the minimum depth allowed for a genotype - any individual failing this threshold is marked as having a missing genotype.
#--maxDP - the maximum depth allowed for a genotype - any individual failing this threshold is marked as having a missing genotype.
#--recode - recode the output - necessary to output a vcf
#--recode-INFO-all (not to loose important information/tags)


def get_sample_names(vcf_file):
    import subprocess
    try:
        result = subprocess.run(
            f"bcftools query -l {vcf_file}",
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        samples = result.stdout.strip().split('\n')
        print(f"Found samples in VCF: {samples}")  # Debug output
        return samples
    except subprocess.CalledProcessError as e:
        print(f"Error getting sample names: {e}")
        return []

rule calculate_per_sample_stats:
    input:
        vcf=VCF_INPUT
    output:
        missing=f"{OUT_PREFIX}.sample_missing.gz",
        het=f"{OUT_PREFIX}.sample_het.gz"
    params:
        out_prefix=lambda wildcards, output: output.missing.replace('.sample_missing.gz', '')
    shell:
        """
        # Calculate missingness per sample
        vcftools --gzvcf {input.vcf} \
            --missing-indv \
            --out {params.out_prefix} && \
        mv {params.out_prefix}.imiss {output.missing}

        # Calculate heterozygosity per sample
        vcftools --gzvcf {input.vcf} \
            --het \
            --out {params.out_prefix} && \
        mv {params.out_prefix}.het {output.het}
        """


rule calculate_stats_per_site:
    input:
        vcf=VCF_INPUT
    output:
        depth=f"{OUT_PREFIX}.depth.gz",
        quality=f"{OUT_PREFIX}.quality.gz",
        missing=f"{OUT_PREFIX}.missing.gz"
    shell:
        """
        # Calculate per-site depth - just INFO/DP
        bcftools query -f '%CHROM\\t%POS\\t%INFO/DP\\n' {input.vcf} | \
        awk '{{
            if($3 != "." && $3 != "") print $0;
            else print $1"\\t"$2"\\t0"
        }}' | gzip > {output.depth}
        
        # Extract site qualities
        bcftools query -f '%CHROM\\t%POS\\t%QUAL\\n' {input.vcf} | gzip > {output.quality}
        
        # Calculate missing data
        bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' {input.vcf} | \
        awk '{{ 
            missing=0; 
            total=NF-2; 
            for(i=3;i<=NF;i++) if($i=="./.") missing++; 
            print $1"\\t"$2"\\t"missing/total 
        }}' | gzip > {output.missing}
        """


rule extract_gq_stats:
    input:
        vcf=VCF_INPUT
    output:
        gq_raw=OUT_PREFIX + "/gq_stats.txt",
        gq_per_sample=OUT_PREFIX + "/gq_per_sample.txt",
        gq_per_site=OUT_PREFIX + "/gq_per_site.txt"
    shell:
        """
        bcftools query -H -f '%CHROM\\t%POS[\\t%GQ]\\n' {input.vcf} > {output.gq_raw}
        Rscript -e '
        data <- read.table("{output.gq_raw}", header=TRUE, check.names=FALSE)
        # Get sample names (excluding CHROM and POS columns)
        samples <- colnames(data)[-(1:2)]
        # Calculate per-sample average GQ
        per_sample_gq <- colMeans(data[,-(1:2)], na.rm=TRUE)
        write.table(data.frame(
            Sample=samples, 
            Mean_GQ=per_sample_gq
        ), "{output.gq_per_sample}", 
        row.names=FALSE, quote=FALSE, sep="\t")
        # Calculate per-site average GQ
        per_site_gq <- data.frame(
            CHROM=data$CHROM,
            POS=data$POS,
            Mean_GQ=rowMeans(data[,-(1:2)], na.rm=TRUE)
        )
        write.table(per_site_gq, "{output.gq_per_site}", 
        row.names=FALSE, quote=FALSE, sep="\t")
        '
        """


rule bcftools_stats:
    input:
        vcf=VCF_INPUT,
        ref=ref_paths['mPanPan1']
    output:
        stats=OUT_PREFIX + ".bcftools.stats"
    threads: 20
    log:
        "logs/bcftools_stats.log"
    shell:
        "bcftools stats --fasta-ref {input.ref} --threads {threads} --depth 0,10000,1 -s - {input.vcf} > {output.stats} 2> {log}"

rule plot_bcf_stats:
    input:
        stats=OUT_PREFIX + ".bcftools.stats"
    output:
        qual_plot=OUT_PREFIX + "/plots/bcftools/quality.png",
        depth_plot=OUT_PREFIX + "/plots/bcftools/depth.png",
        indel_plot=OUT_PREFIX + "/plots/bcftools/indel_distribution.png",
        subst_plot=OUT_PREFIX + "/plots/bcftools/substitution_types.png",
        af_plot=OUT_PREFIX + "/plots/bcftools/allele_frequency.png",
        sample_depth_plot=OUT_PREFIX + "/plots/bcftools/per_sample_depth.png",
        sample_missing_plot=OUT_PREFIX + "/plots/bcftools/per_sample_missing.png",
        sample_het_plot=OUT_PREFIX + "/plots/bcftools/per_sample_heterozygosity.png"
    log:
        "logs/plot_bcf_stats.log"
    shell:
        """
        mkdir -p $(dirname {output.qual_plot})
        Rscript plot_bcf_stats.R {input.stats} {output.qual_plot} {output.depth_plot} \
            {output.indel_plot} {output.subst_plot} {output.af_plot} \
            {output.sample_depth_plot} {output.sample_missing_plot} \
            {output.sample_het_plot} 2> {log}
        """

#--maf - set minor allele frequency - 0.05 

rule pcas_biallelic_chimps:
    input:
        'vcfs_from_reads/panpaniscus_mapped2mPanPan1_wholegenome.BIALLELIC_SNPS.sorted.filtered.vcf.gz'
    output: 
        'plink/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenvec',
        'plink/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenval'
    params: 'plink/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS'
    shell: """
    plink --vcf {input} --pca --maf 0.05 --double-id --allow-extra-chr --make-bed --out {params}
    """


