import numpy as np
import pandas as pd


CHROMS = [str(c) for c in range(1, 23)] + ["X"]

rule all:
    input:
        #"hsa_hgdp.tgp.wholegenome_AFR.sorted.filtered.vcf.gz",
        #"hsa_hgdp.tgp.wholegenome_AFR.sorted.filtered.vcf.gz.tbi",
        "hsa_hgdp.tgp.wholegenome_nonAFR.sorted.filtered.vcf.gz",
        "hsa_hgdp.tgp.wholegenome_nonAFR.sorted.filtered.vcf.gz.tbi"

rule filter_per_chrom:
    input:
        bcf = "hgdp.tgp.gwaspy_phased_haplotypes/hgdp.tgp.gwaspy.merged.chr{chr}.merged.bcf",
        csi = "hgdp.tgp.gwaspy_phased_haplotypes/hgdp.tgp.gwaspy.merged.chr{chr}.merged.bcf.csi",
        #samples = "AFR.txt"
        samples = "nonAFR.txt"
    output:
        vcf = temp("filtered_per_chr/chr{chr}.filtered.vcf.gz"),
        tbi = temp("filtered_per_chr/chr{chr}.filtered.vcf.gz.tbi")
    params:
        tmpdir = "/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/bcftools_tmp/"
    threads: 4
    shell:
        """
        mkdir -p {params.tmpdir}
        # Filter and create VCF.gz
        bcftools view -S {input.samples} \
            --exclude-types indels \
            -m 1 -M 2 \
            {input.bcf} \
            -Oz -o {output.vcf} \
            --threads {threads}
        
        # Index the filtered VCF
        bcftools index -t {output.vcf} \
            --threads {threads}
        """

rule concat_all:
    input:
        vcfs = expand("filtered_per_chr/chr{chr}.filtered.vcf.gz", chr=CHROMS),
        indexes = expand("filtered_per_chr/chr{chr}.filtered.vcf.gz.tbi", chr=CHROMS)
    output:
        temp("hsa_hgdp.tgp.wholegenome.unsorted.vcf.gz")
    shell:
        """
        bcftools concat {input.vcfs} \
            -a \
            -Oz -o {output} \
            --threads 8
        """

rule sort_and_index:
    input:
        "hsa_hgdp.tgp.wholegenome.unsorted.vcf.gz"
    output:
        #"hsa_hgdp.tgp.wholegenome_AFR.sorted.filtered.vcf.gz",
        #"hsa_hgdp.tgp.wholegenome_AFR.sorted.filtered.vcf.gz.tbi"
        "hsa_hgdp.tgp.wholegenome_nonAFR.sorted.filtered.vcf.gz",
        "hsa_hgdp.tgp.wholegenome_nonAFR.sorted.filtered.vcf.gz.tbi"
    shell:
        """
        mkdir -p {params.tmpdir}
        bcftools sort {input} \
            --output-type z \
            -o {output[0]} \
        bcftools index {output[0]} \
        """