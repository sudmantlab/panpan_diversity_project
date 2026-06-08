configfile: "config.yaml"
import pandas as pd
import numpy as np

# Read sample names
SAMPLES = []
with open(config["samples"], "r") as f:
    SAMPLES = [line.strip() for line in f.readlines() if line.strip()]

rule all:
    input:
        "diversity_stats/vcftools_pi/combined_pi_stats.tsv",


# --- VCF Prep ---
rule bgzip_vcf:
    input: config["vcf"]
    output: "diversity_stats/pixy_input/vcf.bgz"
    shell: "mkdir -p diversity_stats/pixy_input && zcat {input} | bgzip > {output}"

rule tabix_index:
    input: "diversity_stats/pixy_input/vcf.bgz"
    output: "diversity_stats/pixy_input/vcf.bgz.tbi"
    shell: "tabix -p vcf {input}"

# --- vcftools ---
rule vcftools_pi_per_sample:
    input: config["vcf"]
    output: "diversity_stats/vcftools_pi/{sample}_pi.sites.pi"
    params: out_prefix = "diversity_stats/vcftools_pi/{sample}_pi"
    log: "logs/vcftools_pi/{sample}.log"
    shell: """
        mkdir -p diversity_stats/vcftools_pi logs/vcftools_pi
        vcftools --gzvcf {input} \
                 --indv {wildcards.sample} \
                 --site-pi \
                 --out {params.out_prefix} \
                 2> {log}
        """

rule combine_vcftools_pi:
    input: expand("diversity_stats/vcftools_pi/{sample}_pi.sites.pi", sample=SAMPLES)
    output: "diversity_stats/vcftools_pi/combined_pi_stats.tsv"
    run:
        stats = []
        for pi_file in input:
            sample = pi_file.split('/')[-1].replace('_pi.sites.pi', '')
            try:
                df = pd.read_csv(pi_file, sep='\t', comment='#')
                if not df.empty:
                    pi_values = df['PI']
                    n_sites = len(pi_values)
                    mean_pi = np.mean(pi_values)
                    std_pi = np.std(pi_values, ddof=1)
                    stderr_pi = std_pi / np.sqrt(n_sites)
                    
                    stats.append({
                        'sample': sample,
                        'mean_pi': mean_pi,
                        'std_pi': std_pi,
                        'stderr_pi': stderr_pi,
                        'n_sites': n_sites
                    })
            except Exception as e:
                print(f"Error processing {sample}: {e}")
        
        pd.DataFrame(stats).to_csv(output[0], sep='\t', index=False, float_format='%.6g')



#rule vcftools_het_all:
#    input: config["vcf"]
#    output: 
#        het_file = "diversity_stats/het_per_indiv/all_samples.het",
#        log_file = "logs/het_all.log"
#    shell: """
#        mkdir -p diversity_stats/het_per_indiv logs
#        vcftools --gzvcf {input} \
#                 --het \
#                 --out diversity_stats/het_per_indiv/all_samples \
#                 2> {output.log_file}
#        """