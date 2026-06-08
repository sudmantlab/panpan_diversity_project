configfile: "config_hsa.yaml"
import pandas as pd
import numpy as np

# Read sample names
SAMPLES = []
with open(config["samples"], "r") as f:
    SAMPLES = [line.strip() for line in f.readlines() if line.strip()]

print(SAMPLES)

rule all:
    input:
        "diversity_stats/vcftools_pi/hsa/combined_pi_stats_nonAFR.tsv",


rule vcftools_pi_per_sample:
    input: config["vcf"]
    output: "diversity_stats/vcftools_pi/hsa/{sample}_pi.sites.pi"
    params: out_prefix = "diversity_stats/vcftools_pi/hsa/{sample}_pi"
    log: "logs/vcftools_pi/hsa/{sample}.log"
    shell: """
        mkdir -p diversity_stats/vcftools_pi/hsa logs/vcftools_pi/hsa
        vcftools --gzvcf {input} \
                 --indv {wildcards.sample} \
                 --site-pi \
                 --out {params.out_prefix} \
                 2> {log}
        """

#rule combine_vcftools_pi_nonAFR:
#    input: expand("diversity_stats/vcftools_pi/hsa/{sample}_pi.sites.pi", sample=SAMPLES)
#    output: "diversity_stats/vcftools_pi/hsa/combined_pi_stats_nonAFR.tsv"
#    run:
#        stats = []
#        for pi_file in input:
#            sample = pi_file.split('/')[-1].replace('_pi.sites.pi', '')
#            try:
#                df = pd.read_csv(pi_file, sep='\t', comment='#')
#                if not df.empty:
#                    pi_values = df['PI']
#                    n_sites = len(pi_values)
#                    mean_pi = np.mean(pi_values)
#                    std_pi = np.std(pi_values, ddof=1)
#                    stderr_pi = std_pi / np.sqrt(n_sites)
#                    
#                    stats.append({
#                        'sample': sample,
#                        'pop': 'nonAFR',  # Added population column
#                        'mean_pi': mean_pi,
#                        'std_pi': std_pi,
#                        'stderr_pi': stderr_pi,
#                        'n_sites': n_sites
#                    })
#            except Exception as e:
#                print(f"Error processing {sample}: {e}")
#        
#        pd.DataFrame(stats).to_csv(output[0], sep='\t', index=False, float_format='%.6g')

rule combine_vcftools_pi_nonAFR:
    input:
        # No input declaration here - we'll handle dynamically in run: section
    output: 
        "diversity_stats/vcftools_pi/hsa/combined_pi_stats_nonAFR.tsv"
    run:
        # Get all possible file paths
        all_files = expand("diversity_stats/vcftools_pi/hsa/{sample}_pi.sites.pi", sample=SAMPLES)
        
        # Filter to only existing files
        existing_files = [f for f in all_files if os.path.exists(f)]
        print(f"Processing {len(existing_files)} existing files out of {len(all_files)} possible")
        
        stats = []
        for pi_file in existing_files:  # Only process existing files
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
                        'pop': 'nonAFR',
                        'mean_pi': mean_pi,
                        'std_pi': std_pi,
                        'stderr_pi': stderr_pi,
                        'n_sites': n_sites
                    })
            except Exception as e:
                print(f"Error processing {sample}: {e}")
        
        pd.DataFrame(stats).to_csv(output[0], sep='\t', index=False, float_format='%.6g')
