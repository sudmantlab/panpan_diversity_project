import subprocess
import pandas as pd
import os
import numpy as np

# --- Configuration ---

# File Paths
#VCF_FILE = "lifted_pantros_mapped2mPanTro3_wholegenome.ALLSITES.sorted.filtered.vcf.gz"
#MAIN_POP_FILE = "pop.txt"
VCF_FILE = "../pantros_mapped2mPanTro3_wholegenome.ALLSITES.sorted.filtered_justSHORTEADS.recode.vcf.gz"
MAIN_POP_FILE = "pop_short.txt"
OUTPUT_DIR = "pixy_results"
ANNOTATION_FILE = "mPanTro3_gene.bed.gz" # Annotation file

# Chromosome List
CHROMS_ALIAS_FILE = "chrom_alias_mPanTro3.tsv"
try:
    # The VCF header confirms that the 'Chrom' column contains the correct names.
    chrom_df = pd.read_csv(CHROMS_ALIAS_FILE, sep='\t')
    CHROMS = chrom_df['Chrom'].tolist()
except FileNotFoundError:
    raise FileNotFoundError(f"Error: The chromosome alias file '{CHROMS_ALIAS_FILE}' was not found.")
except KeyError:
    raise KeyError(f"Error: Could not find a column named 'Chrom' in '{CHROMS_ALIAS_FILE}'. Please check the header.")

# --- Target Rule ---
rule all:
    input:
        #f"{OUTPUT_DIR}/annotated_pi_ratios_final.csv"
        f"{OUTPUT_DIR}/pi/within_pops_pi.txt"




# --- Rule 1: Calculate Pi per chromosome for all populations ---
rule pixy_pi_within_per_chrom:
    input:
        vcf=VCF_FILE,
        pops=MAIN_POP_FILE
    output:
        pi_file=f"{OUTPUT_DIR}/pi/within/{{chrom}}_pi.txt"
    log:
        f"{OUTPUT_DIR}/logs/pixy_pi_within.{{chrom}}.log"
    shell:
        """
        echo "Running pixy for within-population Pi on chromosome {wildcards.chrom}..."
        pixy \
            --stats pi \
            --vcf {input.vcf} \
            --populations {input.pops} \
            --chromosomes {wildcards.chrom} \
            --window_size 10000 \
            --output_folder $(dirname {output.pi_file}) \
            --output_prefix {wildcards.chrom} > {log} 2>&1
        """

# --- Rule 2: Create a temporary pooled population file ---
rule create_pooled_popfile:
    input:
        MAIN_POP_FILE
    output:
        temp(f"{OUTPUT_DIR}/pooled_pop.txt")
    shell:
        """
        awk '{{print $1 "\\tpooled"}}' {input} > {output}
        """

# --- Rule 3: Calculate Pi per chromosome for the pooled population ---
rule pixy_pi_pooled_per_chrom:
    input:
        vcf=VCF_FILE,
        pops=f"{OUTPUT_DIR}/pooled_pop.txt"
    output:
        pi_file=f"{OUTPUT_DIR}/pi/pooled/{{chrom}}_pi.txt"
    log:
        f"{OUTPUT_DIR}/logs/pixy_pi_pooled.{{chrom}}.log"
    shell:
        # CORRECTED: The window size must match the within-population rule for merging.
        """
        echo "Running pixy for pooled population Pi on chromosome {wildcards.chrom}..."
        pixy \
            --stats pi \
            --vcf {input.vcf} \
            --populations {input.pops} \
            --chromosomes {wildcards.chrom} \
            --window_size 10000 \
            --output_folder $(dirname {output.pi_file}) \
            --output_prefix {wildcards.chrom} > {log} 2>&1
        """

# --- Rule 4: Aggregate per-chromosome "within population" results ---
rule aggregate_within_pi:
    input:
        expand(f"{OUTPUT_DIR}/pi/within/{{chrom}}_pi.txt", chrom=CHROMS)
    output:
        f"{OUTPUT_DIR}/pi/within_pops_pi.txt"
    shell:
        """
        # Concatenate files, keeping the header from the first file only
        awk 'FNR==1 && NR!=1 {{ next; }} 1' {input} > {output}
        """

# --- Rule 5: Aggregate per-chromosome "pooled population" results ---
rule aggregate_pooled_pi:
    input:
        expand(f"{OUTPUT_DIR}/pi/pooled/{{chrom}}_pi.txt", chrom=CHROMS)
    output:
        f"{OUTPUT_DIR}/pi/pooled_pi.txt"
    shell:
        """
        # Concatenate files, keeping the header from the first file only
        awk 'FNR==1 && NR!=1 {{ next; }} 1' {input} > {output}
        """

# --- Rule 6: Calculate Pi Ratios using a self-contained pandas script ---
rule calculate_pi_ratio:
    input:
        within_pi=f"{OUTPUT_DIR}/pi/within_pops_pi.txt",
        pooled_pi=f"{OUTPUT_DIR}/pi/pooled_pi.txt"
    output:
        summary_csv=f"{OUTPUT_DIR}/pi_ratios_summary.csv"
    run:
        within_df = pd.read_csv(input.within_pi, sep='\\t', engine='python')
        pooled_df = pd.read_csv(input.pooled_pi, sep='\\t', engine='python')
        
        # Proceed only if there's data to merge
        if not within_df.empty and not pooled_df.empty:
            pooled_subset = pooled_df[['chromosome', 'window_pos_1', 'window_pos_2', 'avg_pi']].rename(
                columns={'avg_pi': 'pi_pooled'}
            )
            merged_df = pd.merge(
                within_df,
                pooled_subset,
                on=['chromosome', 'window_pos_1', 'window_pos_2'],
                how='left'
            )
            merged_df['pi_ratio'] = merged_df['avg_pi'] / merged_df['pi_pooled']
            merged_df.to_csv(output.summary_csv, index=False)
        else:
            # If inputs are empty, create an empty file with just the header
            header = "pop,chromosome,window_pos_1,window_pos_2,avg_pi,no_sites,count_diffs,count_comparisons,count_missing,pi_pooled,pi_ratio\n"
            with open(output.summary_csv, 'w') as f:
                f.write(header)


# --- CORRECTED Rule 7: Annotate windows using a robust pandas + bedtools approach ---
rule annotate_pi_with_bedtools_adapted:
    input:
        metrics=f"{OUTPUT_DIR}/pi_ratios_summary.csv",
        genes=ANNOTATION_FILE
    output:
        f"{OUTPUT_DIR}/annotated_pi_ratios_final.csv"
    run:
        temp_files = [
            "temp_metrics.bed",
            "temp_metrics_sorted.bed",
            "temp_genes_sorted.bed",
            "temp_annotated.bed"
        ]

        try:
            # Step 1: Read the summary CSV with pandas, which handles blank fields correctly.
            metrics_df = pd.read_csv(input.metrics)

            # Step 2: Create a new DataFrame in the correct BED format for bedtools.
            # This ensures a consistent number of columns.
            bed_df = pd.DataFrame({
                'chromosome': metrics_df['chromosome'],
                'start': metrics_df['window_pos_1'] - 1, # Convert to 0-based
                'end': metrics_df['window_pos_2'],
                'pop': metrics_df['pop'],
                'avg_pi': metrics_df['avg_pi'],
                'no_sites': metrics_df['no_sites'],
                'count_diffs': metrics_df['count_diffs'],
                'count_comparisons': metrics_df['count_comparisons'],
                'count_missing': metrics_df['count_missing'],
                'pi_pooled': metrics_df['pi_pooled'],
                'pi_ratio': metrics_df['pi_ratio']
            })
            
            # Step 3: Write the clean, tab-delimited BED file.
            bed_df.to_csv("temp_metrics.bed", sep='\t', header=False, index=False, na_rep='.')

            # Step 4: Sort the metrics and gene files
            shell("bedtools sort -i temp_metrics.bed > temp_metrics_sorted.bed")
            shell("zcat {input.genes} | grep -v '^#' | bedtools sort -i - > temp_genes_sorted.bed")
            
            # Step 5: Annotate with gene name (column 4 from the gene file)
            shell(
                "bedtools map "
                "-a temp_metrics_sorted.bed "
                "-b temp_genes_sorted.bed "
                "-c 4 -o collapse > temp_annotated.bed"
            )

            # Step 6: Create the new, reordered header
            new_header = "chromosome,window_pos_1,window_pos_2,pop,avg_pi,no_sites,count_diffs,count_comparisons,count_missing,pi_pooled,pi_ratio,name"
            
            # Step 7: Write the header to the final output file
            with open(output[0], 'w') as f:
                f.write(new_header + '\n')

            # Step 8: Reformat the annotated data into the new column order and append
            shell(r"""
            awk -F'\t' 'BEGIN {{OFS=","}} 
            {{print $1, $2+1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' temp_annotated.bed >> {output}
            """)

        finally:
            # Cleanup
            for f in temp_files:
                if os.path.exists(f):
                    os.remove(f)