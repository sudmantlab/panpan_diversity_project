from snakemake.io import glob_wildcards
import os
import pandas as pd
import pyranges as pr


NE = "2e4"
M = "1.25e-8"
CHROM_SIZES = "chm13v2.0.sizes"
POPULATIONS = "populations.csv"
SAMPLES = "samples.txt"
GENE_ANNOTATION = "ht2t_CHM13_renamed_chr_gene.bed.gz"
GENOME_FEATURES = "./ht2t_genomefeatures.bed"
NUM_ARGS = 99  # For 100 samples (0-99)


vcf_components = glob_wildcards("vcf_chunks/{chr}_{start}_{end}.vcf")
existing_chunks = list(zip(vcf_components.chr, vcf_components.start, vcf_components.end))
chrs = sorted(list(set(vcf_components.chr))) 

chrom_sizes = {}
with open(CHROM_SIZES) as f:
    for line in f:
        parts = line.strip().split()
        chrom_sizes[parts[0]] = int(parts[1])

# Function to get done files for a chromosome
def get_done_files_for_chrom(wildcards):
    return [f"singer/{c}/{c}_{s}_{e}.done" 
            for c, s, e in existing_chunks if c == wildcards.chr]

rule all:
    input:
        # ARG generation targets
        expand("singer/{chr}/{chr}_{index}.trees", chr=chrs, index=range(100)),
        # Per-chromosome metrics
        #expand("results/{chr}_metrics.csv", chr=chrs),
        expand("results/{chr}_metricsbyPop.csv", chr=chrs),
        
        # Genome-wide results
        #"results/genome-wide_metrics.csv",
        #"results/genome-wide_metricsbyPop.csv",
        
        # Annotated results
        #"results/genome-wide_metrics_annotated.csv",
        #"results/genome-wide_metrics_annotatedbyPop.csv",
        #"results/genome-wide_metrics_annotated_nonTelo.csv",
        #"results/genome-wide_metrics_annotatedbyPop_nonTelo.csv"

rule run_singer:
    input:
        vcf = "vcf_chunks/{chr}_{start}_{end}.vcf"  # Now has .vcf extension
    output:
        done = touch("singer/{chr}/{chr}_{start}_{end}.done"),
        trees = directory("singer/{chr}/{chr}_{start}_{end}_trees")
    params:
        output_prefix = "singer/{chr}/{chr}_{start}_{end}"
    shell:
        """
        mkdir -p singer/{wildcards.chr}
        mkdir -p {output.trees}
        ../singer_master \
            -Ne {NE} \
            -m {M} \
            -vcf {input.vcf} \
            -output {params.output_prefix} \
            -start {wildcards.start} \
            -end {wildcards.end}
        # Move output files to tree directory
        mv {params.output_prefix}_*.txt {output.trees}/
        """

rule merge_args:
    input:
        done_files = get_done_files_for_chrom,
    output:
        trees = expand("singer/{chr}/{chr}_{index}.trees",
                      chr="{chr}", index=range(100))
    params:
        file_list_dir = "file_lists/{chr}",
        missing_log = "logs/missing_files_{chr}.log"
    run:
        chrom = wildcards.chr
        os.makedirs(params.file_list_dir, exist_ok=True)
        
        # Get chunks for this chromosome
        chr_chunks = [(s, e) for c, s, e in existing_chunks if c == chrom]
        chr_chunks_sorted = sorted(chr_chunks, key=lambda x: int(x[0]))
        
        # Process each index (0-99)
        for index in range(100):
            file_list = []
            for s, e in chr_chunks_sorted:
                # Path to tree files in the chunk directory
                node_file = f"singer/{chrom}/{chrom}_{s}_{e}_trees/{chrom}_{s}_{e}_nodes_{index}.txt"
                branch_file = f"singer/{chrom}/{chrom}_{s}_{e}_trees/{chrom}_{s}_{e}_branches_{index}.txt"
                mut_file = f"singer/{chrom}/{chrom}_{s}_{e}_trees/{chrom}_{s}_{e}_muts_{index}.txt"
                
                if all(os.path.exists(f) for f in [node_file, branch_file, mut_file]):
                    file_list.append((node_file, branch_file, mut_file, s))
            
            if file_list:
                list_file = f"{params.file_list_dir}/{chrom}_{index}.txt"
                with open(list_file, "w") as f:
                    for node, branch, mut, s in sorted(file_list, key=lambda x: int(x[3])):
                        f.write(f"{node} {branch} {mut} {s}\n")
                
                # Create output directory if needed
                os.makedirs(os.path.dirname(output.trees[index]), exist_ok=True)
                
                shell(
                    f"python merge_ARGs.py "
                    f"--file_table {list_file} "
                    f"--chrom_sizes {CHROM_SIZES} "
                    f"--output {output.trees[index]}"
                )
            else:
                # Create empty output if no files found
                with open(output.trees[index], "w") as f:
                    f.write("")
                with open(params.missing_log, "a") as log:
                    log.write(f"No chunks found for {chrom} index {index}\n")

rule compute_tmrca_windows:
    input:
        trees = expand("singer/{chr}/{chr}_{index}.trees",
                      chr="{chr}",
                      index=range(NUM_ARGS+1)),
        chrom_sizes = CHROM_SIZES
    output:
        "results/{chr}_metrics.csv"
    log:
        "logs/{chr}_processing.log"
    params:
        chrom = "{chr}",
        num_args = NUM_ARGS
    shell:
        """
        python windowed_tmrca_stats.py \
            --chrom {wildcards.chr} \
            --chrom-sizes {input.chrom_sizes} \
            --output {output} \
            --num-args {params.num_args}
        """

rule compute_tmrca_windows_byPop:
    input:
        tree_dependencies = expand("singer/{chr}/{chr}_{index}.trees",
                                   chr="{chr}",
                                   index=range(NUM_ARGS+1)),
        chrom_sizes_file = CHROM_SIZES,
        pops_file = POPULATIONS,
        samples_file = SAMPLES
    output:
        outfile = "results/{chr}_metricsbyPop.csv"
    log:
        logfile = "logs/{chr}_processing_byPop.log" 
    params:
        chrom_val = "{chr}",  
        num_args_val = NUM_ARGS 
    shell:
        """
        python windowed_tmrca_stats_byPOP.py \
            --chrom {params.chrom_val} \
            --chrom-sizes {input.chrom_sizes_file} \
            --pops {input.pops_file} \
            --samples {input.samples_file} \
            --output {output.outfile} \
            --num-args {params.num_args_val} \
            > {log.logfile} 2>&1
        """


rule concatenate_metrics:
    input:
        expand("results/{chr}_metrics.csv", chr=chrs)
    output:
        "results/genome-wide_metrics.csv"
    run:
        pd.concat([pd.read_csv(f) for f in input]).to_csv(output[0], index=False)


rule concatenate_metrics_byPop:
    input:
        expand("results/{chr}_metricsbyPop.csv", chr=chrs)
    output:
        "results/genome-wide_metricsbyPop.csv"
    run:
        pd.concat([pd.read_csv(f) for f in input]).to_csv(output[0], index=False)


rule annotate_genes_byPop:
    input:
        metrics = "results/genome-wide_metricsbyPop.csv",
        genes = GENE_ANNOTATION
    output:
        "results/genome-wide_metrics_annotatedbyPop.csv"
    run:
        import os
        import pandas as pd
        import numpy as np
        
        # Temporary files
        temp_files = [
            "temp_metrics_pop.bed",
            "temp_metrics_pop_sorted.bed",
            "temp_genes_sorted.bed",
            "temp_annotated_pop.bed"
        ]

        try:
            # 1. Read metrics with proper NA handling
            df = pd.read_csv(input.metrics)
            
            # 2. Identify numeric columns
            numeric_cols = ['start', 'end', 'avg_tmrca', 'T_pooled', 'T_within', 'Tpooled_Twithin_ratio']
            
            # 3. Convert numeric columns, preserving NA values
            for col in numeric_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # 4. Write to BED format with NA as '.'
            df.to_csv(
                'temp_metrics_pop.bed',
                sep='\t',
                columns=['chromosome', 'start', 'end', 'population',
                        'avg_tmrca', 'T_pooled', 'T_within', 'Tpooled_Twithin_ratio'],
                header=False,
                index=False,
                na_rep='.'
            )
            
            # 5. Sort files
            if os.path.getsize('temp_metrics_pop.bed') == 0:
                raise ValueError("Empty BED file created")
                
            sort_cmd = "bedtools sort -i temp_metrics_pop.bed > temp_metrics_pop_sorted.bed"
            if os.system(sort_cmd) != 0:
                raise RuntimeError("bedtools sort failed")
            
            # 6. Sort gene annotations
            sort_genes_cmd = f"zcat {input.genes} | bedtools sort > temp_genes_sorted.bed"
            if os.system(sort_genes_cmd) != 0:
                raise RuntimeError("Gene sorting failed")
            
            # 7. Annotate with genes
            map_cmd = (
                "bedtools map "
                "-a temp_metrics_pop_sorted.bed "
                "-b temp_genes_sorted.bed "
                "-c 4 -o collapse > temp_annotated_pop.bed"
            )
            if os.system(map_cmd) != 0:
                raise RuntimeError("bedtools map failed")
            
            # 8. Create final output
            final_cols = [
                'chromosome', 'start', 'end', 'population',
                'avg_tmrca', 'T_pooled', 'T_within',
                'Tpooled_Twithin_ratio', 'genes'
            ]
            
            # Read annotated data
            annotated = pd.read_csv(
                'temp_annotated_pop.bed',
                sep='\t',
                header=None,
                names=final_cols[:-1] + ['genes'],
                dtype={'genes': str}
            )
            
            # Convert '.' back to empty strings
            annotated.replace('.', np.nan, inplace=True)
            
            # Ensure proper column types
            for col in numeric_cols:
                if col in annotated:
                    annotated[col] = pd.to_numeric(annotated[col], errors='coerce')
            
            # Write final output
            annotated[final_cols].to_csv(output[0], index=False)

        except Exception as e:
            # Clean up on error
            for f in temp_files:
                if os.path.exists(f):
                    os.remove(f)
            raise e
            
        finally:
            # Cleanup temporary files
            for f in temp_files:
                if os.path.exists(f):
                    os.remove(f)
