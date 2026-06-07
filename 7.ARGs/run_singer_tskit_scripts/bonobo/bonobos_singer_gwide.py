from snakemake.io import glob_wildcards
import os
import pandas as pd
import pyranges as pr


NE = "4e4"
M = "1.25e-8"
CHROM_SIZES = "mPanPan1_v2.0_pri.sizes"
POPULATIONS = "populations.csv"
SAMPLES = "samples.txt"
GENE_ANNOTATION = "mPanPan1_gene.bed.gz"
GENOME_FEATURES = "./mPanPan1_genomefeatures.bed"
NUM_ARGS = 99  # For 100 samples (0-99)


vcf_components = glob_wildcards("vcf_chunks/{chr}_{start}_{end}.vcf")
chrs = sorted(list(set(vcf_components.chr))) 

chrom_sizes = {}
with open(CHROM_SIZES) as f:
    for line in f:
        parts = line.strip().split()
        chrom_sizes[parts[0]] = int(parts[1])

rule all:
    input:
        # ARG generation targets
        expand("singer/{chr}/{chr}_{index}.trees", chr=chrs, index=range(100)),
        
        # Per-chromosome metrics
        expand("results/{chr}_metrics.csv", chr=chrs),
        expand("void_diagnostics/{chr}", chr=chrs),  
        #expand("results/{chr}_metricsbyPop.csv", chr=chrs),
        
        # Genome-wide results
        "results/genome-wide_metrics.csv",
        
        # Annotated results
        "results/genome-wide_metrics_annotated.csv",

rule run_singer:
    input:
        "vcf_chunks/{chr}_{start}_{end}.vcf"
    output:
        touch("singer/{chr}/{chr}_{start}_{end}.done")
    params:
        vcf_prefix="vcf_chunks/{chr}_{start}_{end}",
        output_prefix="singer/{chr}/{chr}_{start}_{end}"
    shell:
        """
        mkdir -p singer/{wildcards.chr}
        ./singer_master \
            -Ne {NE} \
            -m {M} \
            -vcf {params.vcf_prefix} \
            -output {params.output_prefix} \
            -start {wildcards.start} \
            -end {wildcards.end}
        """

rule merge_args:
    input:
        done_files = expand("singer/{chr}/{chr}_{start}_{end}.done",
                          chr=chrs,
                          start=vcf_components.start,
                          end=vcf_components.end)
    output:
        trees = expand("singer/{chr}/{chr}_{index}.trees",
                      chr=chrs,
                      index=range(100))
    params:
        file_list_dir = "file_lists",
        missing_log = "logs/missing_files.log"
    run:
        chrom = wildcards.chr
        os.makedirs(params.file_list_dir, exist_ok=True)
        chr_blocks = [(s, e) for c, s, e in zip(vcf_components.chr,
                                                vcf_components.start,
                                                vcf_components.end)
                     if c == chrom]
        chr_blocks_sorted = sorted(chr_blocks, key=lambda x: int(x[0]))

        for index in range(100):
            file_list = []
            for s, e in chr_blocks_sorted:
                node_file = f"singer/{chrom}/{chrom}_{s}_{e}_nodes_{index}.txt"
                branch_file = f"singer/{chrom}/{chrom}_{s}_{e}_branches_{index}.txt"
                mut_file = f"singer/{chrom}/{chrom}_{s}_{e}_muts_{index}.txt"

                if all(os.path.exists(f) for f in [node_file, branch_file, mut_file]):
                    file_list.append((node_file, branch_file, mut_file, s))

            if file_list:
                list_file = f"{params.file_list_dir}/{chrom}_{index}.txt"
                with open(list_file, "w") as f:
                    for node, branch, mut, s in sorted(file_list, key=lambda x: int(x[3])):
                        f.write(f"{node} {branch} {mut} {s}\n")

                shell(
                    f"python merge_ARGs.py "
                    f"--file_table {list_file} "
                    f"--chrom_sizes {CHROM_SIZES} "
                    f"--output singer/{chrom}/{chrom}_{index}.trees"
                )

rule compute_tmrca_windows:
    input:
        trees = expand("singer/{chr}/{chr}_{index}.trees",
                      chr="{chr}",
                      index=range(NUM_ARGS+1)),
        chrom_sizes = CHROM_SIZES
    output:
        metrics="results/{chr}_metrics.csv",
        diag_dir=directory("void_diagnostics/{chr}")
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
            --output {output.metrics} \
            --num-args {params.num_args} \
            --diag-output {output.diag_dir}
        """

rule concatenate_metrics:
    input:
        expand("results/{chr}_metrics.csv", chr=chrs)
    output:
        "results/genome-wide_metrics.csv"
    run:
        pd.concat([pd.read_csv(f) for f in input]).to_csv(output[0], index=False)


rule annotate_genes:
    input:
        metrics = "results/genome-wide_metrics.csv",
        genes = GENE_ANNOTATION
    output:
        "results/genome-wide_metrics_annotated.csv"
    run:
        import os
        
        temp_files = [
            "temp_metrics.bed",
            "temp_metrics_sorted.bed",
            "temp_genes_sorted.bed",
            "temp_annotated.bed"
        ]

        try:
            # Convert CSV to sorted BED
            shell(r"""
            awk -F',' 'BEGIN {{OFS="\t"}} 
            NR==1 {{next}}  # Skip header
            {{print $1, $2, $3, $4, $5}}' {input.metrics} > temp_metrics.bed
            """)
            
            shell("bedtools sort -i temp_metrics.bed > temp_metrics_sorted.bed")
            
            # Process gene annotations
            shell("zcat {input.genes} | bedtools sort > temp_genes_sorted.bed")
            
            # Map genes to metrics
            shell(
                "bedtools map -a temp_metrics_sorted.bed "
                "-b temp_genes_sorted.bed -c 4 -o collapse > temp_annotated.bed"
            )

            # Final CSV formatting with proper AWK syntax
            shell(r"""
            echo "chromosome,start,end,avg_tmrca,avg_pairwise_coalescence_time,genes" > {output}
            awk -F'\t' 'BEGIN {{OFS=","}} 
            {{print $1, $2, $3, $4, $5, $6}}' temp_annotated.bed >> {output}
            """)

        finally:
            # Cleanup temporary files
            for f in temp_files:
                if os.path.exists(f):
                    os.remove(f)
