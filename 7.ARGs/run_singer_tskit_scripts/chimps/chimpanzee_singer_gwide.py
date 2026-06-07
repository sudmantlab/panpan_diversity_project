from snakemake.io import glob_wildcards
import os
import pandas as pd

# --- PARAMETERS ---
NE = "4e4"
M = "1.25e-8"
CHROM_SIZES = "mPanTro3_v2.0_pri.sizes"
POPULATIONS = "populations.csv"
SAMPLES = "samples.txt"
GENE_ANNOTATION = "mPanTro3_gene.bed.gz"
NUM_ARGS = 99  # For 100 samples (0-99)

# --- SETUP ---

# Find all VCF chunks to determine the chromosomes to process
try:
    vcf_components = glob_wildcards("vcf_chunks/{chr}_{start}_{end}.vcf")
    existing_chunks = list(zip(vcf_components.chr, vcf_components.start, vcf_components.end))
    chrs = sorted(list(set(vcf_components.chr)))
except FileNotFoundError:
    chrs = []
    existing_chunks = []
    print("Warning: 'vcf_chunks' directory not found or empty. SINGER pipeline may not run.")

# Read regions of interest for regional tree extraction
try:
    regions_df = pd.read_csv("regions.bed", sep='\\s+', header=None,
                             names=["chr", "start", "end", "name", "peak_pos"]) # Reads 5th column
except FileNotFoundError:
    print("Warning: 'regions.bed' not found. Rules for regional extraction will not be active.")
    regions_df = pd.DataFrame(columns=["chr", "start", "end", "name", "peak_pos"])

# Define the offsets
OFFSETS = [-5, 0, 5]


all_svg_files = []
for index in range(NUM_ARGS + 1):
    for offset in OFFSETS:
        # Add offset suffix to filename (e.g., _minus5, _0, _plus5)
        offset_str = f"_minus{abs(offset)}" if offset < 0 else f"_plus{offset}" if offset > 0 else "_peak"
        all_svg_files.extend(
            expand(
                # Add offset_str to path/filename
                "regional_svgs/{chr}/{start}-{end}_{name}/{peak_pos}/{offset_str}/{chr}_{start}-{end}_{name}_{peak_pos}_{offset_str}_{index}.svg",
                zip,
                chr=regions_df.chr,
                start=regions_df.start,
                end=regions_df.end,
                name=regions_df.name,
                peak_pos=regions_df.peak_pos,
                offset_str=offset_str, # Use the generated string
                index=index
            )
        )




# --- all_nwk_files list generation ---
all_nwk_files = []
for index in range(NUM_ARGS + 1):
    all_nwk_files.extend(
        expand(
            # ADDED {peak_pos} to filename
            "newick_trees/{chr}/{start}-{end}_{name}/{peak_pos}/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.nwk",
            zip,
            chr=regions_df.chr,
            start=regions_df.start,
            end=regions_df.end,
            name=regions_df.name,
            peak_pos=regions_df.peak_pos, # Expand using peak_pos column
            index=index
        )
    )



final_nwk_files = []
for index in range(NUM_ARGS + 1):
    final_nwk_files.extend(
        expand(
            # Path matches the output of relabel_rescale_newick
            "final_trees/{chr}/{start}-{end}_{name}/{peak_pos}/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.nwk",
            zip,
            chr=regions_df.chr,
            start=regions_df.start,
            end=regions_df.end,
            name=regions_df.name,
            peak_pos=regions_df.peak_pos,
            index=index
        )
    )


# Helper function for merge_args rule
def get_done_files_for_chrom(wildcards):
    return [f"singer/{c}/{c}_{s}_{e}.done"
            for c, s, e in existing_chunks if c == wildcards.chr]


# --- FINAL TARGETS ---
rule all:
    input:
        # Your original targets
        #expand("singer/{chr}/{chr}_{index}.trees", chr=chrs, index=range(NUM_ARGS + 1)),
        #expand("results/{chr}_metricsbyPop.csv", chr=chrs),
        #"results/genome-wide_metrics_annotatedbyPop.csv",
        # New target for regional .nwk files
        all_svg_files,
        #all_nwk_files,
        #final_nwk_files
        

      

# --- SINGER & TMRCA PIPELINE RULES ---

rule run_singer:
    input:
        "vcf_chunks/{chr}_{start}_{end}.vcf"
    output:
        done=touch("singer/{chr}/{chr}_{start}_{end}.done"),
        trees=directory("singer/{chr}/{chr}_{start}_{end}_trees")
    params:
        vcf_prefix="vcf_chunks/{chr}_{start}_{end}",
        output_prefix="singer/{chr}/{chr}_{start}_{end}"
    shell:
        """
        mkdir -p singer/{wildcards.chr} {output.trees}
        ../singer_master \\
            -Ne {NE} \\
            -m {M} \\
            -vcf {params.vcf_prefix} \\
            -output {params.output_prefix} \\
            -start {wildcards.start} \\
            -end {wildcards.end}

        mv {params.output_prefix}_*.txt {output.trees}/
        """

rule merge_args:
    input:
        done_files=get_done_files_for_chrom,
    output:
        trees=expand("singer/{chr}/{chr}_{index}.trees",
                      chr="{chr}", index=range(NUM_ARGS + 1))
    params:
        file_list_dir="file_lists/{chr}",
        missing_log="logs/missing_files_{chr}.log"
    run:
        chrom = wildcards.chr
        os.makedirs(params.file_list_dir, exist_ok=True)

        chr_chunks = [(s, e) for c, s, e in existing_chunks if c == chrom]
        chr_chunks_sorted = sorted(chr_chunks, key=lambda x: int(x[0]))

        for index in range(NUM_ARGS + 1):
            file_list = []
            for s, e in chr_chunks_sorted:
                node_file = f"singer/{chrom}/{chrom}_{s}_{e}_trees/{chrom}_{s}_{e}_nodes_{index}.txt"
                branch_file = f"singer/{chrom}/{chrom}_{s}_{e}_trees/{chrom}_{s}_{e}_branches_{index}.txt"
                mut_file = f"singer/{chrom}/{chrom}_{s}_{e}_trees/{chrom}_{s}_{e}_muts_{index}.txt"

                if all(os.path.exists(f) for f in [node_file, branch_file, mut_file]):
                    file_list.append((node_file, branch_file, mut_file, s))

            if file_list:
                list_file = f"{params.file_list_dir}/{chrom}_{index}.txt"
                with open(list_file, "w") as f:
                    for node, branch, mut, s in sorted(file_list, key=lambda x: int(x[3])):
                        f.write(f"{node} {branch} {mut} {s}\\n")

                os.makedirs(os.path.dirname(output.trees[index]), exist_ok=True)

                shell(
                    "python merge_ARGs.py "
                    "--file_table {list_file} "
                    "--chrom_sizes {CHROM_SIZES} "
                    "--output {output.trees[index]}"
                )
            else:
                # --- THIS BLOCK IS NOW CORRECTED ---
                # Create an empty output file
                with open(output.trees[index], "w") as f:
                    pass
                # Log that no chunks were found
                with open(params.missing_log, "a") as log:
                    log.write(f"No chunks found for {chrom} index {index}\\n")

rule compute_tmrca_windows_byPop:
    input:
        tree_dependencies=expand("singer/{chr}/{chr}_{index}.trees",
                                   chr="{chr}",
                                   index=range(NUM_ARGS + 1)),
        chrom_sizes_file=CHROM_SIZES,
        pops_file=POPULATIONS,
        samples_file=SAMPLES
    output:
        outfile="results/{chr}_metricsbyPop.csv",
        diag_dir=directory("pop_void_diagnostics/{chr}")
    log:
        logfile="logs/{chr}_processing_byPop.log"
    params:
        chrom_val="{chr}",
        num_args_val=NUM_ARGS
    shell:
        """
        python windowed_tmrca_stats_byPOP.py \\
            --chrom {params.chrom_val} \\
            --chrom-sizes {input.chrom_sizes_file} \\
            --pops {input.pops_file} \\
            --samples {input.samples_file} \\
            --output {output.outfile} \\
            --diag-output {output.diag_dir} \\
            --num-args {params.num_args_val} > {log.logfile} 2>&1
        """

rule concatenate_metrics_byPop:
    input:
        expand("results/{chr}_metricsbyPop.csv", chr=chrs)
    output:
        "results/genome-wide_metricsbyPop.csv"
    run:
        pd.concat([pd.read_csv(f) for f in input]).to_csv(output[0], index=False)

rule annotate_genes_byPop:
    input:
        metrics="results/genome-wide_metricsbyPop.csv",
        genes=GENE_ANNOTATION
    output:
        "results/genome-wide_metrics_annotatedbyPop.csv"
    run:
        # This is a complex run block, leaving your original implementation as is.
        # It's recommended to move complex logic to external scripts for clarity.
        import os
        import pandas as pd
        import numpy as np
        temp_files = ["temp_metrics_pop.bed", "temp_metrics_pop_sorted.bed", "temp_genes_sorted.bed", "temp_annotated_pop.bed"]
        try:
            df = pd.read_csv(input.metrics)
            numeric_cols = ['start', 'end', 'avg_tmrca', 'T_pooled', 'T_within', 'Tpooled_Twithin_ratio']
            for col in numeric_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            df.to_csv('temp_metrics_pop.bed', sep='\\t', columns=['chromosome', 'start', 'end', 'population', 'avg_tmrca', 'T_pooled', 'T_within', 'Tpooled_Twithin_ratio'], header=False, index=False, na_rep='.')
            if os.path.getsize('temp_metrics_pop.bed') == 0: raise ValueError("Empty BED file created")
            if os.system("bedtools sort -i temp_metrics_pop.bed > temp_metrics_pop_sorted.bed") != 0: raise RuntimeError("bedtools sort failed")
            if os.system(f"zcat {input.genes} | bedtools sort > temp_genes_sorted.bed") != 0: raise RuntimeError("Gene sorting failed")
            if os.system("bedtools map -a temp_metrics_pop_sorted.bed -b temp_genes_sorted.bed -c 4 -o collapse > temp_annotated_pop.bed") != 0: raise RuntimeError("bedtools map failed")
            final_cols = ['chromosome', 'start', 'end', 'population', 'avg_tmrca', 'T_pooled', 'T_within', 'Tpooled_Twithin_ratio', 'genes']
            annotated = pd.read_csv('temp_annotated_pop.bed', sep='\\t', header=None, names=final_cols[:-1] + ['genes'], dtype={'genes': str})
            annotated.replace('.', np.nan, inplace=True)
            for col in numeric_cols:
                if col in annotated: annotated[col] = pd.to_numeric(annotated[col], errors='coerce')
            annotated[final_cols].to_csv(output[0], index=False)
        finally:
            for f in temp_files:
                if os.path.exists(f): os.remove(f)


# --- REGIONAL EXTRACTION AND CONVERSION RULES ---

rule extract_regional_trees:
    input:
        full_trees="singer/{chr}/{chr}_{index}.trees"
    output:
        # NOTE: Output filename doesn't need peak_pos here, as it processes the whole region
        regional_trees="regional_trees/{chr}/{start}-{end}_{name}/{chr}_{start}-{end}_{name}_{index}.trees"
    params:
        start="{start}",
        end="{end}"
    log:
        "logs/extract_regional/{chr}_{start}-{end}_{name}_{index}.log"
    wildcard_constraints:
        start=r"\d+",
        end=r"\d+"
    shell:
        """
        python extract_regional_trees.py \\
            --trees-in {input.full_trees} \\
            --start {params.start} \\
            --end {params.end} \\
            --trees-out {output.regional_trees} > {log} 2>&1
        """

# --- Modify the convert_to_newick rule ---
rule convert_to_newick:
    input:
        regional_trees="regional_trees/{chr}/{start}-{end}_{name}/{chr}_{start}-{end}_{name}_{index}.trees"
    output:
        newick="newick_trees/{chr}/{start}-{end}_{name}/{peak_pos}/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.nwk"
    log:
        "logs/newick_conversion/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.log"
    params:
        region_start="{start}",
        absolute_pos="{peak_pos}"
    wildcard_constraints:
        start=r"\d+",
        end=r"\d+",
        peak_pos=r"\d+"
    shell:
        """
        python ./ts_to_newick.py \\
            --ts-file {input.regional_trees} \\
            --newick-file {output.newick} \\
            --absolute-position {params.absolute_pos} \\
            --region-start {params.region_start} > {log} 2>&1
        """


rule relabel_rescale_newick:
    input:
        nwk_in="newick_trees/{chr}/{start}-{end}_{name}/{peak_pos}/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.nwk",
        samples=SAMPLES,
        populations=POPULATIONS
    output:
        nwk_out="final_trees/{chr}/{start}-{end}_{name}/{peak_pos}/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.nwk"
    log:
        "logs/relabel_rescale/{chr}_{start}-{end}_{name}_{peak_pos}_{index}.log"
    params:
        gen_time=25
    wildcard_constraints:
        start=r"\d+",
        end=r"\d+",
        peak_pos=r"\d+"
    shell:
        """
        python ./relabel_rescale_nwk.py --nwk-in {input.nwk_in} --nwk-out {output.nwk_out} --samples-file {input.samples} --pop-file {input.populations} --gen-time {params.gen_time} > {log} 2>&1
        """

# --- Modify the plot_regional_svg_labeled rule ---
rule plot_regional_svg_labeled:
    input:
        regional_trees="regional_trees/{chr}/{start}-{end}_{name}/{chr}_{start}-{end}_{name}_{index}.trees",
        samples=SAMPLES,
        populations=POPULATIONS
    output:
        # Filename now includes {offset_str} wildcard
        svg="regional_svgs/{chr}/{start}-{end}_{name}/{peak_pos}/{offset_str}/{chr}_{start}-{end}_{name}_{peak_pos}_{offset_str}_{index}.svg"
    log:
        "logs/svg_plotting_labeled/{chr}_{start}-{end}_{name}_{peak_pos}_{offset_str}_{index}.log"
    params:
        region_start="{start}",
        absolute_pos="{peak_pos}",
        # Determine offset based on {offset_str} wildcard
        offset=lambda wildcards: 0 if wildcards.offset_str == "_peak" else -5 if wildcards.offset_str == "_minus5" else 5 if wildcards.offset_str == "_plus5" else 0
    wildcard_constraints:
        start=r"\d+",
        end=r"\d+",
        peak_pos=r"\d+",
        offset_str=r"_peak|_minus5|_plus5" # Constraint for the new wildcard
    shell:
        """
        python ./plot_regional_svg_labeled.py \\
            --ts-file {input.regional_trees} \\
            --svg-file {output.svg} \\
            --samples-file {input.samples} \\
            --pop-file {input.populations} \\
            --absolute-position {params.absolute_pos} \\
            --region-start {params.region_start} \\
            --offset {params.offset} > {log} 2>&1 # Pass the offset
        """