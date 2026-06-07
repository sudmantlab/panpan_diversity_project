from snakemake.io import glob_wildcards
import os
import pandas as pd
import pyranges as pr
## Before running it make sure the following exists:
#downloaded the chromosome sizes file
####wget https://github.com/marbl/T2T-Browser/raw/main/T2Tgenomes/mPanTro3_v2.0_pri/mPanTro3_v2.0_pri.sizes

# Converted to BED format (start=0, end=chromosome_length)
####awk '{print $1 "\t0\t" $2}' mPanTro3_v2.0_pri.sizes > chromosomes_full.bed

#downloaded the telomeres and centromeres
####bigBedToBed https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/pattern/mPanTro3_v2.0.GenomeFeature_v0.9.bb mPanTro3_genomefeatures.bed

#filter to get coordinates for primary hap (hap1)
####cat mPanTro3_genomefeatures.bed | grep 'hap1' > mPanTro3_hap1_genomefeatures.bed
####awk '$4 == "Cen"' mPanTro3_hap1_genomefeatures.bed > centromeres.bed

# Sort chromosome regions and centromere annotations
####sort -k1,1 -k2,2n chromosomes_full.bed > chromosomes_sorted.bed
####sort -k1,1 -k2,2n centromeres.bed > centromeres_sorted.bed 
# manually removed MT from chromosomes_sorted.bed 
# manually added hap2 Y to centromeres_sorted.bed in mPanTro3

#Split chromosomes into regions not overlapping centromeres
##### bedtools subtract \
#####    -a chromosomes_sorted.bed \
#####    -b centromeres_sorted.bed \
#####    > chromosomes_without_centromeres.bed

#### bedtools makewindows -b chromosomes_without_centromeres.bed -w 5000000 > 5mb_windows.bed

##### tabix -p vcf pantros_mapped2mPanTro3.sorted.filtered.vcf.gz 

# Split the VCF into chunks
##### mkdir -p vcf_chunks
##### while IFS=$'\t' read -r chr start end; do
#####  bcftools view \
#####    -r "${chr}:${start}-${end}" \
#####   -o "vcf_chunks/${chr}_${start}_${end}.vcf" \
#####    pantros_mapped2mPanTro3.sorted.filtered.vcf.gz
##### done < 5mb_windows.bed



# Configuration 
NE = "4e4"
M = "1.25e-8"
CHROM_SIZES = "mPanTro3_v2.0_pri.sizes"
POPULATIONS = "populations.csv"
SAMPLES = "samples.txt"
GENE_ANNOTATION = "mPanTro3_gene.bed.gz"

vcf_components = glob_wildcards(
    "vcf_chunks/{chr}_{start}_{end}.vcf",
    wildcard_constraints={  # ← CORRECT SPELLING
        "chr": r"[^_]+_.+",
        "start": r"\d+",
        "end": r"\d+"
    }
)

# Extract components
chrs = vcf_components.chr  
starts = vcf_components.start
ends = vcf_components.end
chunks = list(zip(chrs, starts, ends))

chrom_sizes = {}
with open(CHROM_SIZES) as f:
    for line in f:
        parts = line.strip().split()
        chrom_sizes[parts[0]] = int(parts[1])

# Validate chromosome consistency
vcf_chrs = set(chrs)
size_chrs = set(chrom_sizes.keys())
if vcf_chrs != size_chrs:
    print("ERROR: Chromosome mismatch between VCF and sizes file")
    print("Missing in sizes:", vcf_chrs - size_chrs)
    print("Missing in VCF:", size_chrs - vcf_chrs)
    raise SystemExit(1)

# Create ordered chromosome list from sizes file, filtered by VCF presence
ordered_chrs = [c for c in chrom_sizes if c in vcf_chrs]

rule all:
    input:
        expand("singer/{chr}/{chr}_{index}.trees", chr=set(chrs), index=range(100)),
        expand("results/{chr}_metricsbyPop.csv", chr=set(chrs)),
        expand("results/{chr}_metrics.csv", chr=set(chrs)),
        "results/genome-wide_metrics_annotatedbyPop.csv",
        "results/genome-wide_metrics_annotated.csv"

# Remaining rules identical to your last version (run_singer through annotate_genes)

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
    output:
        expand("singer/{chr}/{chr}_{index}.trees", index=range(100), allow_missing=True)
    params:
        file_list_dir="file_lists",
        missing_log="logs/missing_files.log"
    run:
        chrom = wildcards.chr
        os.makedirs(params.file_list_dir, exist_ok=True)

        chrom_length = chrom_sizes[chrom]
        chrom_blocks = [(s, e) for c, s, e in chunks if c == chrom]
        chrom_blocks_sorted = sorted(chrom_blocks, key=lambda x: int(x[0]))

        for index in range(100):
            file_list = []
            for s, e in chrom_blocks_sorted:
                node_file = f"singer/{chrom}/{chrom}_{s}_{e}_nodes_{index}.txt"
                branch_file = f"singer/{chrom}/{chrom}_{s}_{e}_branches_{index}.txt"
                mut_file = f"singer/{chrom}/{chrom}_{s}_{e}_muts_{index}.txt"

                if all(os.path.exists(f) for f in [node_file, branch_file, mut_file]):
                    file_list.append((node_file, branch_file, mut_file, s))

            if file_list:
                list_path = f"{params.file_list_dir}/{chrom}_{index}.txt"
                with open(list_path, "w") as f:
                    for node, branch, mut, s in sorted(file_list, key=lambda x: int(x[3])):
                        f.write(f"{node} {branch} {mut} {s}\n")

                shell(
                    f"python merge_ARGs.py "
                    f"--file_table {list_path} "
                    f"--chrom_sizes {CHROM_SIZES} "
                    f"--output singer/{chrom}/{chrom}_{index}.trees"
                )


rule compute_tmrca_windows_byPop:
    input:
        trees = expand("singer/{chr}/{chr}_{index}.trees", index=range(100)),
        chrom_sizes = CHROM_SIZES,
        pops = POPULATIONS,
        samples = SAMPLES
    output:
        "results/{chr}_metricsbyPop.csv"
    params:
        chrom = "{chr}"
    script:
        "windowed_tmrca_stats_byPOP.py"


rule compute_tmrca_windows:
    input:
        trees = expand("singer/{chr}/{chr}_{index}.trees", index=range(100)),
        chrom_sizes = CHROM_SIZES,
    output:
        "results/{chr}_metrics.csv"
    params:
        chrom = "{chr}"
    script:
        "windowed_tmrca_stats.py"

rule concatenate_metrics_byPop:
    input:
        expand("results/{chr}_metricsbyPop.csv", chr=ordered_chrs)  # ← uses ordered list
    output:
        "results/genome-wide_metricsbyPop.csv"
    run:
        import pandas as pd
        pd.concat([pd.read_csv(f) for f in input], ignore_index=True).to_csv(output[0], index=False)

rule concatenate_metrics:
    input:
        expand("results/{chr}_metrics.csv", chr=ordered_chrs)
    output:
        "results/genome-wide_metrics.csv"
    run:
        pd.concat([pd.read_csv(f) for f in input], ignore_index=True).to_csv(output[0], index=False)


rule annotate_genes_byPop:
    input:
        metrics = "results/genome-wide_metricsbyPop.csv",
        genes = GENE_ANNOTATION
    output:
        "results/genome-wide_metrics_annotatedbyPop.csv"
    run:
        metrics_df = pd.read_csv(input.metrics)
        genes_df = pd.read_csv(
            input.genes,
            sep='\t',
            header=None,
            names=["Chromosome", "Start", "End", "Gene"]
        )

        metrics_pr = pr.PyRanges(
            metrics_df.rename(columns={
                "chromosome": "Chromosome",
                "start": "Start",
                "end": "End"
            })
        )
        genes_pr = pr.PyRanges(genes_df)

        overlaps = metrics_pr.join(genes_pr, how="left", strandedness=False).df

        annotated_df = overlaps.groupby([
            'Chromosome', 'Start', 'End', 
            'population', 'avg_pairwise_coalescence_time', 
            'avg_tmrca', 'Tpooled_Twithin_ratio'
        ])['Gene'].agg(lambda x: ','.join(x.dropna())).reset_index()

        annotated_df = annotated_df[[
            'Chromosome', 'Start', 'End', 'population',
            'avg_pairwise_coalescence_time', 'avg_tmrca',
            'Tpooled_Twithin_ratio', 'Gene'
        ]].rename(columns={
            "Chromosome": "chromosome",
            "Start": "start",
            "End": "end",
            "Gene": "genes"
        })

        annotated_df.to_csv(output[0], index=False)

# Non-population annotation
rule annotate_genes:
    input:
        metrics = "results/genome-wide_metrics.csv",
        genes = GENE_ANNOTATION
    output:
        "results/genome-wide_metrics_annotated.csv"
    run:
        metrics_df = pd.read_csv(input.metrics)
        genes_df = pd.read_csv(input.genes, sep='\t', header=None, 
                             names=["Chromosome", "Start", "End", "Gene"])

        metrics_pr = pr.PyRanges(
            metrics_df.rename(columns={
                "chromosome": "Chromosome",
                "start": "Start",
                "end": "End"
            })
        )
        genes_pr = pr.PyRanges(genes_df)

        overlaps = metrics_pr.join(genes_pr, how="left", strandedness=False).df

        # Non-pop grouping
        annotated_df = overlaps.groupby([
            'Chromosome', 'Start', 'End',
            'avg_pairwise_coalescence_time', 'avg_tmrca'
        ])['Gene'].agg(lambda x: ','.join(x.dropna())).reset_index()

        annotated_df = annotated_df[[
            'Chromosome', 'Start', 'End',
            'avg_pairwise_coalescence_time', 'avg_tmrca', 'Gene'
        ]].rename(columns={
            "Chromosome": "chromosome",
            "Start": "start",
            "End": "end",
            "Gene": "genes"
        })

        annotated_df.to_csv(output[0], index=False)
