# Snakefile

from snakemake.io import glob_wildcards
import os

# Configuration
NE = "4e4"
M = "1.25e-8"
CHROM_SIZES = "mPanPan1_v2.0_pri.sizes"  # Path to your sizes file

# Get VCF chunks and chromosome sizes
vcf_components = glob_wildcards("vcf_chunks/{chr}_{start}_{end}.vcf")
chrs = vcf_components.chr
starts = vcf_components.start
ends = vcf_components.end
chunks = list(zip(chrs, starts, ends))

# Load chromosome sizes into dictionary
chrom_sizes = {}
with open(CHROM_SIZES) as f:
    for line in f:
        chrom, size = line.strip().split()
        chrom_sizes[chrom] = int(size)

rule all:
    input:
        expand("singer/{chr}/{chr}_{index}.trees", chr=set(chrs), index=range(100))

rule run_singer:
    input:
        "vcf_chunks/{chr}_{start}_{end}.vcf"
    output:
        touch("singer/{chr}/{chr}_{start}_{end}.done")
    params:
        vcf_prefix = "vcf_chunks/{chr}_{start}_{end}",
        output_prefix = "singer/{chr}/{chr}_{start}_{end}"
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
        done_files = expand("singer/{chr}/{chr}_{start}_{end}.done", zip, chr=chrs, start=starts, end=ends),
        chrom_sizes = CHROM_SIZES
    output:
        trees = expand("singer/{chr}/{chr}_{index}.trees", index=range(100), allow_missing=True)
    params:
        file_list_dir = "file_lists",
        missing_log = "logs/missing_files.log"
    run:
        chrom = wildcards.chr
        os.makedirs(params.file_list_dir, exist_ok=True)
        os.makedirs(os.path.dirname(params.missing_log), exist_ok=True)

        # Get chromosome length
        chrom_length = chrom_sizes.get(chrom, 0)
        if chrom_length == 0:
            raise ValueError(f"Chromosome {chrom} not found in {CHROM_SIZES}")

        # Generate file lists for each index (0-99)
        for index in range(100):
            file_list = []
            for s, e in [ (s, e) for c, s, e in chunks if c == chrom ]:
                node_file = f"singer/{chrom}/{chrom}_{s}_{e}_nodes_{index}.txt"
                branch_file = f"singer/{chrom}/{chrom}_{s}_{e}_branches_{index}.txt"
                mut_file = f"singer/{chrom}/{chrom}_{s}_{e}_muts_{index}.txt"

                if os.path.exists(node_file) and os.path.exists(branch_file) and os.path.exists(mut_file):
                    file_list.append((node_file, branch_file, mut_file, s))

            if file_list:
                # Write file list
                list_path = f"{params.file_list_dir}/{chrom}_{index}.txt"
                with open(list_path, "w") as f:
                    for node, branch, mut, s in file_list:
                        f.write(f"{node} {branch} {mut} {s}\n")

                # Run merge script
                output_file = f"singer/{chrom}/{chrom}_{index}.trees"
                shell(
                    f"python merge_ARGs.py "
                    f"--file_table {list_path} "
                    f"--chrom_sizes {CHROM_SIZES} "
                    f"--output {output_file}"
                )