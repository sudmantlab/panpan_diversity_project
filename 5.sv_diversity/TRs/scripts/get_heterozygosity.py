#!/bin/env python3

# ---
# Author: Adam, Carolina de Lima
# Date: 2025
# Purpose: Calculate observed and expected heterozygosity based on TR allele lengths.
# Input:
    # A VCF file with: FORMAT field containing AL (allele length) and TDID (TRGT - Dolzhenko et al. 2024)
# Output:
    # A tab-separated file with columns: locus_id, sample, obs_het, exp_het
# ---

import sys
from pathlib import Path
from cyvcf2 import VCF
from collections import Counter

# check input
if len(sys.argv) != 2:
    print("Usage: python get_heterozygosity.py <input.vcf|input.vcf.gz|input.bcf>")
    sys.exit(1)

input_file = sys.argv[1]

# derive prefix
prefix = Path(input_file).name
prefix = prefix.replace(".vcf.gz", "").replace(".vcf", "").replace(".bcf", "")

# output files
per_sample_file = f"{prefix}_per_sample_het.txt"
per_locus_file = f"{prefix}_mean_het_per_locus.txt"

vcf = VCF(input_file)
samples = vcf.samples

out_sample = open(per_sample_file, "w")
out_locus = open(per_locus_file, "w")

out_sample.write("TRID\tindividual\tHo\tHe\n")
out_locus.write("TRID\tmean_Ho\tmean_He\n")

for variant in vcf:

    trid = variant.INFO.get("TRID")
    al = variant.format("AL")

    if al is None:
        continue

    Ho_list = []
    He_list = []

    for i, sample in enumerate(samples):
        alleles = al[i]

        if alleles is None or len(alleles) < 2:
            continue

        a1, a2 = alleles[0], alleles[1]

        if a1 < 0 or a2 < 0:
            continue

        # observed heterozygosity
        Ho = 1 if a1 != a2 else 0
        Ho_list.append(Ho)

        # expected heterozygosity
        length_counts = Counter([a1, a2])
        total = sum(length_counts.values())
        freqs = [count/total for count in length_counts.values()]
        He = 1 - sum(f**2 for f in freqs)
        He_list.append(He)

        out_sample.write(f"{trid}\t{sample}\t{Ho}\t{He}\n")

    if Ho_list:
        mean_Ho = sum(Ho_list) / len(Ho_list)
        mean_He = sum(He_list) / len(He_list)
        out_locus.write(f"{trid}\t{mean_Ho}\t{mean_He}\n")

out_sample.close()
out_locus.close()
