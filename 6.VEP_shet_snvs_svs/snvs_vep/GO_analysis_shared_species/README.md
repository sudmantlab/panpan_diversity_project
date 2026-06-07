# GO / pathway / disease enrichment — SNV-impacted genes shared across great apes

ShinyGO v0.85 (Geet al., 2020) is run using the same gene lists and the SNV-tested-universe background.


### LOC → HGNC renaming

The raw SNV-tested universe contains **15,465 NCBI provisional `LOC<number>`
entries (≈39% of the 39,882 rows)** — a much larger LOC share than the SV
pipeline. ShinyGO / enrichR / KEGG generally do not recognise LOC numbers. 
The `scripts/_build_loc_map.R` helper queries the local NCBI `Homo_sapiens.gene_info` table and tries to match each LOC's
description against three NCBI fields:

1. `description` (e.g. "alpha-amylase 1B" → AMY1B)
2. `Full_name_from_nomenclature_authority`
3. `Other_designations` (pipe-separated synonyms)

after also stripping `-like` and `protein …` suffixes.

Result: **1,350 of 15,465 LOC entries (8.7%) cleanly renamed** to their
NCBI/HGNC-approved Symbol. The lower mapping rate vs. the SV pipeline
(13.7%) reflects the larger long tail of uncharacterised loci in the SNV
universe. Mapping is cached at `backgrounds/loc_to_hgnc_map.tsv`. The
remaining 14,115 LOC entries are genuinely uncharacterized loci with no
NCBI counterpart — they stay in the universe under their `LOC<number>` name
so the background size still reflects the count of testable genes, even
though ShinyGO will silently ignore them.