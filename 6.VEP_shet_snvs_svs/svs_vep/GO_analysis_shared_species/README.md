# GO / pathway / disease enrichment — SV-impacted genes shared across great apes

Companion analysis for Figure 3 (PANPAN manuscript). Tests whether the gene
sets shared across humans, chimpanzees and bonobos that are high-impact
targets of structural variants are enriched in specific biological processes,
molecular functions, pathways, or human-disease terms.

The pipeline is **enrichR-only** (`run_GO_enrichment.R`); ShinyGO v0.85 (Ge
et al., 2020) is run interactively on the web as an orthogonal cross-check
using the same gene lists and an Ensembl protein-coding background that the
script writes into each per-list folder.

## Gene lists (`gene_lists/`)

| File | n | Definition |
|------|---|------------|
| `cross_species_shared_60.txt` | 60 | All rows of `TableSX_genes_impacted_by_svs_shared_across_species.csv` (genes with impact = 4 in ≥ 2 species). |
| `human_chimp_shared_35.txt`   | 35 | Subset of the above with `human = 4 AND chimpanzee = 4 AND bonobo ≠ 4` (patterns 4,4,0 and 4,4,1). |

Gene IDs are HGNC symbols. The minimal LOC-RNA filter from
[`vep_shet_SVs.R:14-16`](../vep_shet_SVs.R) is applied (drops LOC*
spliceosomal-RNA and LOC* small-RNA entries). No stricter protein-coding
filter is applied because well-known polymorphic genes (e.g. GSTM1) are
missing from the T2T BED and we want to keep them.

> Per-species lists (`human_high_impact`, `chimp_high_impact`, `bonobo_high_impact`)
> are still computed inside the script but excluded from `gene_lists` for the
> current focus. Add them back to the `gene_lists` list near the top of
> `scripts/run_GO_enrichment.R` to re-enable.

## Backgrounds

**Recommended ShinyGO background = the tested universe** — the genes
actually evaluated for SV impact (`gene_impact_wide.tsv`), after LOC-RNA
cleaning **and** an NCBI-based LOC→HGNC renaming pass (see below). This is
the statistically appropriate background because it controls for which
genes the SV-impact pipeline could test in the first place — a whole-genome
background would inflate enrichment by including genes the analysis never
had a chance to call.

### LOC → HGNC renaming

The raw tested universe contained 3,377 NCBI provisional `LOC<number>`
entries (20% of the 16,513 rows). ShinyGO / enrichR / KEGG generally do
not recognise LOC numbers, so they're dead weight in the test. The
`scripts/_build_loc_map.R` helper queries the local NCBI `Homo_sapiens.gene_info`
table (cached at `backgrounds/Homo_sapiens.gene_info.gz`) and tries to match
each LOC's description against three NCBI fields:

1. `description` (e.g. "alpha-amylase 1B" → AMY1B)
2. `Full_name_from_nomenclature_authority`
3. `Other_designations` (pipe-separated synonyms)

after also stripping `-like` and `protein …` suffixes.

Result: **463 of 3,377 LOC entries (13.7%) cleanly renamed** to their
NCBI/HGNC-approved Symbol — including genes like ABO (LOC450164), AMY1B,
FCGR3B (LOC124905743), KLRG1, PTMA, CRLF2, EIF3F, multiple ZNF and KRTAP
paralogs. Mapping is cached at `backgrounds/loc_to_hgnc_map.tsv`. The
remaining 2,914 LOC entries are genuinely uncharacterized loci with no NCBI
counterpart — they stay in the universe under their `LOC<number>` name so
the background size still reflects the count of testable genes, even though
ShinyGO will silently ignore them.

### Files

| File | n | Purpose |
|------|---|---------|
| `backgrounds/Homo_sapiens.gene_info.gz` | — | NCBI gene_info, source for the LOC renaming. |
| `backgrounds/loc_to_hgnc_map.tsv` | 3,377 | Per-LOC mapping (`loc_name`, `new_symbol`, `description`). |
| `backgrounds/tested_universe.txt` | **16,116** | Primary background — VEP-tested genes, LOC-renamed and deduplicated. |
| `backgrounds/ensembl_human_protein_coding.txt` | 19,468 | Alternative whole-genome bg (Ensembl BioMart). Use only if cross-comparing with a wider universe. |
| `<list>/shinygo_background_tested.txt` | 16,116 | Drop into ShinyGO's "Custom" background field — primary. |
| `<list>/shinygo_background_ensembl.txt` | 19,468 | Drop into ShinyGO's "Custom" background field — alternative. |

> Why the count dropped 16,513 → 16,116: 463 LOC entries were renamed; ~397
> of their new symbols collided with HGNC entries already in the universe
> (multiple LOC paralogs collapsing to a single canonical symbol).

The R analysis (enrichR) uses enrichR's internal default background — enrichR
does not accept a custom background through its public API, so the
`tested_universe` set is used for ShinyGO only. Cross-tool agreement between
enrichR (default) and ShinyGO (tested-universe custom bg) is therefore an
informal sanity check, not a strict reproduction. No T2T-derived background
is used anywhere.

## Databases (independent passes)

Each database is queried independently — its own enrichR call, its own
re-FDR pool — so KEGG isn't penalised for being co-tested with thousands of
GO terms, etc.

| Pass | enrichR library | Term-size window | FDR cutoff |
|------|------------------|------------------|------------|
| `GO_BP`    | `GO_Biological_Process_2026` | 1 – 60   | 0.10 |
| `GO_MF`    | `GO_Molecular_Function_2026` | 1 – 60   | 0.10 |
| `KEGG`     | `KEGG_2026`                  | 1 – 60   | 0.10 |
| `DisGeNET` | `DisGeNET`                   | 1 – 1000 | 0.10 |

Per-DB size/FDR rationale:
- GO_BP, GO_MF, KEGG match the ShinyGO screenshot the user is reproducing
  (Pathway size Min = 1, Max = 60).
- DisGeNET keeps Max = 1000 so that **Cooley's anemia** (term_size = 70)
  remains visible — it's a biologically meaningful hit driven by HBA2 + GSTM1
  + LPA + PRH1.
- FDR is bumped from 0.05 → 0.10 ("as relaxed as you can") because the
  strongest GO BP/MF hits sit at FDR ≈ 0.064 under min=1/max=60 and would
  otherwise be lost. KEGG and DisGeNET top hits (Galactose metabolism,
  Laryngitis, Cooley's anemia) pass at FDR < 0.05 as well.

All three thresholds plus the database string are configurable per-pass in
the `passes` list near the top of `scripts/run_GO_enrichment.R`.

## Significance protocol (ShinyGO-style re-FDR)

`enrichR` returns raw `P.value` plus an `Adjusted.P.value` computed over the
entire DB. To match ShinyGO's behaviour (FDR limited to the size-filtered
pathway pool), the script:

1. Filters terms to `[min, max]` term_size.
2. Re-runs Benjamini-Hochberg FDR on the filtered subset using the raw
   `P.value` values.
3. Keeps rows below the per-pass FDR cutoff.

The re-FDR adjusted value is written as the primary `p_value` column; the
raw P.value is preserved in `p_value_raw` so the procedure is auditable.

## Output layout

```
GO_analysis/
├── README.md
├── gene_lists/
│   ├── cross_species_shared_60.txt
│   └── human_chimp_shared_35.txt
├── backgrounds/
│   └── ensembl_human_protein_coding.txt
├── cross_species_shared_60/
│   ├── enrichr_KEGG.csv       enrichr_GO_BP.csv       enrichr_GO_MF.csv       enrichr_DisGeNET.csv
│   ├── plot_KEGG.{pdf,png}    plot_GO_BP.{pdf,png}    plot_GO_MF.{pdf,png}    plot_DisGeNET.{pdf,png}
│   ├── shinygo_input.txt
│   └── shinygo_background_ensembl.txt
├── human_chimp_shared_35/      # same layout as above
├── plots/
│   ├── summary_KEGG.{pdf,png}
│   ├── summary_GO_BP.{pdf,png}
│   ├── summary_GO_MF.{pdf,png}
│   └── summary_DisGeNET.{pdf,png}
├── combined/
│   └── SupplementaryTable_enrichment.xlsx   # one sheet per gene list, all four passes
└── scripts/
    └── run_GO_enrichment.R
```

## Run

```sh
# install once
Rscript -e 'install.packages(c("enrichR", "openxlsx", "tidyverse"))'

cd svs_vep
Rscript GO_analysis/scripts/run_GO_enrichment.R
```

## Results overview (latest run; FDR < 0.10 on size-filtered subset, enrichR default bg)

| Gene list | GO_BP | GO_MF | KEGG | DisGeNET |
|---|---|---|---|---|
| `cross_species_shared_60` | 0 | 0 | **2** (Galactose metabolism, Glutathione metabolism) | **2** (Laryngitis, Cooley's anemia) |
| `human_chimp_shared_35`   | 112 | 39 | **9** (Galactose metabolism, **African trypanosomiasis**, **Malaria**, **Glutathione metabolism**, Starch/Sucrose, Sphingolipid, Glycosphingolipid biosynthesis, Citrate cycle, Carbohydrate digestion) | **1** (Laryngitis) |

Highlights:
- **ShinyGO KEGG hits reproduced.** The three pathways shown in the
  screenshot (African trypanosomiasis, Malaria, Glutathione metabolism) all
  appear in `enrichr_KEGG.csv` for `human_chimp_shared_35`. Their effective
  driver is HBA2 (trypanosomiasis, malaria) and the mucin/glycosylation
  genes (MGAM, NAT8, GSTM1) for the metabolism pathways.
- **Cooley's anemia preserved** in DisGeNET for the 60-gene list
  (HBA2 + GSTM1 + LPA + PRH1, FDR = 0.034).
- The **35-gene list** is much richer than the 60-gene set in GO_BP/MF
  because BH on a smaller query expands the proportion of size-1–6 single-
  gene "specific" terms surviving the size-filtered re-FDR. The most
  biologically interesting GO_BP signal in the 35-gene set is the *immune
  response-inhibiting receptor signalling* axis (LILRA6 + KIR2DL1).

## ShinyGO v0.85 — manual cross-check (optional)

Run for each gene list at <http://bioinformatics.sdstate.edu/go/>. Cite
Ge et al., 2020.

1. **Species** → `Human`.
2. **Gene list** → paste contents of `shinygo_input.txt`.
3. **Background** → "Custom" → paste `shinygo_background_tested.txt`
   (recommended: the 16,513 genes evaluated by the VEP pipeline). Alternative:
   `shinygo_background_ensembl.txt` for a wider whole-genome bg, or leave at
   ShinyGO default to match what enrichR uses internally.
4. **Parameters** (matching the R script):
   - FDR cutoff: **0.05** (or 0.10 if matching the relaxed R cutoff)
   - Min pathway size: **1**
   - Max pathway size: **60** (use **1000** if reproducing DisGeNET so anemia is in scope)
5. **Pathway databases to run** (switch dropdown and rerun for each):
   - KEGG
   - GO Biological Process
   - GO Molecular Function
   - DisGeNET (Disease)
6. **Download each results CSV** into `GO_analysis/<list>/shinygo_<DB>.csv`.

## Reproducibility

`run_GO_enrichment.R` is the single source of truth. Re-running it
regenerates identical gene-list files and enrichR CSVs (enrichR queries are
deterministic for a fixed library snapshot). Database snapshots evolve over
time, so the exact term names may drift if Enrichr updates `KEGG_2026` or
`GO_*_2026`. Pin the library names in `passes` if you need a frozen run.
