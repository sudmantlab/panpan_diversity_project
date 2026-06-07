# Preliminary dataset

A digest of the sample material underlying the PANPAN project — the newly sequenced
samples, the public NCBI re-sequencing panel used for sample selection, and the Coriell
family/pedigree information. Each section links to the detailed file.

## Contents

| File | What it holds |
|---|---|
| [`Sample_information.md`](Sample_information.md) | Samples sequenced for this project (short-read, PacBio Sequel II / Revio, HiC, ONT) |
| `Sample_infoNCBI.py` *(kept locally, not tracked)* | Public NCBI/SRA re-sequencing accessions used for sample selection |
| [`Coriell_apes_families.md`](Coriell_apes_families.md) | Pedigrees and parentage of the Coriell bonobo/chimpanzee samples |
| [`pca/`](pca) | PCA of the low-coverage chimpanzee panel (`pca.R`, covariance + eigen outputs) |

---

## 1. Sample information (sequenced for PANPAN)

See **[`Sample_information.md`](Sample_information.md)**. Covers the samples generated for
this study:

- **Illumina low-coverage WGS (~4×)** — 21 *Pan troglodytes* individuals, used to pick a
  diverse, geographically representative set for long-read sequencing and downstream COSIGT
  genotyping.
- **PacBio Sequel II (Phase I)** — 17 individuals (bonobos + chimpanzee subspecies, incl. a
  western×central hybrid), some with additional HiC / ONT.
- **PacBio Revio (Part 2)** — 13 individuals, all with HiC (one, `AG18352_2`, was dropped
  during cell culture).

Reference-genome cross-links:
- **Bonobo `PR00251`** = T2T reference **mPanPan1**.
- **Chimpanzee `AG18354_5`** = T2T reference **mPanTro3**.

## 2. NCBI public re-sequencing panel

See `Sample_infoNCBI.py` *(kept locally, not tracked in the repo)* — a dictionary mapping
BioSamples to their SRA/ENA run accessions for publicly available great-ape WGS
(Prado-Martinez / de Manuel *et al.* great-ape diversity data), used to contextualize and
select samples.

- **60 samples**, **416 SRA/ENA runs** total.
- By taxon: *P. t. schweinfurthii* (eastern) 19 · *P. t. troglodytes* (central) 18 ·
  *P. t. verus* (western) 11 + 1 verus×troglodytes hybrid · *P. t. ellioti*
  (Nigeria–Cameroon) 10 · *P. paniscus* (bonobo) 1.
- Keys encode `BioSample_subspecies_name_country`; values are the per-sample run lists.

## 3. Coriell ape families

See **[`Coriell_apes_families.md`](Coriell_apes_families.md)** — five Coriell families
(2 bonobo, 3 chimpanzee) reconstructed from hand-drawn pedigrees and cross-referenced to the
Coriell distribution list (names, studbook/KB# IDs, sire/dam, and `PR` catalog IDs):

- **F1 (bonobos):** 3-generation pedigree around sire Bosondjo (trio + 3 paternal half-sibs).
- **F2 (bonobos):** Loretta → Erin (mother–daughter; Erin = mPanPan1).
- **F3 (chimps):** Hanky & Goober, paternal half-sibs.
- **F4 (chimps):** Blackie → Kioja.
- **F5 (chimps):** Joshua & Rachel, maternal half-sibs.

Parents drawn without a `PR` ID are unsampled (local/ISIS/studbook IDs only).
