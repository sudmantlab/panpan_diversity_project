import csv

# =============================================================================
# Mapping: description (from CSV) -> genes_renamed (human gene symbol)
# Sources:
#   - NCBI Gene database lookups (HIGH confidence)
#   - User-provided reference list (cross-checked against NCBI)
#   - Inferred from description (clearly embedded gene names)
# Notes on ambiguous/low-confidence cases are inline.
# =============================================================================

# LOC gene descriptions -> human gene symbol (NCBI-verified)
LOC_MAPPING = {
    # SSU72 family (9 bonobo LOC genes, all highly duplicated from SSU72)
    "RNA polymerase II subunit A C-terminal domain phosphatase SSU72 like protein 2":       "SSU72L2",
    "RNA polymerase II subunit A C-terminal domain phosphatase SSU72 like protein 2-like":  "SSU72L2",

    # Tubulin
    "tubulin beta 8B-like":         "TUBB8B",
    "tubulin beta-8 chain-like":    "TUBB8",
    "tubulin beta chain-like":      "TUBB-like",

    # PGM5
    "phosphoglucomutase-like protein 5": "PGM5",

    # PRAME family
    "PRAME family member 1-like":   "PRAMEF1",
    "PRAME family member 1":        "PRAMEF1",
    "PRAME family member 8":        "PRAMEF8",
    "PRAME family member 8-like":   "PRAMEF8",
    "PRAME family member 11-like":  "PRAMEF11",
    "PRAME family member 15-like":  "PRAMEF15",
    "PRAME family member 15":       "PRAMEF15",

    # RFPL4A family
    "ret finger protein-like 4A":   "RFPL4A",

    # SDHA
    "succinate dehydrogenase [ubiquinone] flavoprotein subunit%2C mitochondrial-like": "SDHA",
    "succinate dehydrogenase [ubiquinone] flavoprotein subunit%2C mitochondrial":      "SDHA",

    # Olfactory receptors
    "olfactory receptor 6A2":   "OR6A2",
    "olfactory receptor 2T12":  "OR2T12",

    # MHC / HLA (gorilla Gogo-* alleles -> human HLA locus)
    "class I histocompatibility antigen%2C Gogo-C*0202 alpha chain": "HLA-C",
    "class I histocompatibility antigen%2C Gogo-B*0101 alpha chain": "HLA-B",
    "HLA class II histocompatibility antigen%2C DP beta 1 chain":    "HLA-DPB1",

    # Golgins
    "golgin subfamily A member 8B":         "GOLGA8B",
    "putative golgin subfamily A member 8I": "GOLGA8IP",
    "putative golgin subfamily A member 8G": "GOLGA8G",

    # WASHC2 (two human paralogs; WASHC2A used here per user's reference)
    "WAS protein family homolog 2": "WASHC2A",

    # Phospholipase A2 — description too vague to assign confidently; keep LOC name
    # "phospholipase A2-like" -> remains as-is (LOC117976337)

    # ADAMTS
    "A disintegrin and metalloproteinase with thrombospondin motifs 7-like": "ADAMTS7",

    # SRRM1
    "serine/arginine repetitive matrix protein 1-like": "SRRM1",

    # EIF5B (translation initiation factor IF-2)
    "translation initiation factor IF-2": "EIF5B",

    # ZNF genes
    "zinc finger protein 69":  "ZNF69",
    "zinc finger protein 440": "ZNF440",

    # LILR family
    "leukocyte immunoglobulin-like receptor subfamily B member 1": "LILRB1",
    "leukocyte immunoglobulin-like receptor subfamily A member 2": "LILRA2",
    "leukocyte immunoglobulin-like receptor subfamily B member 3": "LILRB3",

    # EEF2KMT
    "protein-lysine N-methyltransferase EEF2KMT-like": "EEF2KMT",

    # ANTXRL
    "anthrax toxin receptor-like": "ANTXRL",

    # AGAP14 (pseudogene in humans but best available symbol)
    "arf-GAP with GTPase%2C ANK repeat and PH domain-containing protein 14-like": "AGAP14P",

    # CLASRP
    "CLK4-associating serine/arginine rich protein-like": "CLASRP",

    # PARP4
    "protein mono-ADP-ribosyltransferase PARP4-like": "PARP4",

    # CASP/CUX1 — ambiguous alias; keep as LOC (omitted from mapping)

    # RASA4
    "ras GTPase-activating protein 4-like": "RASA4",
    "ras GTPase-activating protein 4":      "RASA4",

    # CNTN4
    "contactin-4": "CNTN4",

    # MEIG1
    "meiosis expressed gene 1 protein homolog": "MEIG1",

    # SPTBN5
    "spectrin beta chain%2C non-erythrocytic 5-like": "SPTBN5",

    # FAM157A
    "putative protein FAM157A": "FAM157A",

    # NBPF15
    "neuroblastoma breakpoint family member 15-like": "NBPF15",

    # ZNG1A
    "zinc-regulated GTPase metalloprotein activator 1A-like": "ZNG1A",

    # NPIPA7
    "nuclear pore complex-interacting protein family member A7": "NPIPA7",

    # TAF11L2
    "TATA-box binding protein associated factor 11 like protein 2": "TAF11L2",

    # NPIPB8 / B9
    "nuclear pore complex-interacting protein family member B8-like": "NPIPB8",
    "nuclear pore complex-interacting protein family member B9-like": "NPIPB9",

    # CELA3B
    "chymotrypsin-like elastase family member 3B": "CELA3B",

    # HNRNPCL2/3
    "heterogeneous nuclear ribonucleoprotein C-like 2": "HNRNPCL2",
    "heterogeneous nuclear ribonucleoprotein C-like 3": "HNRNPCL3",

    # SPDYE family
    "putative speedy protein E7":   "SPDYE7P",
    "speedy protein E12":           "SPDYE12",
    "speedy protein E1-like":       "SPDYE1",
    "speedy protein E1":            "SPDYE1",
    "speedy protein E18":           "SPDYE18",

    # TEKT1
    "tektin-1-like": "TEKT1",

    # ZDHHC11
    "palmitoyltransferase ZDHHC11":        "ZDHHC11",
    "palmitoyltransferase ZDHHC11-like":   "ZDHHC11",

    # SRP68
    "signal recognition particle subunit SRP68-like": "SRP68",

    # FAM153A
    "protein FAM153A": "FAM153A",

    # PRH1
    "proline-rich protein HaeIII subfamily 1": "PRH1",

    # PKD1L2
    "polycystin-1-like protein 2": "PKD1L2",

    # BICRA (BRD4-interacting chromatin-remodeling complex-associated protein)
    "BRD4-interacting chromatin-remodeling complex-associated protein": "BICRA",

    # CDRT15
    "CMT1A duplicated region transcript 15 protein-like protein": "CDRT15",

    # KPRP
    "keratinocyte proline-rich protein-like": "KPRP",

    # HAGH
    "hydroxyacylglutathione hydrolase%2C mitochondrial": "HAGH",

    # FOXD4L1
    "forkhead box protein D4-like 1": "FOXD4L1",

    # ABCA17P (human ortholog is a pseudogene)
    "ATP-binding cassette sub-family A member 17-like": "ABCA17P",

    # DCLRE1C (Artemis)
    "protein artemis-like": "DCLRE1C",

    # APBA2
    "amyloid-beta A4 precursor protein-binding family A member 2-like": "APBA2",

    # RBMXL2
    "RNA-binding motif protein%2C X-linked-like-2": "RBMXL2",

    # RSPH10B
    "radial spoke head 10 homolog B": "RSPH10B",

    # FAHD1
    "acylpyruvase FAHD1%2C mitochondrial": "FAHD1",

    # PMS2L (putative postmeiotic segregation increased 2-like; no canonical HGNC symbol)
    "putative postmeiotic segregation increased 2-like protein 1": "PMS2L1",

    # PDPR
    "pyruvate dehydrogenase phosphatase regulatory subunit%2C mitochondrial-like": "PDPR",

    # EYS
    "eyes shut homolog": "EYS",

    # AQP12B
    "aquaporin-12B": "AQP12B",

    # BCL9L
    "B-cell CLL/lymphoma 9-like protein": "BCL9L",

    # MBD3L3
    "putative methyl-CpG-binding domain protein 3-like 3": "MBD3L3",

    # UGT2A1
    "UDP-glucuronosyltransferase 2A1": "UGT2A1",

    # FAM231C
    "FAM231A/C-like protein": "FAM231C",

    # FAM90A27P (keep from description; no clear single human ortholog)
    "protein FAM90A27P": "FAM90A27P",

    # SHISA5
    "protein shisa-5": "SHISA5",

    # TNXB
    "tenascin-X-like": "TNXB",

    # Phospholipase A2 — note: ambiguous; left as LOC
    # "phospholipase A2-like": <-- intentionally omitted; insufficient info

    # Papa-A (bonobo MHC -> HLA-A; Patr-A is chimp ortholog, Papa-A is bonobo)
    "patr class I histocompatibility antigen%2C A-2 alpha chain-like": "HLA-A",

    # ----------------------
    # Additional entries from user-provided reference list (non-LOC genes,
    # or supplementary descriptions not covered above)
    # ----------------------
    "alpha/beta hydrolase domain-containing protein 17A":                           "ABHD17A",
    "acyl-coenzyme A synthetase ACSM1%2C mitochondrial-like":                      "ACSM1",
    "A disintegrin and metalloproteinase with thrombospondin motifs 7":             "ADAMTS7",
    "AFG3-like protein 1":                                                          "AFG3L1",
    "chitobiosyldiphosphodolichol beta-mannosyltransferase-like":                   "ALG1",
    "ankyrin repeat domain-containing protein 36B":                                 "ANKRD36B",
    "breast carcinoma-amplified sequence 4-like":                                   "BCAS1",
    "BCL-6 corepressor-like protein 1":                                             "BCORL1",
    "putative protein C3P1":                                                        "C3P1",
    "complement C4-A":                                                              "C4A",
    "protein CASC2%2C isoforms 1/2":                                               "CASC2",
    "protein CASP":                                                                 "CASP",
    "putative coiled-coil domain-containing protein 144B":                          "CCDC144B",
    "dual specificity protein phosphatase CDC14C-like":                             "CDC14C",
    "complement factor H-related protein 4":                                        "CFHR4",
    "collagen alpha-4(VI) chain-like":                                              "COL6A4",
    "chymotrypsinogen B":                                                           "CTRB1",
    "ATP-dependent DNA helicase DDX11-like":                                        "DDX11",
    "eukaryotic initiation factor 4A-III":                                          "EIF4A-III",
    "eukaryotic initiation factor 4A-III-like":                                     "EIF4A-III-like",
    "endogenous retrovirus group K member 16 Rec protein":                          "ERVK-16",
    "embryonic stem cell-related gene protein-like":                                "ESRG",
    "protein eyes shut homolog":                                                    "EYS",
    "putative protein FAM157A":                                                     "FAM157A",
    "protein FAM90A27P-like":                                                       "FAM90A27P",
    "cytosolic beta-glucosidase":                                                   "GBA3",
    "putative golgin subfamily A member 8G":                                        "GOLGA8G",
    "putative golgin subfamily A member 8I":                                        "GOLGA8IP",
    "protein GVQW1-like":                                                           "GVQW1",
    "glycophorin A-like protein":                                                   "GYPx_A-like",
    "glycophorin B (MNS blood group)":                                              "GYPx_B-like",
    "glycophorin-B":                                                                "GYPx_B-like",
    "hemoglobin subunit alpha-3-like":                                              "HBA3-like",
    "protein PBMUCL2":                                                              "HCG22",
    "putative HERC2-like protein 3":                                                "HERC2P3",
    "HLA class I histocompatibility antigen%2C A alpha chain":                     "HLA-A",
    "patr class I histocompatibility antigen%2C A-126 alpha chain":               "HLA-A",
    "major histocompatibility complex%2C class I%2C B":                           "HLA-B",
    "major histocompatibility complex%2C class I%2C C":                           "HLA-C",
    "major histocompatibility complex%2C class II%2C DP alpha 1":                 "HLA-DPA1",
    "HLA class II histocompatibility antigen%2C DP beta 1 chain-like":            "HLA-DPB1",
    "major histocompatibility complex%2C class II%2C DP beta 1":                  "HLA-DPB1",
    "major histocompatibility complex%2C class II%2C DQ alpha 1":                 "HLA-DQA1",
    "major histocompatibility complex%2C class II%2C DQ alpha 2":                 "HLA-DQA2",
    "HLA class II histocompatibility antigen%2C DQ beta 1 chain":                 "HLA-DQB1",
    "major histocompatibility complex%2C class II%2C DQ beta 1":                  "HLA-DQB1",
    "major histocompatibility complex%2C class II%2C DR alpha":                   "HLA-DRA",
    "HLA class II histocompatibility antigen%2C DR beta 3 chain":                 "HLA-DRB3",
    "major histocompatibility complex%2C class I%2C F":                           "HLA-F",
    "heterogeneous nuclear ribonucleoprotein C-like 2":                             "HNRNPCL2",
    "putative HTLV-1-related endogenous sequence":                                  "HRES-1",
    "heat shock protein HSP 90-alpha":                                              "HSP90AA1",
    "alpha-L-iduronidase":                                                          "IDUA",
    "interferon-induced transmembrane protein 3-like":                              "IFITM3",
    "immunoglobulin heavy constant gamma 1-like":                                   "IGHG1",
    "immunoglobulin gamma-1 heavy chain-like":                                      "IGHG1",
    "immunoglobulin heavy variable 3-74-like":                                      "IGHV4-3-74",
    "immunoglobulin heavy variable 4-30-4-like":                                    "IGHV4-30-4",
    "immunoglobulin heavy variable 4-34-like":                                      "IGHV4-34",
    "immunoglobulin heavy variable 4-38-2-like":                                    "IGHV4-38-2",
    "immunoglobulin heavy variable 5-51-like":                                      "IGHV5-51",
    "immunoglobulin superfamily member 1":                                          "IGSF1",
    "interleukin-9 receptor-like":                                                  "IL-9R",
    "bifunctional peptidase and (3S)-lysyl hydroxylase JMJD7-like":               "JMJD7",
    "bifunctional peptidase and (3S)-lysyl hydroxylase JMJD7":                    "JMJD7",
    "kinesin-like protein KIF28P":                                                  "KIF28P",
    "leukocyte immunoglobulin-like receptor subfamily B member 3-like":             "LILRB3",
    "MAM and LDL-receptor class A domain-containing protein 1-like":               "MALRD1",
    "MAM and LDL-receptor class A domain-containing protein 2-like":               "MALRD2",
    "meteorin-like protein":                                                        "METRNL",
    "MHC class I polypeptide-related sequence B-like":                             "MICB",
    "protein MOST-1":                                                               "MOST-1",
    "mucin-17":                                                                     "MUC17",
    "mucin-20":                                                                     "MUC20",
    "DNA mismatch repair protein MutL-like":                                        "MUTL",
    "DNA mismatch repair protein MutL":                                             "MUTL",
    "myosin light chain kinase 2%2C skeletal/cardiac muscle-like":                "MYLK2",
    "neuroblastoma breakpoint family member 15-like":                               "NBPF15",
    "serine/threonine-protein kinase Nek4-like":                                    "NEK4",
    "transcription factor NF-E4-like":                                              "NFE4",
    "nuclear pore complex-interacting protein family member B15-like":              "NPIPB15",
    "nuclear pore complex-interacting protein family member B15":                   "NPIPB15",
    "nuclear pore complex-interacting protein family member B8-like":               "NPIPB8",
    "olfactomedin-4-like":                                                          "OLFM4",
    "putative olfactory receptor 10J6":                                             "OR10J6P",
    "olfactory receptor 1L4-like":                                                  "OR1L4",
    "olfactory receptor 2H1-like":                                                  "OR2H1",
    "olfactory receptor 4F3/4F16/4F29":                                             "OR4F3",
    "olfactory receptor 51F2":                                                      "OR51F2",
    "olfactory receptor 51S1":                                                      "OR51S1",
    "olfactory receptor 52A5":                                                      "OR52A5",
    "olfactory receptor 5H6-like":                                                  "OR5H6",
    "olfactory receptor 8U9-like":                                                  "OR8U9",
    "protein mono-ADP-ribosyltransferase PARP4-like":                              "PARP4",
    "rod cGMP-specific 3'%2C5'-cyclic phosphodiesterase subunit beta":            "PDE6B",
    "cytosolic phospholipase A2 beta":                                              "PLA2G4A",
    "cytosolic phospholipase A2 beta-like":                                         "PLA2G4B",
    "putative postmeiotic segregation increased 2-like protein 2":                  "PMS2",
    "mismatch repair endonuclease PMS2":                                            "PMS2",
    "PRAME family member 1":                                                        "PRAMEF1",
    "PRAME family member 11-like":                                                  "PRAMEF11",
    "PRAME family member 15":                                                       "PRAMEF15",
    "PRAME family member 2":                                                        "PRAMEF2",
    "PRAME family member 8-like":                                                   "PRAMEF8",
    "salivary acidic proline-rich phosphoprotein 1/2":                             "PRH1",
    "proline-rich protein 23D1-like":                                               "PRR23D1",
    "pregnancy-specific beta-1-glycoprotein 11":                                    "PSG11",
    "prothymosin alpha-like":                                                       "PTMA",
    "RANBP2-like and GRIP domain-containing protein 5/6":                          "RANBP2",
    "E3 SUMO-protein ligase RanBP2-like":                                           "RANBP2",
    "ras GTPase-activating protein 4":                                              "RASA4",
    "ret finger protein-like 4A":                                                   "RFPL4A",
    "RH-like protein IC":                                                           "RH-like-IC",
    "RH-like protein IIR":                                                          "RH-like-IIR",
    "U6 spliceosomal RNA":                                                          "RNU6-1",
    "large ribosomal subunit protein eL19":                                         "RPL19",
    "large ribosomal subunit protein uL2":                                          "RPL8",
    "small ribosomal subunit protein uS8-like":                                     "RPS15A",
    "putative RRN3-like protein RRN3P2":                                            "RRN3P2",
    "radial spoke head 10 homolog B":                                               "RSPH10B",
    "ras suppressor protein 1-like":                                                "RSU1",
    "protein salvador homolog 1-like":                                              "SAV1",
    "SUMO-interacting motif-containing protein 1":                                  "SIMC1",
    "putative solute carrier organic anion transporter family member 1B7":          "SLCO1B7",
    "small nucleolar RNA U13":                                                      "SNORD13",
    "spectrin beta chain%2C non-erythrocytic 5-like":                             "SPTBN5",
    "TAFA chemokine like family member 1":                                          "TAFA1",
    "TATA-box-binding protein-associated factor 11-like protein 5":                "TAFA5",
    "taste receptor type 2 member 30-like":                                         "TAS2R14",
    "T-box transcription factor TBX1-like":                                         "TBX1",
    "putative T-complex protein 10A homolog":                                       "TCP10L",
    "T cell receptor alpha variable 14/delta variable 4-like":                     "TRAV14DV4",
    "T cell receptor alpha variable 8-6-like":                                      "TRAV8-6",
    "T cell receptor beta variable 7-6-like":                                       "TRBV7-6",
    "E3 ubiquitin-protein ligase TRIM52":                                           "TRIM52",
    "tubulin beta-8 chain-like":                                                    "TUBB8",
    "S-adenosyl-L-methionine-dependent tRNA 4-demethylwyosine synthase TYW1":     "TYW1",
    "UDP-glucuronosyltransferase 2B15":                                             "UGT2B15",
    "putative UPF0607 protein ENSP00000383144":                                     "UPF0607",
    "WASH complex subunit 2A-like":                                                 "WASHC2A",
    "putative protein ZNF321":                                                      "ZNF321",
    "zinc finger protein 665-like":                                                 "ZNF665",
    "zinc finger protein 701":                                                      "ZNF701",
    "zinc finger protein 813":                                                      "ZNF813",
    "zinc finger protein 83":                                                       "ZNF83",
    "putative COBW domain-containing protein 7":                                    "ZNG1F",
}

# =============================================================================
# Apply mapping to filtered CSV
# =============================================================================
rows_out = []
renamed_count = 0
unresolved = []

with open('/Users/joanocha/Desktop/singer_tmp/bonobo_ranked_avg_pairwise_myr_filtered.csv') as f:
    reader = csv.reader(f)
    header = next(reader)

    # Strip any pre-existing genes_renamed columns to avoid duplicates
    while 'genes_renamed' in header:
        idx = header.index('genes_renamed')
        header.pop(idx)

    # Insert genes_renamed right after genes (index 5)
    genes_idx = header.index('genes')
    header_out = header[:genes_idx + 1] + ['genes_renamed'] + header[genes_idx + 1:]

    for row in reader:
        # Strip stale genes_renamed columns if present
        row = row[:len(header)]

        desc = row[7].strip()
        gene = row[5].strip()
        renamed = LOC_MAPPING.get(desc, gene)
        if renamed != gene:
            renamed_count += 1
        elif gene.startswith('LOC') or gene == 'Papa-A':
            unresolved.append((gene, desc))

        row_out = row[:genes_idx + 1] + [renamed] + row[genes_idx + 1:]
        rows_out.append(row_out)

with open('/Users/joanocha/Desktop/singer_tmp/bonobo_ranked_avg_pairwise_myr_filtered.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header_out)
    writer.writerows(rows_out)

print(f'Done. {renamed_count} genes renamed out of {len(rows_out)} rows.')

if unresolved:
    print(f'\nLOC genes with no mapping found ({len(unresolved)}):')
    for gene, desc in unresolved:
        print(f'  {gene:30s}  {desc[:80]}')
