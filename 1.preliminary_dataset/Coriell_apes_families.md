# Coriell ape families — pedigree information

Five families of Coriell bonobo (*Pan paniscus*) and chimpanzee (*Pan troglodytes*)
samples, originally sketched as hand-drawn pedigrees. The pedigree drawings were
cross-referenced against the Coriell **distribution list** and **families** tables to
correct names, studbook IDs, and parentage.

**Sources**
- Hand-drawn pedigrees (`Coriell_apes_families.pdf`) — family groupings and pedigree topology.
- Coriell sample tables (`Matt_old_tables/Pan_samples_list_relationships_notated.xlsx`,
  sheets `Distribution_list` and `Families`) — names, studbook/KB# IDs, and sire/dam IDs.

*(The PDF and the xlsx source tables are kept locally and are git-ignored; this file is the
shareable transcription.)*

**Pedigree symbol key**

| Symbol | Meaning |
|---|---|
| ○ circle | female |
| □ square | male |
| ◇ diamond | offspring (sex not indicated in drawing) |

Each individual carries a **studbook / KB# ID** and a **Coriell catalog ID** (`PR#####`).
Sire/dam values shown as bare numbers (e.g. `360`, `0835`, `169118`) are local/ISIS/studbook
IDs of **unsampled** parents — they have no Coriell `PR` sample.

---

## Family 1 — yellow (bonobos, *Pan paniscus*)

A three-generation bonobo pedigree centered on the sire **Bosondjo**: one trio
(Bosondjo × Matata → Panbanisha) and three paternal half-sibs (Panbanisha, Akili, Lisala).

| Name | Studbook / KB# | Coriell ID | Sex | Role / parents |
|---|---|---|---|---|
| Bosondjo ("Boso") | 12730 | PR00111 | ♂ | sire; father of Akili, Lisala, Panbanisha |
| Matata | 7189 | PR00367 | ♀ | dam of Panbanisha (× Bosondjo) |
| Panbanisha | 8764 | PR00366 | ♀ | offspring of Bosondjo × Matata |
| Akili | 5274 | PR00236 | ♂ | offspring of Bosondjo (dam P1); father of Jumanji |
| Lisala | 8798 | PR00748 | ♀ | offspring of Bosondjo (dam 169118, unsampled) |
| Jumanji | 10227 | PR00802 | — | offspring of Akili (dam 587376, unsampled); 3rd generation |

## Family 2 — green (bonobos, *Pan paniscus*)

Mother–daughter pair.

| Name | Studbook / KB# | Coriell ID | Sex | Role / parents |
|---|---|---|---|---|
| Loretta | 2535 | PR00235 | ♀ | mother of Erin |
| Erin | 8711 | PR00251 | ♀ | daughter of Loretta (sire ISIS 180343, unsampled) |

Note: **PR00251 (Erin) is the same individual as the T2T reference mPanPan1.**

## Family 3 — orange (chimpanzees, *Pan troglodytes*)

Two paternal half-sibs sharing sire `360`.

| Name | Studbook / KB# | Coriell ID | Sex | Role / parents |
|---|---|---|---|---|
| Hanky | 12796 | PR00826 | — | offspring; sire 360, dam 364 |
| Goober | 13273 | PR00400 | — | offspring; sire 360, dam 387 |

Both parents of each individual (sires 360; dams 364/387) are unsampled.

## Family 4 — blue (chimpanzees, *Pan troglodytes*)

| Name | Studbook / KB# | Coriell ID | Sex | Role / parents |
|---|---|---|---|---|
| Blackie | 13470 | PR00549 | ♀ | mother (ISIS 298) of Kioja |
| Kioja | 13471 | PR00548 | — | offspring; sire 297 (unsampled), dam 298 (Blackie) |

## Family 5 — red (chimpanzees, *Pan troglodytes*)

Two maternal half-sibs sharing dam `0835`.

| Name | Studbook / KB# | Coriell ID | Sex | Role / parents |
|---|---|---|---|---|
| Joshua | 10512 | PR01171 | — | offspring; sire 1021, dam 0835 |
| Rachel | 10790 | PR00818 | — | offspring; sire 0834, dam 0835 |

The three parents (sires `1021`, `0834`; shared dam `0835`) are local/ISIS IDs and are unsampled.

---

## Membership in the PANPAN sequenced set

Of the individuals above, the following are part of the project's sequenced sample set
(present in `Sample_information.md`):

**PR00251, PR00366, PR00400, PR00548, PR00818, PR00826, PR01171.**

The remainder (PR00111, PR00235, PR00236, PR00367, PR00549, PR00748, PR00802) appear in the
Coriell catalog but not in the project's sequenced samples.
