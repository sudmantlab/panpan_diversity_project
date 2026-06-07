# Data Processing

NOTE:
BONOBO PR00251 is the same individual as T2T mPanPan1
CHIMPANZEE AG18354_5 is the same individual as T2T mPanTro3


### STATUS

#### Current Status

| Name     | Species                                        | Data                                        |  STATUS |
| -------- | ---------------------------------------------- | ------------------------------------------- | ---------------------------------------------- |
| PR01227  | Pan troglodytes (verus)                        | PacBIO HiFi, HiC, ONT-UL                    | assembled with hifiasm, hifiasm-hic, verkko 2.0, updated with verkko 2.2.1                         |
| PR01225  | Pan troglodytes (verus)                        | PacBIO HiFi                                 | assembled with hifiasm                                              |
| PR01223  | Pan troglodytes (verus)                        | PacBIO HiFi                                 | assembled with hifiasm                                              |
| PR01228  | Pan troglodytes (verus)                        | PacBIO HiFi                                 | assembled with hifiasm                                             |
| PR00512  | Pan troglodytes (verus)                        | PacBIO HiFi                                 | assembled with hifiasm                                            |
| PR01100  | Pan troglodytes (verus)                        | PacBIO HiFi                                 | assembled with hifiasm                                      |
| PR00834  | Pan troglodytes (verus x troglodytes hybrid)                  | PacBIO HiFi, HiC, ONT-UL                    | assemled with hifiasm, hifiasm-HiC                     |
| PR00838  | Pan troglodytes (verus x troglodytes hybrid)                  | PacBIO HiFi                                 | assemled with hifiasm                                           |
| PR00115  | Pan troglodytes (troglodytes)                  | PacBIO HiFi                                 | assemled with hifiasm                                             |
| PR00251  | Pan paniscus                                   | PacBIO HiFi, HiC \*DUPLICATE w T2T effort\* | assemled with hifiasm, hifiasm-HiC, mPanPan1 T2T assembly            |
| PR00366  | Pan paniscus                                   | PacBIO HiFi, HiC, ONT-UL                    | assembled with hifiasm, hifiasm-hic, verkko 2.1, updated with verkko 2.2.1                    |
| PR00445  | Pan paniscus                                   | PacBIO HiFi                                 | assemled with hifiasm                                              |
| PR00249  | Pan paniscus                                   | PacBIO HiFi                                 | assemled with hifiasm                                              |
| PR01008  | Pan troglodytes (schweinfurthii)               | PacBIO HiFi, HiC, ONT-UL                    | assembled with hifiasm, hifiasm-hic, verkko 2.0, updated with verkko 2.2.1                        |
| PR01009  | Pan troglodytes (schweinfurthii)               | PacBIO HiFi                                 | assemled with hifiasm                                               |
| PR01010  | Pan troglodytes (schweinfurthii)               | PacBIO HiFi, HiC, ONT-UL                    | assembled with hifiasm, hifiasm-hic, verkko 2.0, updated with verkko 2.2.1                                |
| PR00496A | Pan troglodytes (troglodytes) | PacBIO HiFi                                 |  assembled with hifiasm                                             |
| AG05253_1 | Pan paniscus     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                               | 
| AG18352_2 | ??     | PacBIO HiFi + HiC  | FAILED/REMOVED                                               | 
| AG18357_3 | Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                              | 
| AG18358_4 | Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                              | 
| AG18354_5 | Pan troglodytes (verus)     | PacBIO HiFi, HiC \*DUPLICATE w T2T effort\*  | assemled with hifiasm-HiC, mPanTro3 T2T assembly                                              | 
| AG18356_6 | Pan troglodytes (verus x troglodytes hybrid)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                                | 
| AG18359_7 | Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                            | 
| AG18355_8 | Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                               | 
| AG18353_9 | Pan troglodytes (verus)    | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                               | 
| AG18361_10| Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                              |  
| AG18360_11| Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                               |
| AG16618_12| Pan troglodytes (verus x troglodytes hybrid)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                               |  
| AG06939_13| Pan troglodytes (verus)     | PacBIO HiFi + HiC  | assemled with hifiasm-HiC                                               | 


| No. | Description:{sample}                          | Hifiasm Hi-Fi assemblies                                   | Hifiasm Hi-C assemblies                                    | Verkko assemblies |
|-----|-----------------------------------------------|------------------------------------------------------------|------------------------------------------------------------|-------------------|
| 1   | Pan troglodytes verus:PR01228                | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 2   | Pan troglodytes verus:PR01227                | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         | PANPAN/output2/{sample}/output/asm/ |
| 3   | Pan troglodytes verus:PR01225                | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 4   | Pan troglodytes verus:PR01223                | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 5   | Pan troglodytes verus:PR01100                | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 6   | Pan troglodytes schweinfurthii:PR01010       | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         | PANPAN/output2/{sample}/output/asm/ |
| 7   | Pan troglodytes schweinfurthii:PR01009       | PANPAN/output/{sample}/hifiasm-fasta/                     |                                                            |                   |
| 8   | Pan troglodytes schweinfurthii:PR01008       | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         | PANPAN/output2/{sample}/output/asm/ |
| 9   | Pan troglodytes troglodytes:PR00838          | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 10  | Pan troglodytes troglodytes:PR00834          | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 11  | Pan troglodytes verus:PR00512                | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 12  | Pan troglodytes troglodytes?:PR00496A        | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 13  | Pan paniscus:PR00445                         | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 14  | Pan paniscus:PR00366                         | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 15  | Pan paniscus:PR00251                         | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 16  | Pan paniscus:PR00249                         | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 17  | Pan troglodytes troglodytes:PR00115          | PANPAN/output2/{sample}/output/hifiasm-fasta/             |                                                            |                   |
| 18  | Pan paniscus:AG05253_1                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 19  | Sample:AG18352_2                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 20  | Sample:AG18357_3                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 21  | Sample:AG18358_4                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 22  | Sample:AG18354_5                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 23  | Sample:AG18356_6                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 24  | Sample:AG18359_7                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 25  | Sample:AG18355_8                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 26  | Sample:AG18353_9                             | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 27  | Sample:AG18361_10                            | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 28  | Sample:AG18360_11                            | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 29  | Sample:AG16618_12                            | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |
| 30  | Sample:AG06939_13                            | PANPAN/output2/{sample}/output/hifiasm-fasta/             | PANPAN/output2/{sample}/output/hifiasm-fasta-HiC/         |                   |




### COMPLETED (PRIMARY AND PHASED) ASSEMBLIES (example with a sample PR01227)
to run assembly pipeline with hifiasm only:
```
snakemake -prj40  --rerun-incomplete --keep-going --use-conda --conda-frontend mamba output/genomescope2/PR01227/joana_settings_shortcut/chimpPR01227/linear_plot.png output/merqury/PR01227/joana_settings_shortcut/no_opts/PR01227.{haplotig,consensus} output/gt-seqstat/PR01227/joana_settings_shortcut/no_opts/PR01227.p_ctg.stats
```

to run assembly pipeline with hifiasm + hic:
```
snakemake -prj40  --rerun-incomplete --keep-going --use-conda --conda-frontend mamba output/hifiasm-fasta-HiC/PR01227/joana_settings_shortcut/no_opts/PR01227.hap1.p_ctg.hic.fa  output/hifiasm-fasta-HiC/PR01227/joana_settings_shortcut/no_opts/PR01227.hap2.p_ctg.hic.fa output/hifiasm-fasta-HiC/PR01227/joana_settings_shortcut/no_opts/PR01227.p_ctg.hic.fa 
```

to run assembly pipeline with hifiasm + hic + ont (verkko):
```
verkko -d asm  --hifi HiFi-adapterFiltered/PR01227/joana_settings_shortcut/mtSinai/l0000008*/*.ccs.filt.fastq.gz --nano /global/scratch/users/joana_rocha/PANPAN/PANPAN_ont/PR01227/*/*_pass.fastq.gz --hic1 trimmed-hic/PR01227/PR01227_CKDL230002181-1A_HMMF5DSX5_L2-trimmed_R1.fastq.gz --hic2 trimmed-hic/PR01227/PR01227_CKDL230002181-1A_HMMF5DSX5_L2-trimmed_R2.fastq.gz --snakeopts "--cores all" --snakeopts "--rerun-incomplete"
```


### COMPLETED MT GENOMES

to run mitoHIFI:
```
mkdir .singularity
ln -s /global/scratch/users/joana_rocha/.singularity ~/.singularity
```

```
#snakemake -prj40 --keep-going --use-singularity --singularity-args "-B /global/scratch -B $(pwd -P)" --use-conda --conda-frontend mamba output/mitoHiFi/PR01227/fromContig/joana_settings_shortcut/no_opts/final_mitogenome.fasta output/mitoHiFi/PR01227/by_reads/joana_settings_shortcut/final_mitogenome.fasta  
```



# HiFi reads filtering outputs for individuals sequenced with Sequel II still containing HiFi adaptaters


![image](https://user-images.githubusercontent.com/42983167/235263044-930cc524-78e6-476f-89cb-a349adf2c680.png)

