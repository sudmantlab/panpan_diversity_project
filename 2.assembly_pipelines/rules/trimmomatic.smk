rule trimmomatic_pe_HiC:
    input:
        #r1="PANPAN_hic/{specimen}/{sample_name}_1.fq.gz",
        #r2="PANPAN_hic/{specimen}/{sample_name}_2.fq.gz"
        r1="PANPAN_hic/{specimen}/{sample_name}_R1.fastq.gz", # for AG chimps
        r2="PANPAN_hic/{specimen}/{sample_name}_R2.fastq.gz" # for AG chimps
    output:
        r1="output/trimmed-hic/{specimen}/{sample_name}-trimmed_R1.fastq.gz",
        r2="output/trimmed-hic/{specimen}/{sample_name}-trimmed_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{specimen}/{sample_name}-trimmed_R1_unpaired.fastq.gz",
        r2_unpaired="output/trimmed/{specimen}/{sample_name}-trimmed_R2_unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{specimen}/{sample_name}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:PANPAN_hic/adapters/TruSeq2-PE.fa:2:40:15", ###  NexteraPE-PE.fa for novagene. this can vary!
                 "SLIDINGWINDOW:5:20"
                ],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 32
    wrapper:
        "0.74.0/bio/trimmomatic/pe"
