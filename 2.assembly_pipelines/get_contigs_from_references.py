
#genomes = ["GCF_002880755.1_Clint_PTRv2_genomic", "GCA_013052645.3_Mhudiblu_PPA_v2", "GCF_009914755.1_T2T-CHM13v2.0_genomic", "GCF_000001405.40_GRCh38.p14_genomic"]
genomes = ["GCF_000001405.40_GRCh38.p14_genomic"]

rule all:
    input:
        expand("PANPAN_graph/hifiasm-fasta/{genome}.p_ctg.fa", genome = genomes)


#rule samtools_faidx:
 #   input: "data/genomes/{genome}.fa"
  #  output: "data/genomes/{genome}.fa.fai"
   # shell: "samtools faidx {input}"

rule genome_file:
    #input: "reference/{genome}.fna.fai"
    input: "reference/{genome}.fa.fai"
    output: "reference/{genome}.genome"
    shell: "cut -f1,2 {input} > {output}"

#bedtools sort -i GCA_013052645.3_Mhudiblu_PPA_v2_gaps.bed -faidx GCA_013052645.3_Mhudiblu_PPA_v2.fna.fai > GCA_013052645.3_Mhudiblu_PPA_v2_gaps.sorted.bed 

rule bedtools_complement:
    input: 
        chrom = "reference/{genome}.genome",
        gap = "reference/{genome}_gaps.sorted.bed"
    output: "reference/{genome}.contigs.bed"
    shell: "bedtools complement -g {input.chrom} -i {input.gap} > {output}"

rule get_fasta_contigs:
    input: 
        fa = "reference/{genome}.fa",
        fai = "reference/{genome}.fa.fai",
        bed = "reference/{genome}.contigs.bed"
    output: "reference/{genome}.p_ctg.fa"
    shell: "bedtools getfasta -name+ -fullHeader -fi {input.fa} -bed {input.bed} -fo {output}"

rule link_fasta_contigs:
    input: "reference/{genome}.p_ctg.fa"
    output: "PANPAN_graph/hifiasm-fasta/{genome}.p_ctg.fa"
    shell: "PROJDIR=$(pwd -P); ln -sf $PROJDIR/{input} {output}"
