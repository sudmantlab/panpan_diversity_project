import numpy as np

#ref_path = '/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fna'
#ref_path = '/global/scratch/users/joana_rocha/PANPAN/reference/GCA_013052645.3_Mhudiblu_PPA_v2.fna'
ref_path = '/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta'
#ref_path= "/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanPan1/mPanPan1.pri.cur.20231122.fasta"
#make sure you ran the following comands on ref_path:
#bwa index /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta
#samtools faidx /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta
#awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta.fai > /global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.bed
#java -jar /global/scratch/users/joana_rocha/software/picard-2.27.2/picard.jar CreateSequenceDictionary R=/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.fasta O=/global/scratch/users/joana_rocha/PANPAN/reference/primates_T2T/mPanTro3/mPanTro3.pri.cur.20231031.dict



accessions = {
        'SAMN01920525_ellioti_Taweh_Cameroon' : ['SRR748156','SRR748157','SRR748158','SRR748159','SRR748160'],
        'SAMN01920523_ellioti_Koto_Cameroon': ['SRR748147','SRR748148','SRR748149','SRR748150','SRR748151'],
        'SAMN01920521_ellioti_JulieLWC21_Cameroon' : ['SRR748138','SRR748139','SRR748140','SRR748141','SRR748142'],
        'SAMN01920520_ellioti_Damian_Cameroon' : ['SRR748133','SRR748134','SRR748135','SRR748136','SRR748137'],
        'SAMN01920517_elliotiAkwaya-Jean_Cameroon' : ['SRR748121', 'SRR748122', 'SRR748123', 'SRR748124'],
        'SAMN01920526_ellioti_Tobi_Cameroon' : ['SRR748161','SRR748162','SRR748163','SRR748164'],
        'SAMN01920522_ellioti_Kopongo_Cameroon' : ['SRR748143','SRR748144','SRR748145','SRR748146'], 
        'SAMN01920519_ellioti_Basho_Cameroon' :  ['SRR748129','SRR748130','SRR748131','SRR748132'],
        'SAMN01920518_ellioti_Banyo_Cameroon' : ['SRR748125','SRR748126','SRR748127','SRR748128'], 
        'SAMN01920524_ellioti_Paquita_Cameroon' : ['SRR748152','SRR748153','SRR748154','SRR748155'],
        'SAMN01920532_schweinfurthii_Nakuu_DemocraticRepublicCongo' : ['SRR726408','SRR726409','SRR726410','SRR726411','SRR726412','SRR726413','SRR726415','SRR726416'],
        'SAMN01920531_schweinfurthii_Kidongo_DemocraticRepublicCongo' : ['SRR726360','SRR726381','SRR726395','SRR726402','SRR726403','SRR726404','SRR726405','SRR726406','SRR726407'],
        'SAMN01920530_schweinfurthii_Bwambale_Uganda' : ['SRR726352','SRR726353','SRR726354','SRR726355','SRR726356','SRR726357','SRR726358','SRR726359'],
        'SAMN01920528_schweinfurthii_Andromeda_Tanzania' : ['SRR726275','SRR726278','SRR726279','SRR726280'],
        'SAMN01920527_schweinfurthii_Vincent_Tanzania' : ['SRR726233','SRR726241','SRR726242','SRR726243'],
        'SAMN01920529_schweinfurthii_Harriet_Uganda' : ['SRR747955','SRR747956','SRR747957','SRR747958'],
        'SAMEA4374799_schweinfurthii_DemocraticRepublicCongo' : ['ERR1709873', 'ERR1709874', 'ERR1709875', 'ERR1709876', 'ERR1709877', 'ERR1709878'],
        'SAMEA4374798_schweinfurthii_Washu_DemocraticRepublicCongo' : ['ERR1710113','ERR1710114','ERR1710115','ERR1710116'],
        'SAMEA4374797_schweinfurthii_Athanga_DemocraticRepublicCongo' : ['ERR1709883', 'ERR1709884'],
        'SAMEA4374796_schweinfurthii_Coco_Zambia' : ['ERR1710107', 'ERR1710108'],
        'SAMEA4374795_schweinfurthii_Maya_DemocraticRepublicCongo' : ['ERR1710109','ERR1710110','ERR1710111','ERR1710112'],
        'SAMEA4374794_schweinfurthii_Frederike_Rwanda' : ['ERR1710105', 'ERR1710106'],
        'SAMEA4374793_schweinfurthii_Trixie_DemocraticRepublicCongo' : ['ERR1710117', 'ERR1710118', 'ERR1710119', 'ERR1710120'],
        'SAMEA4374792_schweinfurthii_Bihati' : ['ERR1710099','ERR1710100','ERR1710101','ERR1710102'],
        'SAMEA4374791_schweinfurthii_Cleo_DemocraticRepublicCongo' : ['ERR1710024','ERR1710025','ERR1710026','ERR1710027','ERR1710028','ERR1710029'],
        'SAMEA4374790_schweinfurthii_Tongo_DemocraticRepublicCongo' : ['ERR1710095','ERR1710096','ERR1710097','ERR1710098'],
        'SAMEA4374789_schweinfurthii_Ikuru_DemocraticRepublicCongo' : ['ERR1710149','ERR1710150','ERR1710151','ERR1710152'],
        'SAMEA4374788_schweinfurthii_Cindy_Uganda' : ['ERR1710085','ERR1710086','ERR1710087','ERR1710088'],
        'SAMEA4374787_schweinfurthii_Padda_DemocraticRepublicCongo' : ['ERR1709990','ERR1709991','ERR1709992','ERR1709993'],
        'SAMN01920536_troglodytes_Clara_Gabon' : ['SRR748071', 'SRR748072','SRR748073','SRR748074','SRR748075'],
        'SAMN01920535_troglodytes_Julie_Gabon' : ['SRR748066', 'SRR748067','SRR748068','SRR748069','SRR748070'],
        'SAMN01920534_troglodytes_Doris_Gabon' : ['SRR748061', 'SRR748062', 'SRR748063','SRR748064','SRR748065'],
        'SAMN01920533_troglodytes_Vaillant_Gabon' : ['SRR748056','SRR748057','SRR748058','SRR748059','SRR748060'],
        'SAMEA4374785_troglodytes_Alfred' : ['ERR1710061','ERR1710062','ERR1710063','ERR1710064','ERR1710065','ERR1710066','ERR1710067','ERR1710068','ERR1710069','ERR1710070','ERR1710071','ERR1710072','ERR1710073','ERR1710074','ERR1710075','ERR1710076','ERR1710077','ERR1710078','ERR1710079','ERR1710080','ERR1710081','ERR1710082','ERR1710083','ERR1710084'],
        'SAMEA4374782_troglodytes_Brigitta' : ['ERR1709936','ERR1709937','ERR1709938','ERR1709939','ERR1709940','ERR1709941','ERR1709942','ERR1709943','ERR1709944','ERR1709945','ERR1709946','ERR1709947','ERR1709948','ERR1709949','ERR1709950','ERR1709951','ERR1709952','ERR1709953','ERR1709954','ERR1709955','ERR1709956','ERR1709957','ERR1709958','ERR1709959'],
        'SAMEA4374781_troglodytes_Gamin' : ['ERR1709966','ERR1709967','ERR1709968','ERR1709969','ERR1709970','ERR1709971','ERR1709972','ERR1709973','ERR1709974','ERR1709975','ERR1709976','ERR1709977','ERR1709978','ERR1709979','ERR1709980','ERR1709981','ERR1709982','ERR1709983','ERR1709984','ERR1709985','ERR1709986','ERR1709987','ERR1709988','ERR1709989'],
        'SAMEA4374780_troglodytes_Luky_EquatorialGuinea' : ['ERR1709912','ERR1709913','ERR1709914','ERR1709915','ERR1709916','ERR1709917','ERR1709918','ERR1709919','ERR1709920','ERR1709921','ERR1709922','ERR1709923','ERR1709924','ERR1709925','ERR1709926','ERR1709927','ERR1709928','ERR1709929','ERR1709930','ERR1709931','ERR1709932','ERR1709933','ERR1709934','ERR1709935'],
        'SAMEA4374779_troglodytes_Lara_EquatorialGuinea' : ['ERR1710125','ERR1710126','ERR1710127','ERR1710128','ERR1710129','ERR1710130','ERR1710131','ERR1710132','ERR1710133','ERR1710134','ERR1710135','ERR1710136','ERR1710137','ERR1710138','ERR1710139','ERR1710140','ERR1710141','ERR1710142','ERR1710143','ERR1710144','ERR1710145','ERR1710146','ERR1710147','ERR1710148'],
        'SAMEA4374778_troglodytes_Ula_EquatorialGuinea' : ['ERR1709885','ERR1709886','ERR1709887','ERR1709888','ERR1709889','ERR1709890','ERR1709891','ERR1709892','ERR1709893','ERR1709894','ERR1709895','ERR1709896','ERR1709897','ERR1709898','ERR1709899','ERR1709900','ERR1709901','ERR1709902','ERR1709903','ERR1709904','ERR1709905','ERR1709906','ERR1709907','ERR1709908'],
        'SAMEA4374777_troglodytes_Mirinda' : ['ERR1709994','ERR1709995','ERR1709996','ERR1709997','ERR1709998','ERR1709999','ERR1710000','ERR1710001','ERR1710002','ERR1710003','ERR1710004','ERR1710005','ERR1710006','ERR1710007','ERR1710008','ERR1710009','ERR1710010','ERR1710011','ERR1710012','ERR1710013','ERR1710014','ERR1710015','ERR1710016','ERR1710017'],
        'SAMEA4374776_troglodytes_Cindy' : ['ERR1710038','ERR1710039','ERR1710040','ERR1710041','ERR1710042','ERR1710043','ERR1710044','ERR1710045','ERR1710046','ERR1710047','ERR1710048','ERR1710049','ERR1710050','ERR1710051','ERR1710052','ERR1710053','ERR1710054','ERR1710055','ERR1710056','ERR1710057','ERR1710058','ERR1710059','ERR1710060'],
        'SAMEA4374775_troglodytes_Noemie_EquatorialGuinea' : ['ERR1709909','ERR1709910','ERR1709911'],
        'SAMEA4374774_troglodytes_Yogui_EquatorialGuinea' : ['ERR1710103','ERR1710104'],
        'SAMEA4374773_troglodytes_Tibe_EquatorialGuinea' : ['ERR1709871', 'ERR1709872'],
        'SAMEA4374772_troglodytes_Blanquita_EquatorialGuinea' : ['ERR1709960','ERR1709961'],
        'SAMEA4374771_troglodytes_Negrita_EquatorialGuinea' : ['ERR1710093','ERR1710094'],
        'SAMEA4374770_troglodytes_Marlin' : ['ERR1710018', 'ERR1710019'],
        'SAMN01920538_verus_Jimmie' : ['SRR748051','SRR748052','SRR748053','SRR748054','SRR748055'],
        'SAMN01920537_verus_Bosco' : ['SRR747951','SRR747952','SRR747953','SRR747954'],
        'SAMN01920541_verusxtroglodytes_Donald' : ['SRR747947','SRR747948','SRR747949','SRR747950'],
        'SAMN01920540_verus_Clint' : ['SRR748179','SRR748180','SRR748181','SRR748182','SRR748183','SRR748184','SRR748185','SRR748186'],
        'SAMN01920539_verus_Koby' : ['SRR748081', 'SRR748082','SRR748083','SRR748084'],
        'SAMEA4374769_verus_Alice_CotedIvoire' : ['ERR1710034','ERR1710035','ERR1710036','ERR1710037'],  
        'SAMEA4374768_verus_Cindy_CotedIvoire' : ['ERR1710121','ERR1710122','ERR1710123','ERR1710124'],      
        'SAMEA4374767_verus_Linda_CotedIvoire' : ['ERR1710020','ERR1710021','ERR1710022','ERR1710023'],
        'SAMEA4374766_verus_SepToni_Liberia' : ['ERR1709879', 'ERR1709880', 'ERR1709881', 'ERR1709882'],
        'SAMEA4374765_verus_Mike_Guinea' : ['ERR1710030','ERR1710031','ERR1710032','ERR1710033'],
        'SAMEA4374764_verus_Annie_Guinea' : ['ERR1710089','ERR1710090','ERR1710091','ERR1710092'],
        'SAMEA4374763_verus_Berta_CotedIvoire' : ['ERR1709962', 'ERR1709963', 'ERR1709964','ERR1709965'],
        'SAMN01920505_Hortense_paniscus': ['SRR726612', 'SRR726613', 'SRR726614', 'SRR726615', 'SRR726607', 'SRR726608', 'SRR726609', 'SRR726610'],
        'SAMN01920506_Kosana_paniscus' : ['SRR740768', 'SRR740769', 'SRR740770', 'SRR740771', 'SRR740772', 'SRR740773', 'SRR740774', 'SRR740775'],
        'SAMN01920507_Dzeeta_paniscus': ['SRR740781', 'SRR740782', 'SRR740783', 'SRR740784','SRR740787', 'SRR740790', 'SRR740792', 'SRR740793'],
        'SAMN01920508_Hermien_paniscus': ['SRR740794', 'SRR740795', 'SRR740796','SRR740797','SRR740798', 'SRR740799', 'SRR740800', 'SRR740801'],
        'SAMN01920509_Desmond_paniscus': ['SRR740802', 'SRR740803', 'SRR740804', 'SRR740805', 'SRR740806', 'SRR740807', 'SRR740808', 'SRR740809'],
        'SAMN01920510_Catherine_paniscus': ['SRR740816', 'SRR740817', 'SRR740818','SRR740819', 'SRR740820', 'SRR740821'],
        'SAMN01920511_Kombote_paniscus': ['SRR740822','SRR740823', 'SRR740824','SRR740825', 'SRR740827', 'SRR740828'],
        'SAMN01920512_Chipita_paniscus': ['SRR740831', 'SRR740832', 'SRR740833', 'SRR740834', 'SRR740835'],
        'SAMN01920513_Bono_paniscus': ['SRR740980', 'SRR740941', 'SRR740905' , 'SRR740911' , 'SRR740853', 'SRR740857'], 
        'SAMN01920514_Natalie_paniscus' : ['SRR741390', 'SRR741327', 'SRR741254', 'SRR741276', 'SRR741194', 'SRR741205'],
        'SAMN01920515_Salonga_paniscus' : ['SRR741813', 'SRR741824', 'SRR741768', 'SRR741770', 'SRR741785'],
        'SAMN01920516_paniscus_Kumbuka' : ['SRR747643','SRR747644','SRR747645','SRR747646','SRR747647','SRR747648'],
        'SAMN01920504_LB502_paniscus' : ['SRR747929', 'SRR747930', 'SRR747927', 'SRR747928'],
        }


accessions_ppa = {
        'SAMN01920505_Hortense_paniscus': ['SRR726612', 'SRR726613', 'SRR726614', 'SRR726615', 'SRR726607', 'SRR726608', 'SRR726609', 'SRR726610'],
        'SAMN01920506_Kosana_paniscus' : ['SRR740768', 'SRR740769', 'SRR740770', 'SRR740771', 'SRR740772', 'SRR740773', 'SRR740774', 'SRR740775'],
        'SAMN01920507_Dzeeta_paniscus': ['SRR740781', 'SRR740782', 'SRR740783', 'SRR740784','SRR740787', 'SRR740790', 'SRR740792', 'SRR740793'],
        'SAMN01920508_Hermien_paniscus': ['SRR740794', 'SRR740795', 'SRR740796','SRR740797','SRR740798', 'SRR740799', 'SRR740800', 'SRR740801'],
        'SAMN01920509_Desmond_paniscus': ['SRR740802', 'SRR740803', 'SRR740804', 'SRR740805', 'SRR740806', 'SRR740807', 'SRR740808', 'SRR740809'],
        'SAMN01920510_Catherine_paniscus': ['SRR740816', 'SRR740817', 'SRR740818','SRR740819', 'SRR740820', 'SRR740821'],
        'SAMN01920511_Kombote_paniscus': ['SRR740822','SRR740823', 'SRR740824','SRR740825', 'SRR740827', 'SRR740828'],
        'SAMN01920512_Chipita_paniscus': ['SRR740831', 'SRR740832', 'SRR740833', 'SRR740834', 'SRR740835'],
        'SAMN01920513_Bono_paniscus': ['SRR740980', 'SRR740941', 'SRR740905' , 'SRR740911' , 'SRR740853', 'SRR740857'], 
        'SAMN01920514_Natalie_paniscus' : ['SRR741390', 'SRR741327', 'SRR741254', 'SRR741276', 'SRR741194', 'SRR741205'],
        'SAMN01920515_Salonga_paniscus' : ['SRR741813', 'SRR741824', 'SRR741768', 'SRR741770', 'SRR741785'],
        'SAMN01920516_paniscus_Kumbuka' : ['SRR747643','SRR747644','SRR747645','SRR747646','SRR747647','SRR747648'],
        'SAMN01920504_LB502_paniscus' : ['SRR747929', 'SRR747930', 'SRR747927', 'SRR747928'],
        }


samples = accessions.keys()
samples_ppa = accessions_ppa.keys()


rule all:
    input:
        #expand('bam_final_mapped2mPanPan1/{sample}.sorted.bam', sample=samples_ppa),
        #'coverage_stats/bam_final_mapped2mPanPan1_coverage_plot.pdf'
        'coverage_stats/bam_final_mapped2mPanTro3_coverage_plot.pdf'
        #expand('bams/{sample}.flagfilt.sorted.bam', sample=samples_ppa),
        #expand('bams/{sample}.flagfilt.rmdup.sorted.rg.bam', sample=samples)
        #expand('fastqs/{sample}_1.fastq.gz', sample=samples),
        #expand('fastqs/{sample}_2.fastq.gz', sample=samples),
        #expand('fastqs/{sample}_1.trim1P.fastq.gz', sample=samples),
        #expand('fastqs/{sample}_2.trim2P.fastq.gz', sample=samples)
    
rule prefetch_samples:
    output: 
        temp('accessions/{accession}_1.fastq.gz'),
        temp('accessions/{accession}_2.fastq.gz')
    shell: """
    cd accessions && ../ffq --ftp {wildcards.accession} | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs -I{{}} sh -c 'curl -s -L -o $(basename "{{}}" ) "{{}}"'
#    """
#cd accessions && ../ffq --ftp {wildcards.accession} | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O 
   
def get_input_merge_fastq(wildcards):
    return expand('accessions/{{accession}}_{number}.fastq.gz'.format(number=wildcards.number),
            accession=accessions[wildcards.sample])

rule merge_fastq:
    input: get_input_merge_fastq
    output: 'fastqs/{sample}_{number}.fastq.gz'
    shell: 'cat {input} > {output}'

#rule cut_adapt:
 #   input:
  #      'fastqs/{sample}_1.fastq.gz',
   #     'fastqs/{sample}_2.fastq.gz'
    #output:
     #   temp('fastqs/{sample}_1.trim.fastq.gz'),
      #  temp('fastqs/{sample}_2.trim.fastq.gz')
    #shell: 'cutadapt  -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -n 8 -e 0.1 -O 1 -m 30 -q 20,20 --max-n 0.5 --pair-filter any --format=sra-fast -o {output[0]} -p {output[1]} {input[0]} {input[1]}'


rule trim_clean:
    input:
        'fastqs/{sample}_1.fastq.gz',
        'fastqs/{sample}_2.fastq.gz'
    output:
       'fastqs/{sample}_1.trim1P.fastq.gz',
        temp('fastqs/{sample}_1.trim1U.fastq.gz'),
        'fastqs/{sample}_2.trim2P.fastq.gz',
        temp('fastqs/{sample}_2.trim2U.fastq.gz'),
    shell:'java -jar /global/scratch/users/joana_rocha/software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -quiet -threads 8 -phred33 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} SLIDINGWINDOW:4:20 MINLEN:30'

rule map_filt:
    input:
        'fastqs/{sample}_1.trim1P.fastq.gz',
        'fastqs/{sample}_2.trim2P.fastq.gz'
    output: temp('bams/{sample}.flagfilt.bam')
    shell: 'bwa mem -t 20 {ref_path} {input[0]} {input[1]} | samtools view -q 15 -bT {ref_path} -F 780 -o {output}'

rule map_filt_sort:
    input: 'bams/{sample}.flagfilt.bam'
    output: 'bams/{sample}.flagfilt.sorted.bam'
    shell: 'samtools sort -o {output} {input}'

rule picard_rmdup:
    input: 'bams/{sample}.flagfilt.sorted.bam'
    output: 'bams/{sample}.flagfilt.rmdup.bam'
    shell: """
    mkdir bams/tmp_0_{wildcards.sample};
    java -Xmx31G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=bams/tmp_0_{wildcards.sample} -jar /global/scratch/users/joana_rocha/software/picard-2.27.2/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I={input}  O={output} M=bams/{wildcards.sample}.rmdup.metrics
    """

rule rmdup_sort:
    input: 'bams/{sample}.flagfilt.rmdup.bam'
    output: 'bams/{sample}.flagfilt.rmdup.sorted.bam'
    shell: 'samtools sort -o {output} {input}'

rule add_readgroup:
    input: 'bams/{sample}.flagfilt.rmdup.sorted.bam'
    output: 'bams/{sample}.flagfilt.rmdup.sorted.rg.bam'
    shell: 'samtools addreplacerg -r ID:S1 -r SM:{wildcards.sample} -r PU:L1 -r PL:ILLUMINA -o {output} {input}'

rule indel_remap:
    input: 'bams/{sample}.flagfilt.rmdup.sorted.rg.bam'
    output:
        'bam_final/{sample}.flagfilt.rmdup.sorted.rg.intervals',
        'bam_final/{sample}.flagfilt.rmdup.sorted.realign.rg.bam'
    shell: """
    samtools index {input} &&
    mkdir bam_final/tmp_0_{wildcards.sample} &&
    /global/scratch/users/joana_rocha/software/jre1.8.0_211/bin/java -Xmx5G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_final/tmp_0_{wildcards.sample} -jar /global/scratch/users/joana_rocha/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {ref_path} -I {input} -o {output[0]} &&
    /global/scratch/users/joana_rocha/software/jre1.8.0_211/bin/java -Xmx4g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_final/tmp_0_{wildcards.sample} -jar /global/scratch/users/joana_rocha/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T IndelRealigner -I {input} -R {ref_path} -targetIntervals {output[0]} -o {output[1]}
    """


#  Create final files for snpCleaner
rule bam_final:
    input: 'bam_final/{sample}.flagfilt.rmdup.sorted.realign.rg.bam',
    output: 'bam_final_mapped2mPanTro3/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """

rule plotcoverage:
    input: 
        bams = expand('bam_final_mapped2mPanTro3/{sample}.sorted.bam', sample=samples)
    output: 
        plot = 'coverage_stats/bam_final_mapped2mPanTro3_coverage_plot.pdf'
    params:
        bam_list = lambda wildcards, input: ' '.join(input.bams)
    shell: """
    plotCoverage --bamfiles {params.bam_list} --plotFile {output.plot} -n 1000000 --plotTitle "Chimpanzees (short-reads) mapped to mPanTro3" --outRawCounts coverage.tab --ignoreDuplicates --plotFileFormat pdf --plotHeight 20 --plotWidth 25 --minMappingQuality 15 -p 24 
    """