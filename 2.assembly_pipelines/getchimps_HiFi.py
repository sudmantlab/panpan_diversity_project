import pandas as pd 
#samples = pd.read_csv('haib22PS6954.2022_05_18.urlsB.csv')
#samples = pd.read_csv('haib22PS6954.2022_08_23.urlsB.csv')
#samples = pd.read_csv('haib22PS6954.2022_08_30.urls.csv')
#samples = pd.read_csv('haib22PS6954.2022_09_06.urls.csv')
#samples = pd.read_csv('haib22PS6954.2022_09_06.urls.repeat.csv')
#samples = pd.read_csv('haib22PS6954.2022_09_12.urls.csv')
#samples = pd.read_csv('haib22PS6954.2022_09_26.urls.csv')
#samples = pd.read_csv('haib22PS6954.2022_10_03.urls.csv')
#samples = pd.read_csv('haib22PS6954.2022_10_17.urls.csv')
#samples = pd.read_csv('haib22PS6954.2023_01_10.urls.csv')
#samples = pd.read_csv('haib22PS6954.2023_01_19.urls.csv')
#samples = pd.read_csv('haib22PS6954.2022_09_06.urls.20230125.csv')
samples = pd.read_csv('haib22PS6954.2023_02_02.urls.csv') #Feb 3 (last batcg)
print(samples)


rule all:
    input:
        expand('PANPAN_ccs/{values[0]}/new/hudsonalpha/l{values[1]}/{values[2]}.l{values[1]}.ccs.bam', values=samples[["chimpname", "lane", "sampleID"]].values),
        
### Feb 3 changed true and False in params because it was switching pbi and bams outnames
rule download_bam:

    output:
        directory('PANPAN_ccs/{chimpname}/')
    params:
        lambda wildcards: samples.loc[(wildcards.sampleID, wildcards.lane, False), "url"]
    shell:
        'wget "{params}" -O {output}'

