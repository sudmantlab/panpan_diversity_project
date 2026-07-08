#!/usr/bin/env python3

#
# script to find all telomeres in a collection of genomes
# Outpus a working file with telo lengths for all found telomeres, and a summary file per sppecies
# Summary file contains numner of telomeres at the begining of chromosomes, end of chromosomes, and number of chromosomes with both telomeres present
#
# Requires seqtk
#

import subprocess
import pandas as pd
import io
import glob as gb

fnaDir = "/XXX/YYY"

first = True
sex = False  # True if you want X & Y counts
teloMotif = "CCCTAA" # primate motif

results_list = []
summary_list = []
chromDf = None

for fnaFile in gb.glob(fnaDir + '/*.fasta'):
    sample = fnaFile.split('/')[-1]
    sample = sample.replace(".chr.fasta", "")

    if first:  # get list of chroms for species
        with open(fnaFile) as f:
            chroms = [line[1:].split()[0] for line in f if line.startswith('>')]
        chromDf = pd.DataFrame({"chrom": chroms})
        if not sex:
            chromDf = chromDf[~chromDf['chrom'].str.contains("chrX")]
            chromDf = chromDf[~chromDf['chrom'].str.contains("chrY")]
        first = False

    result = subprocess.run(['seqtk', 'telo', '-m', teloMotif, fnaFile], stdout=subprocess.PIPE, text=True)
    data = io.StringIO(result.stdout)
    df = pd.read_csv(data, sep='\t', header=None)
    df.columns = ["chrom", "start", "end", "chromLength"]

    tmp1 = df.loc[df.start == 0, ['chrom', 'end']].copy()
    tmp1.columns = ['chrom', 'startLength']

    tmp2 = df.loc[df.start != 0].copy()
    tmp2['len'] = tmp2['chromLength'] - tmp2['start']
    tmp2 = tmp2[['chrom', 'len']]
    tmp2.columns = ['chrom', 'endLength']

    mergedDf = pd.merge(tmp1, chromDf, on='chrom', how='right')
    mergedDf = pd.merge(tmp2, mergedDf, on='chrom', how='right')
    mergedDf.fillna(0, inplace=True)
    mergedDf['startLength'] = mergedDf['startLength'].astype(int)
    mergedDf['endLength'] = mergedDf['endLength'].astype(int)
    mergedDf["sample"] = sample

    results_list.append(mergedDf)
    summary_list.append({
        'sample': sample,
        'teloBegin': int((mergedDf['startLength'] != 0).sum()),
        'teloEnd': int((mergedDf['endLength'] != 0).sum()),
        'T2T': int(((mergedDf['startLength'] != 0) & (mergedDf['endLength'] != 0)).sum()),
    })

resultsDf = pd.concat(results_list, ignore_index=True)
summaryDf = pd.DataFrame(summary_list)

saveName = fnaDir.split('/')[-1]
resultsDf.to_csv(saveName + ".csv", index=False)
summaryDf.to_csv(saveName + "~summary.csv", index=False)


