#!/usr/bin/env python3

# conda activate fasta
# want to use the latest version of IntervalTree (pip install has no merge_neighbours)
#   pip install git+https://github.com/chaimleib/intervaltree.git@01a30a8df977a02e1b6e25eee3838b8d6b340995

from Bio import SeqIO
import pandas as pd
from intervaltree import Interval, IntervalTree
import sys

if len(sys.argv) != 4:
    print("{} syri.out saveFile refName".format(sys.argv[0]))
    exit()

syriFile = sys.argv[1]
saveFile = sys.argv[2]
refName  = sys.argv[3] # Added for vcf header
saveFile = saveFile.replace('.vcf','')
sample = saveFile.split('/')[-1]

print("syriFile = {}".format(syriFile))
print("saveFile = {}.vcf".format(saveFile))

low_memory=False

syriHeader = ['refChr','refStart','refStop','1','2','qryChr','qryStart','qryStop','uid','parentId','type','dupType']
syriTypes = {'refChr'    :'str',
            'refStart'  :'str',
            'refStop'   :'str',
            '1'         :'str',
            '2'         :'str',
            'qryChr'    :'str',
            'qryStart'  :'str',
            'qryStop'   :'str',
            'uid'       :'str',
            'parentId'  :'str',
            'type'      :'str',
            'dupType'   :'str'}
minSyriSize = 10000
mergeDist   = 1000

def dupMerge(a, b):
    return f"{a}.{b}" 

def mergeIntervals(df, refQry):
    treeDict = {}
    if refQry == "ref":
        for index, row in df.iterrows():
            if row.refChr == "-":
                continue
            if row.refChr not in treeDict:
                treeDict[row.refChr] = IntervalTree()
            tree = treeDict[row.refChr]
            tree.addi(int(row.refStart), int(row.refStop), row.uid)
    else:
        return None
    newDf = []
    for tree in treeDict:
      treeDict[tree].merge_overlaps(dupMerge)
      treeDict[tree].merge_neighbors(data_reducer=dupMerge, distance=mergeDist)
      for interval in treeDict[tree]:
          newDf.append({
              'chr'  : tree,
              'start': interval.begin,
              'end'  : interval.end,
              'uid'  : interval.data})
    if not newDf:
        return pd.DataFrame(columns=['chr', 'start', 'end', 'uid'])
    newDf = pd.DataFrame(newDf)
    newDf.reset_index(inplace=True)
    newDf = newDf[['chr', 'start', 'end', 'uid']]
    return newDf

def writeVcfHeader(file, sample):
    file.writelines([
        "##fileformat=VCFv4.2\n",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
        "##source=SF_from_SyRI\n",
        "##reference=", refName, "\n",
        "##ALT=<ID=DEL,Description=\"Deletion\">\n",
        "##ALT=<ID=DUP,Description=\"Duplication, Region of elevated or reduced copy number\">\n",
        "##ALT=<ID=INV,Description=\"Inversion\">\n",
        "##ALT=<ID=TRANS,Description=\"Translocation\">\n",
        "##ALT=<ID=HDR,Description=\"Homologous Deletion Region information\">\n",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n",
        "##INFO=<ID=SyRI,Number=1,Type=String,Description=\"SyRI annotation category\">\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample, "\n"])
    return None

# read in syri
syriDf = pd.read_csv(syriFile, sep="\t",names=syriHeader, dtype=syriTypes)
syriDf.type = syriDf.type.replace(['INVTR', 'INVDP'], ["TRANS", "DUP"])

vcfRecords = []

# Gather non-duplications (Skipping reference-absent insertions)
for syri in syriDf[syriDf.type != "DUP"].itertuples(index=False):
    if syri.refChr == "-":
        continue
        
    length = int(syri.refStop) - int(syri.refStart) + 1
    if length < minSyriSize:
        稳continue

    if syri.type == "INV":
        vcfRecords.append([syri.refChr, int(syri.refStart), sample+"~"+syri.uid, "N", "<INV>", "60", "PASS", "SVTYPE=INV;SVLEN="+str(length)+";END="+str(syri.refStop)+";SyRI=INV", "GT", "1/1"])
    if syri.type == "TRANS":
        vcfRecords.append([syri.refChr, int(syri.refStart), sample+"~"+syri.uid, "N", "<TRANS>", "60", "PASS", "SVTYPE=TRANS;SVLEN="+str(length)+";END="+str(syri.refStop)+";SyRI=TRANS", "GT", "1/1"])
    if syri.type == "NOTAL":
        if syri.qryChr == "-": # DEL
            vcfRecords.append([syri.refChr, int(syri.refStart), sample+"~"+syri.uid, "N", "<DEL>", "60", "PASS", "SVTYPE=DEL;SVLEN=-"+str(length)+";END="+str(syri.refStop)+";SyRI=NOTAL.ref", "GT", "1/1"])
    if syri.type == "HDR":
        vcfRecords.append([syri.refChr, int(syri.refStart), sample+"~"+syri.uid, "N", "<HDR>", "60", "PASS", "SVTYPE=HDR;SVLEN=-"+str(length)+";END="+str(syri.refStop)+";SyRI=HDR", "GT", "1/1"])

# Gather duplications - copyloss
dupsDf = syriDf[(syriDf.type=="DUP") & (syriDf.dupType=="copyloss")] 
for index, dup in mergeIntervals(dupsDf, 'ref').iterrows():    
    if (dup.end-dup.start) < minSyriSize:
        continue
    vcfRecords.append([dup.chr, int(dup.start), sample+"~"+dup.uid, "N", "<DUP>", "60", "PASS", "SVTYPE=DUP;SVLEN="+str(dup.end-dup.start)+";END="+str(dup.end)+";SyRI=DUP:copyloss", "GT", "1/1"])

# sort and save vcf
if vcfRecords:
    df = pd.DataFrame(vcfRecords, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample])
    df['POS'] = df['POS'].astype(int)
    df = df.sort_values(by=['CHROM', 'POS']).reset_index(drop=True)
    
    with open(saveFile+".vcf", "w") as file:
        writeVcfHeader(file, sample)
        for index, row in df.iterrows():
            file.writelines([str(row['CHROM']), "\t", str(row['POS']), "\t", str(row['ID']), "\t", str(row['REF']), "\t", str(row['ALT']), "\t", str(row['QUAL']), "\t", str(row['FILTER']), "\t", str(row['INFO']), "\t", str(row['FORMAT']), "\t", str(row[sample]), "\n"])
