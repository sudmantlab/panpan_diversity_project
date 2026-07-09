#!/usr/bin/env python3

# conda activate fasta
# want to use the latest version of IntervalTree (pip install has no merge_neighbours)
#   pip install git+https://github.com/chaimleib/intervaltree.git@01a30a8df977a02e1b6e25eee3838b8d6b340995


from Bio import SeqIO
import pandas as pd
import pafpy as PafFile
from intervaltree import Interval, IntervalTree
import sys

len(sys.argv)
if len(sys.argv) != 6:
    print("{} refFasta qryFasta paf syri.out saveFile".format(sys.argv[0]))
    exit()

refFile  = sys.argv[1]
qryFile  = sys.argv[2]
pafFile  = sys.argv[3]
syriFile = sys.argv[4]
saveFile = sys.argv[5]
saveFile = saveFile.replace('.vcf','')
sample = saveFile.split('/')[-1]

print("refFile  = {}".format(refFile))
print("qryFile  = {}".format(qryFile))
print("pafFile  = {}".format(pafFile))
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
    return f"{a}.{b}" # ":".join([str(a),str(b)])

def mergeIntervals(df, refQry):
    '''
    ref = copyloss
    qry = copygain
    '''
    treeDict = {}
    if refQry == "ref":
        for index, row in df.iterrows():
            if row.refChr not in treeDict:
                treeDict[row.refChr] = IntervalTree()
            tree = treeDict[row.refChr]
            tree.addi(int(row.refStart), int(row.refStop), row.uid)
    elif refQry == "qry":
        for index, row in df.iterrows():
            if row.qryChr not in treeDict:
                treeDict[row.qryChr] = IntervalTree()
            tree = treeDict[row.qryChr]
            tree.addi(int(row.qryStart), int(row.qryStop), row.uid)
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
    newDf = pd.DataFrame(newDf)
    #newDf['len'] = newDf.end - newDf.start +1
    #newDf = newDf[newDf.len >= minSyriSize]
    newDf.reset_index(inplace=True)
    newDf = newDf[['chr', 'start', 'end', 'uid']]
    return newDf

def writeVcfHeader(file, sample):
    file.writelines([
        "##fileformat=VCFv4.2\n",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
        "##source=SF_from_SyRI\n",
        "##reference=", refFile.split("/")[-1].replace(".fasta",""), "\n",
        lenHeader,
        "##ALT=<ID=DEL,Description=\"Deletion\">\n",
        "##ALT=<ID=DUP,Description=\"Duplication, Region of elevated or reduced copy number\">\n",
        "##ALT=<ID=INV,Description=\"Inversion\">\n",
        "##ALT=<ID=TRANS,Description=\"Translocation\">\n",
        "##ALT=<ID=INS,Description=\"Insertion\">\n",
        "##ALT=<ID=HDR,Description=\"Homologous Deletion Region information\">\n",
        "##INFO=<ID=DEST,Number=0,Type=String,Description=\"Destination chromosome and position\">\n",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n",
        "##INFO=<ID=SyRI,Number=1,Type=String,Description=\"SyRI annotation category\">\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample, "\n"])
    return None

# Read in fasta files
refFasta = {}
for record in SeqIO.parse(refFile, "fasta"):
    refFasta[record.id] = record.seq

qryFasta = {}
for record in SeqIO.parse(qryFile, "fasta"):
    qryFasta[record.id] = record.seq


# read in paf file
pafDict = {}
with PafFile.PafFile(pafFile) as ppp:
    for record in ppp:
        if record.qname not in pafDict:
            pafDict[record.qname] = {}
        a = pafDict[record.qname]
        a[record.qstart] = str(record.tname) + ":" + str(record.tstart)
        a[record.qend]   = str(record.tname) + ":" + str(record.tend)


# read in syri
syriDf = pd.read_csv(syriFile, sep="\t",names=syriHeader, dtype=syriTypes)
syriDf.type = syriDf.type.replace(['INVTR', 'INVDP'], ["TRANS", "DUP"])


### #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMATAG06939_13.hap1
# get ref seq lengths for vcf header
lenHeader=""
for x in refFasta:
    lenHeader = lenHeader + "##contig=<ID=" + x + ",length=" + str(len(refFasta[x])) + ">\n"

## ##fileDate="2024/01/19 12:23:18"
with open(saveFile+".vcf", "w") as file, open(saveFile + "~seq.vcf", "w") as seqFile:
    writeVcfHeader(file, sample)
    writeVcfHeader(seqFile, sample)
    for syri in syriDf[syriDf.type != "DUP"].itertuples(index=False):
        length = int(syri.qryStop if syri.refChr == "-" else syri.refStop) - int(syri.qryStart if syri.refChr == "-" else syri.refStart) + 1
        if length < minSyriSize:
            continue
        if syri.type == "INV":
            file.writelines([syri.refChr,"\t", str(syri.refStart), "\t",sample,"~",syri.uid,"\tN\t", "<INV>", "\t60\tPASS\tSVTYPE=INV;SVLEN=", str(length), ";END=", str(syri.refStop), ";SyRI=INV\tGT\t1/1\n"])
        if syri.type == "TRANS":
            file.writelines([syri.refChr,"\t", str(syri.refStart), "\t",sample,"~",syri.uid,"\tN\t", "<TRANS>", "\t60\tPASS\tSVTYPE=TRANS;SVLEN=", str(length), ";END=", str(syri.refStop), ";SyRI=TRANS;DEST=",syri.qryChr,".",syri.qryStart,"\tGT\t1/1\n"])
        if syri.type == "NOTAL":
            if syri.qryChr == "-": # DEL
                file.writelines([syri.refChr,"\t", str(syri.refStart), "\t",sample,"~",syri.uid,"\tN\t", "<DEL>", "\t60\tPASS\tSVTYPE=DEL;SVLEN=", str(length), ";END=", str(syri.refStop), ";SyRI=NOTAL.ref\tGT\t1/1\n"])
            if syri.refChr == "-": # INS, insert seq, output to seqFile
                seq = str(qryFasta[syri.qryChr][int(syri.qryStart)-1:int(syri.qryStop)])
                if seq == "":
                    continue
                seqFile.writelines([syri.qryChr,"\t", str(syri.qryStart), "\t",sample,"~",syri.uid,"\tN\t", seq, "\t60\tPASS\tSVTYPE=INS;SVLEN=", str(length), ";END=", str(syri.qryStop), ";SyRI=NOTAL.qry\tGT\t1/1\n"])
        if syri.type == "HDR":
            file.writelines([syri.refChr,"\t", str(syri.refStart), "\t",sample,"~",syri.uid,"\tN\t", "<HDR>", "\t60\tPASS\tSVTYPE=HDR;SVLEN=", str(length), ";END=", str(syri.refStop), ";SyRI=HDR\tGT\t1/1\n"])
    dupsDf = syriDf[(syriDf.type=="DUP") & (syriDf.dupType=="copyloss")] # dups - copyloss (ref has extra; qry has lost; -> alt = DEL)
    for index, dup in mergeIntervals(dupsDf, 'ref').iterrows():    # merge ref DUPs and iterate
        if (dup.end-dup.start) < minSyriSize:
            continue
        file.writelines([dup.chr, "\t", str(dup.start), "\t", sample,"~", dup.uid, "\tN\t", "<DUP>", "\t60\tPASS\tSVTYPE=DUP;SVLEN=", str(dup.end-dup.start), ";END=", str(dup.end), ";SyRI=DUP:copyloss\tGT\t1/1\n"])
    dupsDf = syriDf[(syriDf.type=="DUP") & (syriDf.dupType=="copygain")] # dups - copygain (ref has lost; qry has extra; -> alt = INS)
    for index, dup in mergeIntervals(dupsDf, 'ref').iterrows():    # merge ref DUPs and iterate
        if (dup.end-dup.start) < minSyriSize:
            continue
        file.writelines([dup.chr, "\t", str(dup.start), "\t", sample,"~", dup.uid, "\tN\t", "<DUP>", "\t60\tPASS\tSVTYPE=DUP;SVLEN=", str(dup.end-dup.start), ";END=", str(dup.end), ";SyRI=DUP:copygain\tGT\t1/1\n"])
