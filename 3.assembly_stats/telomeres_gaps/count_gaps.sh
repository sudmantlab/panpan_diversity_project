#!/bin/bash
#
# Count the number of gaps within a genome
# Outps a per fasta, per chromosome file and a summary file
# Requires bioawk
#

genomeDir="XXX YYY"

for dir in $genomeDir
do
    echo $dir
    echo "sample,sequence,count" > $(basename $dir)~gaps.csv
    for x in $dir/*.fasta
    do
        bioawk -c fastx -v SAMPLE=$(basename $x .fasta) '{n = split($seq, arr, /N+/); print SAMPLE","$name","n-1}' $x
    done >> $(basename $dir)~gaps.csv

    awk -v FS="," -v OFS="," 'BEGIN{print "sample,totalCount"} NR>1{counts[$1] += $3} END{for(i in counts){print i,counts[i]}}'  $(basename $dir)~gaps.csv > $(basename $dir)~summary.csv
    awk -v FS="," -v OFS="," 'BEGIN{print "sample,totalCount"} NR>1{if($2!~/X/ || $2!~/Y/){counts[$1] += $3}} END{for(i in counts){print i,counts[i]}}'  $(basename $dir)~gaps.csv > $(basename $dir)~summary-autosomes.csv
done
