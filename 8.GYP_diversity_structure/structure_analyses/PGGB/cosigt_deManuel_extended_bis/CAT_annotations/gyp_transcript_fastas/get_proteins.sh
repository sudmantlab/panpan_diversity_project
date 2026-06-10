cat gyp_cds_chimp.fa | sed "s/ transcript=/@/g" >gyp_cds_chimp.formatted.fa
cat gyp_cds_human.fa | sed "s/ transcript=/@/g" >gyp_cds_human.formatted.fa
transeq -sequence gyp_cds_chimp.formatted.fa -outseq gyp_cds_chimp.formatted.translated.fa
transeq -sequence gyp_cds_human.formatted.fa -outseq gyp_cds_human.formatted.translated.fa
sed -i ""  "s/*//g"  gyp_cds_chimp.formatted.translated.fa
sed -i ""  "s/*//g"  gyp_cds_human.formatted.translated.fa

java -cp /Users/petersudmant/Documents/science/programs/homologhmm_1.05/biojava.jar:/Users/petersudmant/Documents/science/programs/homologhmm_1.05/homologhmm.jar:/Users/petersudmant/Documents/science/programs/homologhmm_1.05 se.ki.cgb.hmmdecode.Phobius gyp_cds_chimp.formatted.translated.fa  >gyp_cds_chimp.formatted.translated.phobius.output.txt
echo -e "ID\tFT\tType\tStart\tEnd\tComment" > gyp_cds_chimp.formatted.translated.phobius.output.tidy.tsv && awk '/^ID/ {id=$2} /^FT/ {c=""; for(i=5;i<=NF;i++) c=(c==""?$i:c" "$i); printf "%s\t%s\t%s\t%s\t%s\t%s\n", id, $1, $2, $3, $4, (c==""?"NA":c)}' gyp_cds_chimp.formatted.translated.phobius.output.txt >> gyp_cds_chimp.formatted.translated.phobius.output.tidy.tsv


java -cp /Users/petersudmant/Documents/science/programs/homologhmm_1.05/biojava.jar:/Users/petersudmant/Documents/science/programs/homologhmm_1.05/homologhmm.jar:/Users/petersudmant/Documents/science/programs/homologhmm_1.05 se.ki.cgb.hmmdecode.Phobius gyp_cds_human.formatted.translated.fa  >gyp_cds_human.formatted.translated.phobius.output.txt
echo -e "ID\tFT\tType\tStart\tEnd\tComment" > gyp_cds_human.formatted.translated.phobius.output.tidy.tsv && awk '/^ID/ {id=$2} /^FT/ {c=""; for(i=5;i<=NF;i++) c=(c==""?$i:c" "$i); printf "%s\t%s\t%s\t%s\t%s\t%s\n", id, $1, $2, $3, $4, (c==""?"NA":c)}' gyp_cds_human.formatted.translated.phobius.output.txt >> gyp_cds_human.formatted.translated.phobius.output.tidy.tsv
