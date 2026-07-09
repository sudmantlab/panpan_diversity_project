label=$(basename $1 .vcf.gz)

truvari collapse                     \
    -i $1                            \
    -o $label~truvari_merge.vcf      \
    -c $label~truvari_collapsed.vcf  \
    --sizemin   10000                \
    --sizemax   1000000000           \
    --pctovl 0.5                     \
    -k maxqual                       \
    --chain                          \
    --pctseq 0                       \
    --refdist 5000                   \
    --pctsize .5 &> $label.LOG
