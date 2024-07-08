awk '/^>/ {OUT=substr($1,2) ".fa";print " ">OUT}; OUT{print >OUT}' ../clean_mito_contigs.fasta


echo "" > all.50.repeats
for filename in ./*.fa; do
    echo "SOURCE:$filename" >> all.50.repeats
    repeat-match -n 50 $filename >> all.50.repeats
done