#!/bin/bash

awk '
BEGIN {{ OFS = "\t" }}

$2 == "miniprot" {{
    split($1, array, /[:-]/);
    contigName = array[1];
    subseqStart = array[2];

    # Adjusting start and end coordinates
    $1 = contigName;
    $4 = $4 + subseqStart;
    $5 = $5 + subseqStart;

    # Reconstruct the line without columns 10 and 11
    print $1, $2, $3, $4, $5, $6, $7, $8, $9;
}}
' intergenic_miniprot.out > intergenic_miniprot.coord_adjusted.gff

python3 get_best_miniprot.py intergenic_miniprot.coord_adjusted.gff > intergenic_miniprot.best.coord_adjusted.gff

awk '$3 == "mRNA" {split($9, array, /[=;]/); print(array[2])}' intergenic_miniprot.best.coord_adjusted.gff > intergenic_miniprot.best.list

awk '/^##STA/{sequence=$2; getline; if ($1 ~ /^tig/){split($9, idArray, ";"); gsub("ID=", "", idArray[1]); print ">" idArray[1] "\n" sequence}}' intergenic_miniprot.out > intergenic_miniprot.raw_peptide.fasta

seqtk subseq intergenic_miniprot.raw_peptide.fasta intergenic_miniprot.best.list > intergenic_miniprot.best_peptide.fasta

gff3_to_fasta -g intergenic_miniprot.best.coord_adjusted.gff -f ../../../genome_assembly/a_subset_10_26.fasta -u mRNA CDS -st user_defined -d complex -o intergenic_miniprot.best
