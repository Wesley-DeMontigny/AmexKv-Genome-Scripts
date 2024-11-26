#!/bin/bash

MAX_JOBS=10

do_search() {
    orthogroup=$1
    echo "Debug: Received Orthogroup = '$orthogroup'"
    orthogroupfile="./OrthoFinder/Orthogroup_Sequences/$orthogroup.fa"
    blastp -subject $orthogroupfile -query $orthogroupfile -outfmt '6 qseqid sseqid qlen slen bitscore' > ./OrthoBlast/$orthogroup.blastp
}

export -f do_search

awk '$9 != "None" && $9 != "Stops" {print $9}' protein_encoding_genes.tsv | sort -u | parallel -j $MAX_JOBS do_search