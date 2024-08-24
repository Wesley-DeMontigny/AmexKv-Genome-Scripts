#!/bin/bash

cat <(awk '{if(/^##FASTA/){exit} else {print $0}}' ./maker_miniprot_files/para_10_26.run4.miniprot.all.gff) ./maker_miniprot_files/all_stop_containing_alignments.gff > raw_maker_miniprot.gff

echo "Adding final gene entries"
python add_final_entries.py

cat final_genes.gff raw_maker_miniprot.gff > final.gff
