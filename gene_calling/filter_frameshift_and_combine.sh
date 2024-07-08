awk '$9 !~ "Frameshift" && $3 == "mRNA"' all_stop_containing_alignments.gff | awk -F'[=;]' '{print $2}' > non_fs_ids.list

awk '$9 !~ "Frameshift" && $3 == "mRNA"' all_stop_containing_alignments.gff > miniprot_mRNAs.gff

awk '$3 == "mRNA"' run4_maker_only.filtered.gff > maker_mRNAs.gff

seqtk subseq all_stop_containing_alignments.fasta non_fs_ids.list > non_fs.fasta and seqtk subseq all_stop_containing_alignments.clean.fasta non_fs_ids.list > non_fs.clean.fasta

cat run4_maker.filtered.proteins.fasta non_fs.clean.fasta > combined_proteins.fasta

awk '/^>/ {print substr($0, 2)}' combined_proteins.fasta > combined_proteins.list