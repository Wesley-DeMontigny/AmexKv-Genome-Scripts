minimap2 -ax asm5 /nas2/wesleyd/previous_amoebophrya_project/AmexKv_Old/scaffolds_amkv_spades.fasta ../genome_assembly/a_subset_10_26.fasta > new_old_alignment.sam

samtools view -Sb new_old_alignment.sam | samtools sort -o new_old_alignment.sorted.bam

python3 compare_assemblies.py