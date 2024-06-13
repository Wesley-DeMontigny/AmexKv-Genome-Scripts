minimap2 -ax splice ../genome_assembly/canu_718.contigs.fasta amexkv_trinity_dn.fasta > trinity_dn_alignment.canu718.sam

samtools view -Sb trinity_dn_alignment.canu718.sam -o trinity_dn_alignment.canu718.bam

samtools sort trinity_dn_alignment.canu718.bam -o trinity_dn_alignment.canu718.sorted.bam

samtools index trinity_dn_alignment.canu718.sorted.bam

nice -n 15 nohup ~/tools/regtools/build/regtools junctions extract -s RF -a 30 -m 40 -M 10000 trinity_dn_alignment.canu718.sorted.bam -o regtools_scan.trinity_aln.bed > regtools.trinity_aln.log &

python3 extract_introns_trinity.py regtools_scan.trinity_aln.bed ../../genome_assembly/a_subset_10_26.fasta introns.trinity.para_10_26.fasta ../../genome_assembly/a_subset_10_26.list

python3 extract_intron_boundary_regions.py regtools_scan.para_10_26.no_strand.a30.bed ../../genome_assembly/a_subset_10_26.fasta intron_bounded_region.para_subset_10_26.no_strand.a30.fasta

weblogo --format eps -A DNA -t "Intron Consensus Sequence" -c "classic" --size large --annotate 1,2,3,4,5,6,-6,-5,-4,-3,-2,-1 -S 1 < intron_bounded_region.para_subset_10_26.no_strand.a30.fasta > intron_consensus.eps