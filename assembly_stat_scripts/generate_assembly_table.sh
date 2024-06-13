infoseq canu_718.contigs.fasta -only -desc -name -pgc > canu_718_gc.out

minimap2 -ax sr canu_718.contigs.fasta ../transcriptome_assembly/SRR1610332_1.fastq ../transcriptome_assembly/SRR1610332_2.fastq > canu_718_SRR1610332.aln.sam

minimap2 -ax sr canu_718.contigs.fasta ../transcriptome_assembly/SRR1610333_1.fastq ../transcriptome_assembly/SRR1610333_2.fastq > canu_718_SRR1610333.aln.sam

samtools merge -o canu_718_merged.aln.sorted.bam canu_718_SRR1610333.aln.sorted.bam canu_718_SRR1610332.aln.sorted.bam

samtools coverage canu_718_merged.aln.sorted.bam -o canu_718_merged.aln.coverage

awk '{print $1 "\t" $6}' canu_718_merged.aln.coverage > canu_718_merged.aln.coverage.tsv

python3 generate_tig_csv.py