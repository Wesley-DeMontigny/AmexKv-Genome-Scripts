import pysam

def is_contig_consumed(aln_file, threshold=0.9):
    bamfile = pysam.AlignmentFile(aln_file, "rb")
    consumed_contigs = []

    for contig in bamfile.header.references:
        contig_length = bamfile.header.get_reference_length(contig)
        aligned_length = 0

        for read in bamfile.fetch(contig):
            if not read.is_unmapped:
                aligned_length += read.reference_length

        if (aligned_length / contig_length) >= threshold:
            consumed_contigs.append(contig)
            print(f"{contig} is Consumed")
		else:
			print(f"{contig} is not Consumed")

    bamfile.close()
    return consumed_contigs

consumed = is_contig_consumed("new_old_alignment.sorted.bam")
print(len(consumed))