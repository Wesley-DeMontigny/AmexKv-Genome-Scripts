import pyranges as pr
import pandas as pd

# Function to process each line of the BLAST output
def process_blast_line(line):
    # Splitting the line into columns
    cols = line.split('\t')
    try:
        # Extracting the sequence name and hit indices

        sequence_name, _, _, _, _, _, hit_start, hit_end, *_ = cols
        contig, indices = sequence_name.split(':')
        start_index, end_index = map(int, indices.split('-'))

        # Converting hit indices to int
        hit_start, hit_end = int(hit_start), int(hit_end)
        if(hit_start > hit_end):
            temp = hit_start
            hit_start = hit_end
            hit_end = temp

        # Translating hit indices to larger contig coordinates
        contig_hit_start = start_index + hit_start - 1
        contig_hit_end = start_index + hit_end - 1

        # Extending indices by 2500 bp in each direction, ensuring not to exceed subsequence boundaries
        extended_start = max(start_index, contig_hit_start - 2500)
        extended_end = min(end_index, contig_hit_end + 2500)

        return [contig, extended_start, extended_end]
    except:
        pass

with open("./maker_intergenic.fasta.split/parallel_blast_intergenic.out", "r") as file:
    df = pd.DataFrame([], columns=["Chromosome", "Start", "End"])

    for line in file:
        data = process_blast_line(line)
        if(data):
            df = pd.concat([df, pd.DataFrame([data], columns=df.columns)], ignore_index=True)

    gr = pr.PyRanges(df)

    clusters = gr.cluster()
    reduced = clusters.merge()


    print(reduced)
    print("Exporting to BED (gene_ranges.bed)...")
    reduced.to_bed("gene_ranges.bed")