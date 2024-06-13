from Bio import SeqIO, Seq
import sys

### python extract_introns.py BED_File Assembly_File Output_File [Subset]

allowed = []
if len(sys.argv) > 4:
    with open(sys.argv[4]) as subset_file:
        lines = subset_file.readlines()
        for l in lines:
            allowed.append(l.rstrip())


junction_data = []
with open(sys.argv[1]) as junction_file:
    lines = junction_file.readlines()
    if len(allowed) == 0:
        for l in lines[1:]:
            data = l.rstrip().split("\t")
            junction_data.append([data[0], int(data[1]), int(data[2]), data[3], int(data[4]), data[10], data[11], data[5]]) ### contig, start, end, junction_id, score, block_sizes, block_locations
    else:
        for l in lines[1:]:
            data = l.rstrip().split("\t")
            if data[0] in allowed:
                junction_data.append([data[0], int(data[1]), int(data[2]), data[3], int(data[4]), data[10], data[11], data[5]])


fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
with open(sys.argv[3], "w") as output_introns:
    for fasta in fasta_sequences:
        name, sequence, seq_obj = fasta.id, str(fasta.seq), fasta.seq
        for j in junction_data:
            if j[0] == name:
                if j[4] >= 75 and j[7] != "?":
                    junction_seq = sequence[j[1]:j[2]]
                    if j[7] == "-":
                        junction_seq = str(seq_obj.complement())[j[1]:j[2]]
                    else:
                        junction_seq = sequence[j[1]:j[2]]
                    block_sizes = j[5].split(",")
                    block_locations = j[6].split(",")
                    intron_seq = junction_seq[int(block_locations[0])+int(block_sizes[0])-4 : int(block_locations[1])+4]
                    output_introns.write(">" + j[3].replace("JUNC", "INTRON_BOUND") + " Score: " + str(j[4]) +"\n")
                    if j[7] == "-":
                        output_introns.write(intron_seq[::-1][0:8]+intron_seq[::-1][-8:] + "\n")
                    else:
                        output_introns.write(intron_seq[0:8]+intron_seq[-8:] + "\n")