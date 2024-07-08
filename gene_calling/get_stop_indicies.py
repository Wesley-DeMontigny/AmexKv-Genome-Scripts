import pandas as pd 
from Bio import SeqIO
import math

df = pd.DataFrame(columns=["ID", "Index", "Codon"])

with open('./maker_miniprot_files/non_fs_cds.fasta', mode='r') as fasta: 
	for record in SeqIO.parse(fasta, 'fasta'):
		name, sequence, seq_obj = record.id, str(record.seq), record.seq
		UGA = 0
		UAA = 0
		UAG = 0
		for i in range(math.floor(len(sequence)/3)):
			if sequence[i*3:i*3+3] in ["TGA", "TAA", "TAG"]:
				df = df.append({"ID": name, "Index": (i*3)+1, "Codon": sequence[i*3:i*3+3]}, ignore_index=True)
df.to_csv("stop_codon_indicies.tsv", sep=',', index=False)