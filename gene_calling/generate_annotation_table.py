import pandas as pd 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import math

columns = [ 'ID', 'Internal ID', 'Contig', 'Range', 'Protein Length', 'UAA Stops', 'UAG Stops', 'UGA Stops', 'Orthogroup', 'A120 Ortholog', 'A25 Ortholog', 'Perkinsus Ortholog', "Symbiodinium Ortholog", "Fugacium Ortholog", "Breviolum Ortholog", "Plasmodium Ortholog", "Toxoplasma Ortholog", "Chromera Ortholog", "Vitrella Ortholog", "Theileria Ortholog", "Cryptosporidium Ortholog", 'Pfam Domains', 'Superfamily', 'PANTHER Family', 'Notes' ] 
df = pd.DataFrame(columns=columns)

### Assign new IDs and initialize the table
counter = 1
with open("./maker_miniprot_files/combined_proteins.list") as ids:
	for l in ids.readlines():
		df = df.append({ 'ID': f'AmexKv{counter}', 'Internal ID': l.rstrip().split(" ")[0], 'Contig': '', 'Range': '', 'Protein Length': 0, 'UAA Stops': 0, 'UAG Stops': 0, 'UGA Stops': 0, 'Orthogroup': '', 'A120 Ortholog': 'None', 'A25 Ortholog': 'None', 'Perkinsus Ortholog':'None', "Symbiodinium Ortholog":'None', "Fugacium Ortholog":'None', "Breviolum Ortholog":'None', "Plasmodium Ortholog":'None', "Toxoplasma Ortholog":'None', "Chromera Ortholog": 'None', "Vitrella Ortholog": 'None', "Theileria Ortholog": 'None', "Cryptosporidium Ortholog": 'None', 'Pfam Domains': 'None', 'Superfamily': 'None', 'PANTHER Family': 'None', 'Notes': ''}, ignore_index=True)
		counter += 1

### Assign contigs and ranges
with open("./maker_miniprot_files/miniprot_mRNAs.gff") as gff:
	for l in gff.readlines():
		ID = l.split("ID=")[1].split(";")[0]
		tab_delim = l.split("\t")
		df.loc[df['Internal ID'] == ID, 'Contig'] = tab_delim[0]
		df.loc[df['Internal ID'] == ID, 'Range'] = f"{tab_delim[3]}-{tab_delim[4]}"
with open("./maker_miniprot_files/maker_mRNAs.gff") as gff:
	for l in gff.readlines():
		ID = l.split("ID=")[1].split(";")[0]
		tab_delim = l.split("\t")
		df.loc[df['Internal ID'] == ID, 'Contig'] = tab_delim[0]
		df.loc[df['Internal ID'] == ID, 'Range'] = f"{tab_delim[3]}-{tab_delim[4]}"

### Get FASTA information for miniprot
with open('./maker_miniprot_files/non_fs_cds.fasta', mode='r') as fasta: 
	for record in SeqIO.parse(fasta, 'fasta'):
		name, sequence, seq_obj = record.id, str(record.seq), record.seq
		UGA = 0
		UAA = 0
		UAG = 0
		for i in range(math.floor(len(sequence)/3)):
			if sequence[i*3:i*3+3] == "TGA":
				UGA += 1
			elif sequence[i*3:i*3+3] == "TAA":
				UAA += 1
			elif sequence[i*3:i*3+3] == "TAG":
				UAG += 1
		
		df.loc[df['Internal ID'] == name, 'Protein Length'] = int(len(sequence)/3)
		df.loc[df['Internal ID'] == name, 'UGA Stops'] = UGA
		df.loc[df['Internal ID'] == name, 'UAA Stops'] = UAA
		df.loc[df['Internal ID'] == name, 'UAG Stops'] = UAG

### Get FASTA information for maker
with open('./maker_miniprot_files/run4_maker.filtered.proteins.fasta', mode='r') as fasta: 
	for record in SeqIO.parse(fasta, 'fasta'):
		name, sequence, seq_obj = record.id, str(record.seq), record.seq
		df.loc[df['Internal ID'] == name, 'Protein Length'] = len(sequence)

### Get all of the ortholog information
orthologs = pd.read_csv('./OrthoFinder/Orthologues/combined_proteins.tsv', sep='\t')
for index, row in orthologs.iterrows():
	row_orthos = row["combined_proteins"].split(", ")
	for ID in row_orthos:
		df.loc[df['Internal ID'] == ID, 'Orthogroup'] = row['Orthogroup']
		if row['Species'] == "Smic.genome.annotation.pep.longest":
			df.loc[df['Internal ID'] == ID, 'Symbiodinium Ortholog'] = row['Orthologs']
		elif row['Species'] == "symbB.v1.2.augustus.prot":
			df.loc[df['Internal ID'] == ID, 'Breviolum Ortholog'] = row['Orthologs']
		elif row['Species'] == "Symbiodinium_kawagutii.0819.final.gene":
			df.loc[df['Internal ID'] == ID, 'Fugacium Ortholog'] = row['Orthologs']
		elif row['Species'] == "A25_protein":
			df.loc[df['Internal ID'] == ID, 'A25 Ortholog'] = row['Orthologs']
		elif row['Species'] == "A120_protein":
			df.loc[df['Internal ID'] == ID, 'A120 Ortholog'] = row['Orthologs']
		elif row['Species'] == "Perkinsus_marinus_atcc_50983_gca_000006405.JCVI_PMG_1.0.pep.all":
			df.loc[df['Internal ID'] == ID, 'Perkinsus Ortholog'] = row['Orthologs']
		elif row['Species'] == "ToxoDB-9.0_TgondiiME49_AnnotatedProteins":
			df.loc[df['Internal ID'] == ID, 'Toxoplasma Ortholog'] = row['Orthologs']
		elif row['Species'] == "PiroplasmaDB-64_TequiWA_AnnotatedProteins":
			df.loc[df['Internal ID'] == ID, 'Theileria Ortholog'] = row['Orthologs']
		elif row['Species'] == "PlasmoDB-9.0_Pfalciparum3D7_AnnotatedProteins":
			df.loc[df['Internal ID'] == ID, 'Plasmodium Ortholog'] = row['Orthologs']
		elif row['Species'] == "CryptoDB-64_CveliaCCMP2878_AnnotatedProteins":
			df.loc[df['Internal ID'] == ID, 'Chromera Ortholog'] = row['Orthologs']
		elif row['Species'] == "CryptoDB-64_VbrassicaformisCCMP3155_AnnotatedProteins":
			df.loc[df['Internal ID'] == ID, 'Vitrella Ortholog'] = row['Orthologs']
		elif row['Species'] == "CryptoDB-8.0_CparvumIowaII_AnnotatedProteins":
			df.loc[df['Internal ID'] == ID, 'Cryptosporidium Ortholog'] = row['Orthologs']

### Load in InterProScan annotations
interpro = pd.read_csv('./interproscan/interproscan_miniprot.tsv', sep='\t', header=None)
for index, row in df.iterrows():
	cognate_rows = interpro[interpro[0] == row['Internal ID']]
	sf = cognate_rows[cognate_rows[3] == "SUPERFAMILY"]
	panther = cognate_rows[cognate_rows[3] == "PANTHER"]
	pfam = cognate_rows[cognate_rows[3] == "Pfam"]
	sf_str = ""
	panther_str = ""
	pfam_str = ""
	for i, r in sf.iterrows():
		sf_str += f"{r[5]} ({r[4]}, {r[6]}-{r[7]}); "
	for i, r in panther.iterrows():
		panther_str += f"{r[5]} ({r[4]}, {r[6]}-{r[7]}); "
	for i, r in pfam.iterrows():
		pfam_str += f"{r[5]} ({r[4]}, {r[6]}-{r[7]}); "
	if sf_str != "":
		row['Superfamily'] = sf_str[:-2]
	if panther_str != "":
		row['PANTHER Family'] = panther_str[:-2]
	if pfam_str != "":
		row['Pfam Domains'] = pfam_str[:-2]

df.to_csv("protein_encoding_genes.tsv", sep='\t', index=False)