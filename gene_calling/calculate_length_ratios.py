import pandas as pd
import os

blast_files = [f for f in os.listdir(".") if f.endswith(".blastp")]

ratios = pd.DataFrame(columns = ["og", "gene", "best_hit_ratio"])

for f in blast_files:
	df = pd.read_csv(f, sep="\t", header=None)
	filtered_df = df[df[0].apply(lambda x: "maker" in str(x) or "MP" in str(x)) & df[1].apply(lambda x: "maker" not in str(x) and "MP" not in str(x))]
	genes = filtered_df.groupby(filtered_df.columns[0])

	for gene_name, gene_data in genes:
		best_hit = gene_data[4].idxmax()
		ratio = df.loc[best_hit, df.columns[2]] / df.loc[best_hit, df.columns[3]]
		ratios.loc[len(ratios)] = {"og":f.split(".")[0], "gene": df.loc[best_hit, df.columns[0]], "best_hit_ratio": ratio}

ratios.to_csv("best_hit_ratios.tsv", sep="\t", index=False)