### Generates a CSV I can use for determining the parasite subset of the assembly
import pandas as pd

tig_info = pd.read_csv("canu_718.contigs.layout.tigInfo", delimiter="\t")
gc_info = pd.read_csv("canu_718_gc.tsv", delimiter="\t")
rna_cov_info = pd.read_csv("canu_718_merged.aln.coverage.tsv", delimiter="\t")

sequence_ids = []
at_content = []
rna_coverage = []


for index, row in tig_info.iterrows():
        missing = 8 - int(len(str(row['#tigID'])))
        seq_id = "tig" +  (missing * "0") + str(row['#tigID'])
        sequence_ids.append(seq_id)
        gc_vals = gc_info.loc[gc_info['Name'] == seq_id]["%GC"].values
        rna_vals = rna_cov_info.loc[rna_cov_info['#rname'] == seq_id]["coverage"].values
        if len(gc_vals) > 0:
                at_content.append(100-gc_vals[0])
        else:
                at_content.append(None)
        if len(rna_vals) > 0:
                rna_coverage.append(rna_vals[0])
        else:
                rna_coverage.append(None)

tig_info["at_content"] = at_content
tig_info["rna_coverage"] = rna_coverage
tig_info["sequence_id"] = sequence_ids

print(tig_info)

tig_info.to_csv("canu_718.csv", index=False)