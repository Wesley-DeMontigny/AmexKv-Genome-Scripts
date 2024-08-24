import csv

def read_tsv(tsv_file):
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return [row for row in reader]

def parse_gff(gff_file):
    gff_entries = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            attributes = {}
            for attribute in parts[8].split(";"):
                key_value = attribute.split("=")
                if len(key_value) == 2:
                    attributes[key_value[0]] = key_value[1]
            gff_entries.append({
                'seqid': parts[0],
                'source': parts[1],
                'type': parts[2],
                'start': parts[3],
                'end': parts[4],
                'score': parts[5],
                'strand': parts[6],
                'phase': parts[7],
                'attributes': attributes
            })
    return gff_entries

def index_gff_entries(gff_data):
    indexed_gff = {}
    for entry in gff_data:
        gff_id = entry['attributes'].get('ID')
        if gff_id:
            indexed_gff[gff_id] = entry
    return indexed_gff

def write_gff(new_gff_file, tsv_data, indexed_gff):
    with open(new_gff_file, 'w') as f:
        for tsv_entry in tsv_data:
            internal_id = tsv_entry['Internal ID']
            gff_entry = indexed_gff.get(internal_id)
            if not gff_entry:
                continue

            new_attributes = {'ID': tsv_entry['ID']}
            for key, value in tsv_entry.items():
                if key not in ['ID', 'Internal ID', 'Contig', 'Range']:
                    new_attributes[key.replace(" ", "_")] = value

            attribute_str = ";".join([f"{k}={v}" for k, v in new_attributes.items()])

            f.write("\t".join([
                gff_entry['seqid'],
                "geneAnnotations",
                "gene",
                gff_entry['start'],
                gff_entry['end'],
                gff_entry['score'],
                gff_entry['strand'],
                gff_entry['phase'],
                attribute_str
            ]) + "\n")

def main(tsv_file, gff_file, output_gff):
    tsv_data = read_tsv(tsv_file)
    gff_data = parse_gff(gff_file)
    indexed_gff = index_gff_entries(gff_data)
    write_gff(output_gff, tsv_data, indexed_gff)

if __name__ == "__main__":
    tsv_file = './protein_encoding_genes.tsv'
    gff_file = './raw_maker_miniprot.gff'
    output_gff = './final_genes.gff'

    main(tsv_file, gff_file, output_gff)
