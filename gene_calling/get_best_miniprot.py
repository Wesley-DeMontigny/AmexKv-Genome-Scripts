import sys

def parse_gff(gff_file_path):
    entries = {'mRNA': [], 'CDS': []}
    with open(gff_file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9 or "Frameshift=" in line or (int(parts[4]) - int(parts[3])) < 150:
                continue

            feature_type = parts[2]
            if feature_type in ['mRNA', 'CDS']:
                entries[feature_type].append(parts)
    return entries

def find_overlapping_mRNAs(mrna_entries):
    mrna_entries.sort(key=lambda x: (x[0], int(x[3])))
    overlapping_groups = []
    current_group = []

    for entry in mrna_entries:
        if current_group and int(entry[3]) <= int(current_group[-1][4]):
            current_group.append(entry)
        else:
            if current_group:
                overlapping_groups.append(current_group)
            current_group = [entry]

    if current_group:
        overlapping_groups.append(current_group)

    return overlapping_groups

def filter_highest_score(overlapping_groups, cds_entries):
    filtered_entries = []

    for group in overlapping_groups:
        highest_score_entry = max(group, key=lambda x: float(x[5]))
        filtered_entries.append(highest_score_entry)

        for cds in cds_entries:
            if cds[8].find(f"Parent={highest_score_entry[8].split(';')[0].split('=')[1]}") != -1:
                filtered_entries.append(cds)

    return filtered_entries

def main(gff_file_path):
    entries = parse_gff(gff_file_path)
    overlapping_groups = find_overlapping_mRNAs(entries['mRNA'])
    filtered_entries = filter_highest_score(overlapping_groups, entries['CDS'])

    for entry in filtered_entries:
        print('\t'.join(entry))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python get_best_miniprot.py <GFF_file_path>")
        sys.exit(1)

    gff_file_path = sys.argv[1]
    main(gff_file_path)