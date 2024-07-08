#!/usr/bin/awk -f

# BEGIN block to initialize the dictionary
BEGIN {
    FS = "\t"
    while ((getline < "protein_encoding_genes.tsv") > 0) {
        dict[$2] = $1
    }
    close("protein_encoding_genes.tsv")
}

# Function to process the FASTA file
function process_fasta() {
    if (substr($0, 1, 1) == ">") {
        # Extract the key from the sequence name
        key = substr($0, 2)
        split(key, parts, " ")
        key = parts[1]
        
        # Rename the sequence name according to the dictionary
        if (key in dict) {
            print ">" dict[key]
        } else {
            print $0
        }
    } else {
        print $0
    }
}

# Main block to process each line of the FASTA file
{
    process_fasta()
}

