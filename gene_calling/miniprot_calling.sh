# Parameters
GENOME_FILE="a_subset_10_26.genome"
GENOME_FASTA="../../../genome_assembly/a_subset_10_26.fasta"
BLAST_DB="/data1/users/wdem/canu_718/annotation/maker/Amoebophrya_Species_Peptides.fasta"
MAKER_PREFIX="para_10_26.FINAL_run3"
MAX_JOBS=10

# Function to perform BLAST
do_blast() {
    # I DON'T KNOW WHY - IT ONLY LIKES WHEN I PUT THE BLAST DB IN MANUALLY
    fasta_file=$1
    output_file="$fasta_file.blastout"
    blastx -query "$fasta_file" -subject /data1/users/wdem/canu_718/annotation/maker/Amoebophrya_Species_Peptides.fasta -outfmt "6 std" -evalue 1e-5 > $output_file
    echo "$fasta_file Searched."
}
# Export function so it's available to parallel
export -f do_blast


# Get Gene Ranges
awk '$3 == "mRNA"' $MAKER_PREFIX.all.gff > $MAKER_PREFIX.mRNA.gff
sort -k1,1 -k4,4n $MAKER_PREFIX.mRNA.gff > $MAKER_PREFIX.mRNA.sorted.gff

# Get Intergenic Ranges
bedtools complement -i $MAKER_PREFIX.mRNA.sorted.gff -g $GENOME_FILE > maker_intergenic.bed
bedtools getfasta -fi $GENOME_FASTA -bed maker_intergenic.bed > maker_intergenic.fasta

#Split Intergenic Ranges
~/tools/seqkit split -i maker_intergenic.fasta
pushd maker_intergenic.fasta.split

# Find all FASTA files and run BLAST in parallel
find ./ -name "*.fasta" | parallel -j $MAX_JOBS do_blast

#Concatenate Output
cat *.blastout > parallel_blast_intergenic.out
popd

# Extend the blast hits
python3 process_ranges.py

# Get FASTA of extended blast hits
bedtools getfasta -fi $GENOME_FASTA -bed gene_ranges.bed > gene_ranges.fasta

# Run MiniProt
miniprot -t16 -d gene_ranges.mpi gene_ranges.fasta
miniprot -Iut16 --trans --gff -P MP2 gene_ranges.mpi $BLAST_DB > intergenic_miniprot.out