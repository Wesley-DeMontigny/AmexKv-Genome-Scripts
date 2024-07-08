library(dplyr)
library(GenomicRanges)
library(tidyverse)

blast <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/maker/final_annotation/miniprot_transcriptome_search.p90.blastn", sep="\t", header=FALSE)
gene_table <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/Supplemental Data/annotated_protein_encoding_genes.tsv", sep="\t", header=TRUE)

# Make sure you get the bounds right in the blast hits
for (i in seq(1, length(blast$V1)))
{
  s = 0
  e = 0
  if (blast$V7[i] > blast$V8[i])
  {
    s <- blast$V7[i]
    e <- blast$V8[i]
  }
  else
  {
    s <- blast$V7[i]
    e <- blast$V8[i]
  }
  blast$V7[i] <- s
  blast$V8[i] <- e
}

# Load blast hits as a range
r <- GRanges(seqnames = blast$V1, ranges = IRanges(start = blast$V7, end = blast$V8))
r$transcript_name <- blast$V2

# Get a list of everything represented in the genome. Set the CDS as a "chromosome" with its total length
unique_ids <- unique(as.factor(seqnames(r)))
cds_lengths <- c()
for (i in seq(1, length(unique_ids)))
{
  seqlengths(r)[unique_ids[i]] <- blast$V13[blast$V1 == unique_ids[i]][1]
}

### Reduce the overlaps in each CDS "chromosome" and calculate total coverage
no_overlaps <- GenomicRanges::reduce(r)
coverages <- data.frame(matrix(nrow=0, ncol=4))
colnames(coverages) <- c("ID", "Internal ID", "Cumulative Coverage", "Corresponding Transcripts")
for (i in seq(1, length(unique_ids)))
{
  cov_length <- sum(no_overlaps[no_overlaps@seqnames == unique_ids[i]]@ranges@width)
  total_length <- as.integer(seqlengths(no_overlaps)[unique_ids[i]])
  g_id = gene_table$ID[gene_table$Internal.ID == unique_ids[i]]
  # I don't know why, but if I want the CDS ID as a string, I need to paste() it
  corresponding_transcripts = paste(unique(blast$V2[blast$V1 == paste(unique_ids[i])]), collapse="; ")
  coverages[nrow(coverages) + 1,] <- c(g_id, paste(unique_ids[i]), cov_length / total_length, corresponding_transcripts)
}

# Number of CDS sequences with at least 90% coverage
length(coverages[coverages$`Cumulative Coverage` >= 0.80,]$ID)

write.csv(coverages, "C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/Supplemental Data/transcribed_stopped_orfs.csv")
