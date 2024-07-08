library(gggenomes)
library(ggplot2)
library(tidyverse)

# The indices of the stop codons do not account for the introns - here I manually adjust the indices according to intron size.
adjust_for_introns <- function(IDs, indices)
{
  corrected = c()
  for(i in seq(1, length(IDs)))
  {
    cor_i = indices[i]
    if(IDs[i] == "MP007004")
    {
      if(indices[i] > 952)
      {
        cor_i = cor_i + (1115-952)
        if(cor_i > 1289)
        {
          cor_i = cor_i + (1443-1289)
        }
      }
    }
    else if(IDs[i] == "MP008582")
    {
      if(indices[i] > 597)
      {
        cor_i = cor_i + (742-597)
        if(cor_i > 1243)
        {
          cor_i = cor_i + (1435-1243)
          if(cor_i > 1664)
          {
            cor_i = cor_i + (1837-1664)
            if(cor_i > 3201)
            {
              cor_i = cor_i + (3459-3201)
            }
          }
        }
      }
    }
    else if(IDs[i] == "MP000988")
    {
      if(indices[i] > 601)
      {
        cor_i = cor_i + (856-601)
        if(cor_i > 1311)
        {
          cor_i = cor_i + (1624-1311)
          if(cor_i > 2082)
          {
            cor_i = cor_i + (2291-2082)
            if(cor_i > 2836)
            {
              cor_i = cor_i + (3195-2836)
            }
          }
        }
      }
    }
    corrected <- append(corrected, cor_i)
  }
  return(corrected)
}

stop_indices <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/final_annotation/stop_codon_indices.tsv", sep=",", header=TRUE)


# Actin UGA End
genes <- tibble(
  seq_id = c("Actin (AmexKv886)", "Spo11 (AmexKv7515)", 
             "DNA polymerase delta, catalytic subunit (AmexKv8001)",
             "MCM7 (AmexKv5560)"),
  start = c(1, 1, 1, 1),
  end = c(5409, 1551, 4262, 4219),
  introns = list(c(114,4244,4535,4721), 
                 c(952, 1115, 1289, 1443),
                 c(597, 742, 1243, 1435, 1664, 1837, 3201, 3459),
                 c(601, 856, 1311, 1624, 2082, 2291, 2836, 3195))
)

# Only get genes I'm using
used_gene_stops <- stop_indices %>% filter(ID == "MP007004" | ID == "MP008582" | ID == "MP000988")
used_gene_stops <- used_gene_stops %>% mutate(seq_id = ifelse(ID=="MP007004", "Spo11 (AmexKv7515)", 
                                                              ifelse(ID == "MP008582", "DNA polymerase delta, catalytic subunit (AmexKv8001)", "MCM7 (AmexKv5560)")))
# Shift the index to the middle of the codon for presentation
used_gene_stops <- used_gene_stops %>% mutate(start = Index+1, end=Index+1)
# Now I need to shift the starts according to the introns
used_gene_stops <- used_gene_stops %>% mutate(start=adjust_for_introns(ID, start), end=adjust_for_introns(ID, end))
# Use variants to mark stops https://thackl.github.io/gggenomes/reference/geom_variant.html
gggenomes(genes, feat=used_gene_stops) + geom_seq() + geom_seq_label() +
  geom_gene(fill="#3bac58", color="black") + geom_variant(mapping=aes(shape=Codon), offset = 0.2)
