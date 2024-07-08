library(dplyr)
library(GenomicRanges)
library(gggenomes)
library(ggplot2)
library(tidyverse)
library(Biostrings)
library(karyoploteR)

select_mitos <- data.frame(seq_id=c("tig00002434", "tig00002369", "tig00002716", "tig00000395", "tig00000418", "tig00002685"), 
                           start=c(1, 1, 1, 1, 1, 1), 
                           end=c(41137, 56915, 29558, 18542, 31184, 14140),
                           length = c(41137, 56915, 29558, 18542, 31184, 14140)-c(1, 1, 1, 1, 1, 1))

proteins <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/protein_hit.bed", header=FALSE, sep="\t")
proteins <- proteins %>% dplyr::rename(seq_id = V1, start = V2, end = V3, feat_id = V4, strand = V6) %>%
  mutate(is_fragment = grepl("fragment", feat_id, fixed = TRUE), 
         feat_id = ifelse(feat_id == "AmexKv_CYTB_fragment" | feat_id == "AmexKv_CYTB", "cytB", 
                          ifelse(feat_id == "AmexKv_COX1_fragment" | feat_id == "AmexKv_COX1", "coxI", "coxIII")))
inverted_repeats <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/clean_inverted_repeats.bed", header=FALSE, sep="\t")
inverted_repeats <- inverted_repeats %>% dplyr::rename(seq_id = V1, start = V2, end = V3, feat_id = V4, strand = V6)

gggenomes(seqs = select_mitos, genes=proteins, feats=inverted_repeats) + geom_seq() + geom_bin_label() + 
  geom_gene(mapping = aes(fill=feat_id), stroke=0.5) + geom_feat(alpha=0.35) + 
  scale_fill_manual(values=c("#3bac58", "#02d2de", "#ebeb02")) + guides(fill=guide_legend(feat_id="Gene Product"))
