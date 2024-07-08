library(dplyr)
library(ggplot2)
library(cowplot)

assembly <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/assembly/canu_718.csv")

hits <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/clean_mito_contig_search.tblastn", sep="\t", header=FALSE)
hits <- hits %>% rename(query=V1, contig=V2, pident=V3, len=V4, mismatch=V5, 
                        gap=V6, query_start=V7, query_end=V8, contig_start=V9, 
                        contig_end=V10, evalue=V11,bitscore=V12, contig_seq=V13)


hits_strict <- hits %>% filter(pident > 95 & evalue < 1e-75)
unique_strict_contig_names <- unique(hits_strict$contig)
unique_strict_contigs <- assembly %>% filter(sequence_id %in% unique_strict_contig_names)

at_hist <- ggplot(data = unique_strict_contigs, mapping = aes(x=at_content)) + geom_histogram(bins=20, color="white", fill="#6E6B6C") +
  theme_bw() + labs(y="Contig Count", x="AT Proportion")

size_hist <- ggplot(data = unique_strict_contigs, mapping = aes(x=tigLen/1000)) + geom_histogram(bins=20, color="white", fill="#6E6B6C") +
  theme_bw() + labs(y="Contig Count", x="Contig Length (Kb)")

coverage_hist <- ggplot(data = unique_strict_contigs, mapping = aes(x=coverage)) + geom_histogram(bins=20, color="white", fill="#6E6B6C") +
  theme_bw() + labs(y="Contig Count", x="Coverage")

plot_grid(plot_grid(at_hist, size_hist, 
                    align="v", ncol=1, rel_heights = c(1, 1)), coverage_hist, rel_heights = c(1,1))

cytb_len = 342
cox3_len = 197
cox1_len = 465
hits_strict <- hits_strict %>% mutate(fragment = 
                                        ifelse(query == "AmexKv_CYTB", len/cytb_len < 0.9,
                                        ifelse(query == "AmexKv_COX3", len/cox3_len < 0.9, len/cox1_len < 0.9)))


write.csv(unique_strict_contigs, "C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/mito_contigs.csv")
write.csv(hits_strict, "C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/mito_stritct_hits.csv")
write.csv(hits, "C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/hits.csv")


bed_out <- hits_strict %>% select(c(contig, contig_start, contig_end, query, bitscore, fragment))
bed_out <- bed_out %>% mutate(query = paste(query, ifelse(fragment, "_fragment", ""), sep = ""))
bed_out <- bed_out %>% mutate(begin = ifelse(contig_start > contig_end, contig_end, contig_start), end = ifelse(contig_start > contig_end, contig_start, contig_end),
                              strand = ifelse(contig_start > contig_end, "-", "+"))

bed_out <- bed_out %>% select(contig, begin, end, query, bitscore, strand)
write.csv(bed_out, "C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/protein_hit_bed.csv")
