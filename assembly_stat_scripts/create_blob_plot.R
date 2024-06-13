library(ggplot2)
library(dplyr)
library(scales)
library(ggpointdensity)
library(cowplot)

data <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/assembly/canu_718.csv")
mitos <- read.csv("C:/Users/wescd/Desktop/Current Projects/IMET 2023/Canu718/mitochondria/Clean_Run/mito_contigs.csv")

contigs <- data %>% filter(coverage >= 1)

contigs <- contigs %>% filter(at_content > 25) %>% filter(at_content < 70) %>% filter(coverage <= 1e3)
contigs <- contigs %>% mutate(tig_owner= ifelse(at_content > 63 & at_content < 67 & coverage >=20 & rna_coverage > 50, "Amoebophrya", ifelse(sequence_id %in% mitos$sequence_id,"Mito","Other")))

a_subset <- contigs %>% filter(tig_owner == "Amoebophrya")

mainPlot <- ggplot(contigs, mapping=aes(at_content, coverage)) + geom_point(aes(size=tigLen/1000, color=tig_owner), alpha=0.3)
mainPlot <- mainPlot + scale_y_continuous(trans='log10', breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(10^.x))) +
  theme_bw() + guides(shape = guide_legend(override.aes = list(size = 1))) +
  guides(color = "none") + scale_color_manual(values = c("#3bac58", "#eec069", "#de425b")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.box = "horizontal", legend.margin=margin(c(0,0,0,0)),
        axis.text = element_text(size = 8)) +
  geom_rug(mapping=aes(color=tig_owner), alpha=0.05) + labs(size="Contig Length (kb)", y="Coverage",
                                               x="AT Content")

contigs_kb <- contigs %>% mutate(tigLen = tigLen/1000000)

topPlot <- ggplot(contigs_kb) + geom_histogram(aes(x=at_content, weight=tigLen), bins=100, color="white", fill="#6E6B6C") + theme_bw() + 
  theme(legend.position="none", axis.text = element_text(size = 8)) + xlab(NULL) + scale_x_continuous(labels = NULL, 
                                                                                                      breaks=c(0.3, 0.4, 0.5, 0.6, 0.7)) +
  labs(y="Span (Mb)")

sidePlot <- ggplot(contigs_kb) + geom_histogram(aes(x=coverage, weight=tigLen), bins=100, color="white", fill="#6E6B6C") + 
  coord_flip() + theme_bw() + theme(legend.position="none", axis.text = element_text(size = 8), font="arial") + 
  xlab(NULL) + scale_x_continuous(trans='log10', breaks=trans_breaks('log10', function(x) 10^x), labels=NULL) +
  labs(y="Span (Mb)")

plot_grid(plot_grid(topPlot, mainPlot + theme(legend.position="none"), 
                    align="v", ncol=1, rel_heights = c(20, 80)), 
          plot_grid(get_legend(mainPlot), sidePlot, rel_heights =c (20,80),
                    align="", ncol=1), rel_widths = c(80, 20))
