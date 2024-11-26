library(ggplot2)
library(dplyr)
library(gridExtra)

data <- read.csv("./best_hit_ratios.tsv", sep="\t")

data <- data %>% mutate(maker_predicted = grepl("maker", gene, fixed = TRUE))

data <- data %>% group_by(maker_predicted)

data %>% summarize(Mean = mean(best_hit_ratio))

plot1 <- ggplot(data = data, mapping = aes(best_hit_ratio)) + geom_histogram(aes(fill=maker_predicted), bins=50) + 
         scale_x_continuous(limits = c(0, 5)) + theme_classic() + ylab("Count") + 
         guides(fill="none") + xlab("Length Ratio")

plot2 <- ggplot(data = data, mapping = aes(best_hit_ratio)) + geom_histogram(aes(fill=maker_predicted), bins=50) + 
  theme_classic() + ylab("Count") + guides(fill=guide_legend(title="MAKER Predicted")) + xlab("Length Ratio")

grid.arrange(plot1, plot2, ncol=2)
