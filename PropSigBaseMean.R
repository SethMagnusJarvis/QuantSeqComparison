library(tidyverse)
library(ggplot2)
library(patchwork)

CountGraph <- function(BaseMeans){
  BaseMeans$category <- cut(BaseMeans$baseMean, breaks=c(-Inf,10, 20, 40, 100, 1000,Inf), 
                                    labels=c("<10","10-20", "20-39", "40-99", "100-999", ">1000"))
  Counted<- BaseMeans %>% 
    group_by(category) %>%
    summarise(sig_genes = sum(pvalue<0.05, na.rm = TRUE),
              total_genes = n(),
              prop_sig = sig_genes / total_genes)
  
  PropGraph <- ggplot(Counted, aes(x=category, y=prop_sig)) + 
    geom_bar(stat="identity") + labs(x="Mean Reads of genes") + 
    theme(text = element_text(size=20)) +
    ylab("Proportoin of significant genes")
  return(PropGraph)
}

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")

d14QuantProp <- CountGraph(d14Quant) + ggtitle("QuantSeq d14 proportoin of significant reads \nseperated by the number of reads") + theme_bw() + geom_bar(stat="identity", fill="#99CCFF")
d14UnsampledProp <- CountGraph(d14Unsampled) + ggtitle("RNA-seq d14 proportoin of genes that are significant \nseperated by the number of reads") + theme_bw() + geom_bar(stat="identity", fill='#99FF99')
KOQuantProp <- CountGraph(KOQuant) + ggtitle("QuantSeq KO proportoin of significant reads \nseperated by the number of reads") + theme_bw() + geom_bar(stat="identity", fill="#0080FF")
KOUnsampledProp <- CountGraph(KOUnsampled) + ggtitle("RNA-seq KO proportoin of significant reads \nseperated by the number of reads") + theme_bw() + geom_bar(stat="identity", fill="#33FF33")

d14QuantProp + d14UnsampledProp + KOQuantProp + KOUnsampledProp
