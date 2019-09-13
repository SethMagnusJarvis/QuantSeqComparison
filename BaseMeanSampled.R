library(tidyverse)
library(DESeq2)

d14Sampled <- read_csv("DESeqD14RNASeqSampled.csv")
KOSampled <- read_csv("DESeqKORNASeqSampled.csv")
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")


PlotScatter <- function(Quant, Total)
{
  Join <- full_join(Total, Quant, by = "EnsemblID" , suffix = c(".RNASeq", ".Quant")) %>%
    na.omit
  GCPlot <- ggplot(Join, aes(x=log(baseMean.RNASeq), y=log(baseMean.Quant))) + 
    geom_point() +
    theme(text = element_text(size=20)) +
    xlab("Log Mean Reads in RNA-seq") +
    ylab("Log Mean Reads in QuantSeq")
  return(GCPlot)
}

d14Scatter <- PlotScatter(d14Sampled, d14Quant) + ggtitle("d14 RNA-seq vs QuantSeq \n mean reads using downsampled RNA-seq")
KOScatter <- PlotScatter(KOSampled, KOQuant) + ggtitle("KO RNA-seq vs QuantSeq \n mean reads using downsampled RNA-seq")

d14Scatter + KOScatter
