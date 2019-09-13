library(tidyverse)
library(DESeq2)
library(patchwork)

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")
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
    ylab("Log Mean Reads in QuantSeq") +
    expand_limits(x = 0, y = 0)
  return(GCPlot)
}

d14Scatter <- PlotScatter(d14Unsampled, d14Quant) + ggtitle("d14 RNA-seq vs QuantSeq mean reads")
KOScatter <- PlotScatter(KOUnsampled, KOQuant) + ggtitle("KO RNA-seq vs QuantSeq mean reads")

d14Scatter + KOScatter
