library(tidyverse)
library(ggplot2)

LinearMaker <- function(RNA, Quant, GC){
  Join <- full_join(RNA, Quant, by = "EnsemblID" , suffix = c(".RNASeq", ".Quant")) %>%
    left_join(GC, by = c("EnsemblID" = "ensemblID")) %>%
    mutate('Diff' = log2FoldChange.Quant - log2FoldChange.RNASeq)
  Lm <- summary(lm(Join$Diff ~ Join$GCount))
}

LinearMakerSplit <- function(Seq, GC){
  Join <- left_join(Seq, GC, by = c("EnsemblID" = "ensemblID"))
  Lm <- summary(lm(Join$log2FoldChange ~ Join$GCount))
}

PlotScatter <- function(Quant, Total, GC)
{
  Join <- full_join(Total, Quant, by = "EnsemblID" , suffix = c(".RNASeq", ".Quant")) %>%
    left_join(GC, by = c("EnsemblID" = "ensemblID")) %>%
    mutate('Compare' = log2FoldChange.Quant - log2FoldChange.RNASeq)
  Join <- Join[complete.cases(Join),]
  GCPlot <- ggplot(Join, aes(x=log(Compare), y=GCount)) + 
    geom_point() +
    theme(text = element_text(size=20)) +
    xlab("Fold change difference") +
    ylab("Number of GC")# +
  #geom_smooth(method = "lm")
  return(GCPlot)
}

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")
d14Sampled <- read_csv("DESeqD14RNASeqSampled.csv")
KOSampled <- read_csv("DESeqKORNASeqSampled.csv") 
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")

AllGC20 <- read_csv("AllNamedGenes20WithNames.csv")

SummarisedGC <- AllGC20 %>%
  group_by(ensemblID) %>%
  dplyr::summarise(GCount = mean(GCount))

KOSampleLM <- LinearMaker(KOSampled, KOQuant, AllGC20)
d14SampledLM <- LinearMaker(d14Sampled, d14Quant, AllGC20)
KOUnsampledLM <- LinearMaker(KOUnsampled, KOQuant, AllGC20)
d14UnsampledLM <- LinearMaker(d14Unsampled, d14Quant, AllGC20)

KOSampleLMMean <- LinearMaker(KOSampled, KOQuant, SummarisedGC)
d14SampledLMMean <- LinearMaker(d14Sampled, d14Quant, SummarisedGC)
KOUnsampledLMMean <- LinearMaker(KOUnsampled, KOQuant, SummarisedGC)
d14UnsampledLMMean <- LinearMaker(d14Unsampled, d14Quant, SummarisedGC)


KOSampleOnlyLM <- LinearMakerSplit(KOSampled, AllGC20)
d14SampledOnlyLM <- LinearMakerSplit(d14Sampled, AllGC20)
KOUnsampledOnlyLM <- LinearMakerSplit(KOUnsampled, AllGC20)
d14UnsampledOnlyLM <- LinearMakerSplit(d14Unsampled, AllGC20)
KOQuantOnlyLM <- LinearMakerSplit(KOQuant, AllGC20)
d14UQuantOnlyLM <- LinearMakerSplit(d14Quant, AllGC20)

KOSampleScatter <- PlotScatter(KOSampled, KOQuant, AllGC20) + ggtitle("Sampled KO and Quant log fold change vs GC") 
d14SampledScatter <- PlotScatter(d14Sampled, d14Quant, AllGC20) + ggtitle("Sampled d14 and Quant log fold change vs GC") 
KOUnsampledScatter <- PlotScatter(KOUnsampled, KOQuant, AllGC20)  + ggtitle("Unsampled KO and Quant log fold change vs GC") 
d14UnsampledScatter <- PlotScatter(d14Unsampled, d14Quant, AllGC20) + ggtitle("Unsampled d14 and Quant log fold change vs GC") 

d14SampledScatter+KOSampleScatter+d14UnsampledScatter+KOUnsampledScatter
