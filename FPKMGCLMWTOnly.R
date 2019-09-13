library(tidyverse)
library(DESeq2)


GetStats <- function(RPKM){
  RPKM <- RPKM %>%
    column_to_rownames("ensemblID")
  WT <- select(RPKM, contains("WT"))
  WT$WTMean <- rowMeans(WT)
  WT <- WT %>%
    rownames_to_column("ensemblID")
  return(WT)
}

LinearMaker <- function(Quant, Total, GC)
{
  Join <- full_join(Total, Quant, by = "ensemblID" , suffix = c(".RNASeq", ".Quant")) %>%
    left_join(GC, by = "ensemblID") %>%
    mutate('Compare' = WTMean.Quant - WTMean.RNASeq)
  Join <- Join[complete.cases(Join),]
  Lm <- summary(lm(Join$Compare ~ Join$GCount))
  return(Lm)
}

d14UnsampledRPKM <- read_csv("d14UnsampledRPKM.csv")
d14UnsampledMean <- GetStats(d14UnsampledRPKM)
KOUnsampledRPKM <- read_csv("KOUnsampledRPKM.csv")
KOUnsampledMean <- GetStats(KOUnsampledRPKM)
d14SampledRPKM <- read_csv("d14SampledRPKM.csv")
d14SampledMean <- GetStats(d14SampledRPKM)
KOSampledFPKM <- read_csv("KOSampledRPKM.csv")
KOSampledMean <- GetStats(KOSampledFPKM)
d14QuantFPKM <- read_csv("d14QuantRPKM.csv")
d14QuantMean <- GetStats(d14QuantFPKM)
KOQuantFPKM <- read_csv("KOQuantRPKM.csv")
KOQuantMean<- GetStats(KOQuantFPKM)

AllGC20 <- read_csv("AllNamedGenes20WithNames.csv")

SummarisedGC <- AllGC20 %>%
  group_by(ensemblID) %>%
  dplyr::summarise(GCount = mean(GCount))

KOSampledLM <- LinearMaker(KOSampledMean, KOQuantMean, AllGC20)
d14SampledLM <- LinearMaker(d14SampledMean, d14QuantMean, AllGC20)
KOUnsampledLM <- LinearMaker(KOUnsampledMean, KOQuantMean, AllGC20)
d14UnsampledLM <- LinearMaker(d14UnsampledMean, d14QuantMean, AllGC20)

KOSampleLMMean <- LinearMaker(KOSampledMean, KOQuantMean, SummarisedGC)
d14SampledLMMean <- LinearMaker(d14SampledMean, d14QuantMean, SummarisedGC)
KOUnsampledLMMean <- LinearMaker(KOUnsampledMean, KOQuantMean, SummarisedGC)
d14UnsampledLMMean <- LinearMaker(d14UnsampledMean, d14QuantMean, SummarisedGC)


PlotScatter <- function(Quant, Total, GC)
  {
  Join <- full_join(Total, Quant, by = "ensemblID" , suffix = c(".RNASeq", ".Quant")) %>%
    left_join(GC, by = "ensemblID") %>%
    mutate('Compare' = WTMean.Quant - WTMean.RNASeq)
  Join <- Join[complete.cases(Join),]
  GCPlot <- ggplot(Join, aes(x=log(Compare), y=GCount)) + 
    geom_point() +
    theme(text = element_text(size=20)) +
    xlab("Log of difference in RPKM") +
    ylab("Number of GC")# +
    #geom_smooth(method = "lm")
  return(GCPlot)
}

KOSampleScatter <- PlotScatter(KOSampledMean, KOQuantMean, AllGC20) + ggtitle("Sampled KO and Quant RPKM vs GC WT only") 
d14SampledScatter <- PlotScatter(d14SampledMean, d14QuantMean, AllGC20) + ggtitle("Sampled d14 and Quant RPKM vs GC WT only") 
KOUnsampledScatter <- PlotScatter(KOUnsampledMean, KOQuantMean, AllGC20)  + ggtitle("Unsampled KO and Quant RPKM vs GC WT only") 
d14UnsampledScatter <- PlotScatter(d14UnsampledMean, d14QuantMean, AllGC20) + ggtitle("Unsampled d14 and Quant RPKM vs GC WT only") 

d14SampledScatter+KOSampleScatter+d14UnsampledScatter+KOUnsampledScatter
