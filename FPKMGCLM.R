library(tidyverse)
library(DESeq2)

GetStats <- function(RPKM){
  RPKM <- RPKM %>%
    column_to_rownames("ensemblID")
  WT <- select(RPKM, contains("WT"))
  WT$WTMean <- rowMeans(WT)
  WT <- WT %>%
    rownames_to_column("ensemblID")
  HOM <- select(RPKM, -contains("WT"))
  HOM$HOMMean <- rowMeans(HOM)
  HOM <- HOM %>%
    rownames_to_column("ensemblID")
  Means <- full_join(WT, HOM, by="ensemblID")
  Means <- mutate(Means, "Diff"= WTMean - HOMMean)
  return(Means)
}

LinearMaker <- function(Quant, Total, GC)
{
  Join <- full_join(Total, Quant, by = "ensemblID" , suffix = c(".RNASeq", ".Quant")) %>%
    left_join(GC, by = "ensemblID") %>%
    mutate('Compare' = Diff.Quant - Diff.RNASeq)
  Join <- Join[complete.cases(Join),]
  Lm <- summary(lm(Join$Compare ~ Join$GCount))
}

LMBreakup <- LinearMaker <- function(RPKM, GC){
  Join <- left_join(RPKM, GC, by = "ensemblID") 
  Join <- Join[complete.cases(Join),]
  Lm <- summary(lm(Join$Diff ~ Join$GCount))
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


KOSampleLM <- LinearMaker(KOSampledMean, KOQuantMean, AllGC20)
d14SampleLM <- LinearMaker(d14SampledMean, d14QuantMean, AllGC20)
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
    mutate('Compare' = Diff.Quant - Diff.RNASeq)
  Join <- Join[complete.cases(Join),]
  GCPlot <- ggplot(Join, aes(x=log(Compare), y=GCount)) + 
    geom_point() +
    theme(text = element_text(size=20)) +
    xlab("Log of difference in RPKM") +
    ylab("Number of GC")
  return(GCPlot)
}

KOSampleScatter <- PlotScatter(KOSampledMean, KOQuantMean, AllGC20) + ggtitle("Sampled KO and Quant RPKM vs GC") 
d14SampledScatter <- PlotScatter(d14SampledMean, d14QuantMean, AllGC20) + ggtitle("Sampled d14 and Quant RPKM vs GC") 
KOUnsampledScatter <- PlotScatter(KOUnsampledMean, KOQuantMean, AllGC20)  + ggtitle("Unsampled KO and Quant RPKM vs GC") 
d14UnsampledScatter <- PlotScatter(d14UnsampledMean, d14QuantMean, AllGC20) + ggtitle("Unsampled d14 and Quant RPKM vs GC") 

d14SampledScatter+KOSampleScatter+d14UnsampledScatter+KOUnsampledScatter

KOSampleOnlyLM <- LMBreakup(KOSampledMean, AllGC20)
d14SampledOnlyLM <- LMBreakup(d14SampledMean, AllGC20)
KOUnsampledOnlyLM <- LMBreakup(KOUnsampledMean, AllGC20)
d14UnsampledOnlyLM <- LMBreakup(d14UnsampledMean, AllGC20)
KOQuantOnlyLM <- LMBreakup(KOQuantMean, AllGC20)
d14UQuantOnlyLM <- LMBreakup(d14QuantMean, AllGC20)
