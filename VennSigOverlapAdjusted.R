library(tidyverse)
library(ggplot2)
library(limma)

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")
d14Sampled <- read_csv("DESeqD14RNASeqSampled.csv")
KOSampled <- read_csv("DESeqKORNASeqSampled.csv") 
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")

SigGenes <- function(DESeqResults){
  QuantDat <- as.data.frame(DESeqResults) %>%
    drop_na
  SigGenes <- QuantDat[QuantDat$padj < 0.05,]
  SigList <- SigGenes$EnsemblID
  return(SigList)
}

ListCompare <- function(Quant, Total){
  QuantSigList <- SigGenes(Quant)
  TotalSigList <- SigGenes(Total)
  ids = sort(unique(c(as.character(QuantSigList),  as.character(TotalSigList) )))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% QuantSigList
    counts[i, 2] = ids[i] %in% TotalSigList
  }
  colnames(counts) = c("QuantSeq", "RNA-Seq")
  row.names(counts) = ids
  
  return(as.data.frame(counts))
}


d14UnsampledComp <- ListCompare(d14Quant, d14Unsampled)
d14SampledComp <- ListCompare(d14Quant, d14Sampled)
KOUnsampledComp <- ListCompare(KOQuant, KOUnsampled)
KOSampledComp <- ListCompare(KOQuant, KOSampled)

par(mfrow=c(2,2))
vennDiagram(d14UnsampledComp) + title("Unsampled d14 Venn")
vennDiagram(KOUnsampledComp) + title("Unsampled KO Venn")
vennDiagram(d14SampledComp) + title("Sampled d14 Venn")
vennDiagram(KOSampledComp) + title("Sampled KO Venn")
