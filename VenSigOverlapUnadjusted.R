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

SigGenesUnadjusted <- function(DESeqResults){
  QuantDat <- as.data.frame(DESeqResults) %>%
    drop_na
  SigGenes <- QuantDat[QuantDat$pvalue < 0.05,]
  SigList <- SigGenes$EnsemblID
  return(SigList)
}

ListCompareTotalUnadjusted <- function(Quant, Total){
  QuantSigList <- SigGenes(Quant)
  TotalSigList <- SigGenesUnadjusted(Total)
  ids = sort(unique(c(as.character(QuantSigList),  as.character(TotalSigList) )))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% QuantSigList
    counts[i, 2] = ids[i] %in% TotalSigList
  }
  colnames(counts) = c("QuantSeq", "RNA-seq")
  row.names(counts) = ids
  
  return(as.data.frame(counts))
}

ListCompareQuantUnsampled <- function(Quant, Total){
  QuantSigList <- SigGenesUnadjusted(Quant)
  TotalSigList <- SigGenes(Total)
  ids = sort(unique(c(as.character(QuantSigList),  as.character(TotalSigList) )))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% QuantSigList
    counts[i, 2] = ids[i] %in% TotalSigList
  }
  colnames(counts) = c("QuantSeq", "RNA-seq")
  row.names(counts) = ids
  
  return(as.data.frame(counts))
}

d14UnsampledTotalUn <- ListCompareTotalUnadjusted(d14Quant, d14Unsampled)
d14SampledTotalUn <- ListCompareTotalUnadjusted(d14Quant, d14Sampled)

d14UnsampledQuantUn <- ListCompareQuantUnsampled(d14Quant, d14Unsampled)
d14SampledQuantUn <- ListCompareQuantUnsampled(d14Quant, d14Sampled)

KOUnsampledTotalUn <- ListCompareTotalUnadjusted(KOQuant, KOUnsampled)
KOSampledTotalUn <-  ListCompareTotalUnadjusted(KOQuant, KOSampled)

KOUnsampledQuantUn <- ListCompareQuantUnsampled(KOQuant, KOUnsampled)
KOSampledQuantUn <-  ListCompareQuantUnsampled(KOQuant, KOSampled)

par(mfrow=c(2,2))
vennDiagram(d14UnsampledTotalUn) + title("Unsampled d14 with unadjusted RNA-seq")
vennDiagram(KOUnsampledTotalUn) + title("Unsampled KO with unadjusted RNA-seq")
vennDiagram(d14UnsampledQuantUn) + title("Unsampled d14 with unadjusted QuantSeq")
vennDiagram(KOUnsampledQuantUn) + title("Unsampled KO with unadjusted QuantSeq")

vennDiagram(d14SampledTotalUn) + title("Sampled d14 with unadjusted RNA-seq")
vennDiagram(KOSampledTotalUn) + title("Sampled KO with unadjusted RNA-seq")
vennDiagram(d14SampledQuantUn) + title("Sampled d14 with unadjusted QuantSeq")
vennDiagram(KOSampledQuantUn) + title("Sampled KO with unadjusted QuantSeq")
