library(tidyverse)
library(DESeq2)

GetFPKM <- function (ReadList){
  coldata<-data.frame(sample = colnames(ReadList))
  coldata$Condition[1:4] <- "WT"
  coldata$Condition[5:8] <- "HOM"
  row.names(coldata) <- coldata$sample
  coldata$sample <- NULL
  coldata$Condition <- factor(coldata$Condition,levels=c("WT","HOM"))
  
  #run DeSeq on QuantSeq
  formula0 = as.formula("~ 1") 
  ReadList <- ReadList[complete.cases(ReadList),]
  CDS <- DESeqDataSetFromMatrix(countData = ReadList, colData = coldata, design = as.formula("~ Condition"))
  CDS <- DESeq(CDS, test = "LRT", reduced = formula0, minReplicatesForReplace = 5 ) 
  FPKM <- data.frame(counts(CDS, normalized=T)) %>%
    rownames_to_column("ensemblID")
  return(FPKM)
}

d14UnsampledRPKM <- read_csv("FUS_d14_rpkms.csv") %>%
  select(-ends_with("HET")) %>%
  select(ensemblID, ends_with("WT"), ends_with("HOM"))

KOUnsampledRPKM <- read_csv("FUS_KO_rpkms.csv") %>%
  select(-ends_with("HET")) %>%
  select(ensemblID, ends_with("WT"), ends_with("KO"))

d14Sampled <- read_csv("d14TotalRawSampled.csv") %>%
  select(-ends_with("HET")) %>%
  select(gene_id, ends_with("WT"), ends_with("HOM")) %>%
  column_to_rownames("gene_id")
d14SampledRPKM <- GetFPKM(d14Sampled)


KOSampled <- read_csv("KOTotalRawSampled.csv")%>%
  select(-ends_with("HET")) %>%
  select(gene_id, ends_with("WT"), ends_with("KO"))  %>%
  column_to_rownames("gene_id")
KOSampledFPKM <- GetFPKM(KOSampled)

d14QuantRaw <- as.data.frame(read.table("jack4.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
d14QuantSeq <- select(d14QuantRaw, "gene_id", "c1.d14_WT1":"t4.d14_HOM4") %>% 
  column_to_rownames("gene_id")
d14QuantFPKM <- GetFPKM(d14QuantSeq)

KOQuantRaw <- as.data.frame(read.table("jack2.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
KOQuantSeq <- select(KOQuantRaw, "gene_id", "c1.KO_WT1":"t4.KO_HOM4") %>% 
  column_to_rownames("gene_id")
KOQuantFPKM <- GetFPKM(KOQuantSeq)

write_csv(d14UnsampledRPKM, "d14UnsampledRPKM.csv")
write_csv(KOUnsampledRPKM, "KOUnsampledRPKM.csv")
write_csv(d14SampledRPKM, "d14SampledRPKM.csv")
write_csv(KOSampledFPKM, "KOSampledRPKM.csv")
write_csv(d14QuantFPKM, "d14QuantRPKM.csv")
write_csv(KOQuantFPKM, "KOQuantRPKM.csv")

