library(tidyverse)
library(DESeq2)

RunDEseq <- function (ReadList){
  coldata<-data.frame(sample = colnames(ReadList))
  coldata$Condition[1:4] <- "WT"
  coldata$Condition[5:8] <- "HOM"
  row.names(coldata) <- coldata$sample
  coldata$sample <- NULL
  coldata$Condition <- factor(coldata$Condition,levels=c("WT","HOM"))
  
  #run DeSeq on QuantSeq
  formula0 = as.formula("~ 1") 
  CDS <- DESeqDataSetFromMatrix(countData = ReadList, colData = coldata, design = as.formula("~ Condition"))
  CDS <- DESeq(CDS, test = "LRT", reduced = formula0, minReplicatesForReplace = 5 ) 
  deseq.res <- as.data.frame(results(CDS)) %>% 
    rownames_to_column("EnsemblID")
  return(deseq.res)
}

d14Sampled <- read_csv("d14TotalRawSampled.csv") %>%
  select(-ends_with("HET")) %>%
  select(gene_id, ends_with("WT"), ends_with("HOM")) %>%
  column_to_rownames("gene_id")

KOSampled <- read_csv("KOTotalRawSampled.csv")%>%
  select(-ends_with("HET")) %>%
  select(gene_id, ends_with("WT"), ends_with("KO"))  %>%
  column_to_rownames("gene_id")

d14Sampled[is.na(d14Sampled)] <- 0
KOSampled[is.na(KOSampled)] <- 0

d14QuantRaw <- as.data.frame(read.table("jack4.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
d14QuantSeq <- select(d14QuantRaw, "gene_id", "c1.d14_WT1":"t4.d14_HOM4") %>% 
  column_to_rownames("gene_id")

KOQuantRaw <- as.data.frame(read.table("jack2.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
KOQuantSeq <- select(KOQuantRaw, "gene_id", "c1.KO_WT1":"t4.KO_HOM4") %>% 
  column_to_rownames("gene_id")

d14TotalSEQSampled <- RunDEseq(d14Sampled)
KOTotalSEQSampled <- RunDEseq(KOSampled)

d14QuantSEQ <- RunDEseq(d14QuantSeq)
KOQuantSEQ <- RunDEseq(KOQuantSeq)

d14TotalSEQUnsampled <- as.data.frame(read.table("deseq_Nicol_d14_differential_expression.tab", sep="\t",header=TRUE))[,1:7]
KOTotalSEQUnsampled <- as.data.frame(read.table("deseq_Nicol_KO_differential_expression.tab", sep="\t",header=TRUE))[,1:7]

write_csv(d14TotalSEQSampled, "DESeqD14RNASeqSampled.csv")
write_csv(KOTotalSEQSampled, "DESeqKORNASeqSampled.csv")
write_csv(d14QuantSEQ, "DESeqD14QuantSeq.csv")
write_csv(KOQuantSEQ, "DESeqKOQuantSeq.csv")
write_csv(d14TotalSEQUnsampled, "DESeqD14RNASeqUnsampled.csv")
write_csv(KOTotalSEQUnsampled, "DESeqKORNASeqUnsampled.csv")
