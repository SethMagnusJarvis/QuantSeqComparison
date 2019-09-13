library(tidyverse)
library(ggplot2)
library(tibble)
library(patchwork)

#Read in bed file of introns
setwd("~/Google Drive/Work/UCL/Seth/TestingQuantSeq")
bed <- read_tsv("mouse_mm10_all_introns.bed",col_names = FALSE)
colnames(bed) <- c("chrom", "start", "end", "humanname", "geneid", "Strand", "transcriptid", "frame", "type", "misc")

#Sum total intron length of gene
bed$geneid <- str_replace(bed$geneid ,"\\..*","")
bed$transcriptid <- str_replace(bed$transcriptid ,"\\..*","")
bed$length  <- bed$end - bed$start 
Length <- bed %>% group_by(geneid) %>% summarise(unique(geneid), length = sum(length), introns = n())


#Read in differential expression stuff d14
d14QuantDE <- read_csv("DESeqD14QuantSeq.csv")
d14TotalDE <-  read_csv("DESeqD14RNASeqUnsampled.csv")

KOQuantDE <- read_csv("DESeqKOQuantSeq.csv")
KOTotalDE <- read_csv("DESeqKORNASeqUnsampled.csv")

#Count PolyA signals used found
d14PolyANum <- read_tsv("d14_WT_KO_expression_sites_v2.tab") %>%
  group_by(gene_id) %>% 
  summarise(polyASignals = n())
KOPolyANum <- read_tsv("fus_WT_KO_expression_sites_v2.tab") %>%
  group_by(gene_id) %>% 
  summarise(polyASignals = n())

proportioncounting <- function(DE, SignalCount) {
  Merge <- left_join(DE, SignalCount, by=c("EnsemblID"="gene_id"))
  Merge <- Merge[complete.cases(Merge),]
  Merge$category <- cut(Merge$polyASignals, breaks=c(-Inf,1, 2, 3, 4, 5,Inf), labels=c("1","2","3", "4", "5", ">5"))
  Merged <- Merge %>%  group_by(category) %>%
    summarise(sig_genes = sum(pvalue<0.05, na.rm = TRUE),
              total_genes = n(),
              prop_sig = sig_genes / total_genes)
}

d14Quant <- proportioncounting(d14QuantDE, d14PolyANum)
d14Total <- proportioncounting(d14TotalDE, d14PolyANum)                     
KOQuant <- proportioncounting(KOQuantDE, KOPolyANum)                     
KOTotal <- proportioncounting(KOTotalDE, KOPolyANum)                     

sigdifference <- function(DE, SignalCount) {
  Merge <- left_join(DE, SignalCount, by=c("EnsemblID"="gene_id"))
  Merge <- Merge[complete.cases(Merge),]
  Merge$category <- cut(Merge$polyASignals, breaks=c(-Inf,1, 2, 3, 4, 5,Inf), labels=c("1","2","3", "4", "5", ">5"))
  Merged <- Merge %>%  group_by(category)
}

d14QuantMerge <- sigdifference(d14QuantDE, d14PolyANum)
d14TotalMerge <- sigdifference(d14TotalDE, d14PolyANum)                     
KOQuantMerge <- sigdifference(KOQuantDE, KOPolyANum)                     
KOTotalMerge <- sigdifference(KOTotalDE, KOPolyANum)   

d14QuantPlot <- ggplot(d14Quant, aes(x=category, y=prop_sig)) + labs(x="Number Of PolyA Signals", y="Proportion of Significant Genes") + ggtitle("Proportion of significant genes  by number \nof polyA signals in QuantSeq d14") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20)) + theme(text = element_text(size=20)) + theme_bw() + geom_bar(stat="identity", fill="Blue")
d14TotaltPlot <- ggplot(d14Total, aes(x=category, y=prop_sig)) + labs(x="Number Of PolyA Signals", y="Proportion of Significant Genes") + ggtitle("Proportion of significant genes  by number \nof polyA signals in RNA-seq d14") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20)) + theme(text = element_text(size=20)) + theme_bw() + geom_bar(stat="identity", fill="springgreen1")
KOQuantPlot <- ggplot(KOQuant, aes(x=category, y=prop_sig)) + labs(x="Number Of PolyA Signals", y="Proportion of Significant Genes") + ggtitle("Proportion of significant genes  by number \nof polyA signals in QuantSeq KO") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20)) + theme(text = element_text(size=20)) + theme_bw() + geom_bar(stat="identity", fill="navyblue")
KOTotaltPlot <- ggplot(KOTotal, aes(x=category, y=prop_sig)) + labs(x="Number Of PolyA Signals", y="Proportion of Significant Genes") + ggtitle("Proportion of significant genes  by number \nof polyA signals in RNA-seq KO") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20)) + theme(text = element_text(size=20))+ theme_bw() + geom_bar(stat="identity", fill="forestgreen")
d14QuantPlot+d14TotaltPlot+KOQuantPlot+KOTotaltPlot

