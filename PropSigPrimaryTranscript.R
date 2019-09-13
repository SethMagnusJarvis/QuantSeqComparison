library(tidyverse)
library(ggplot2)

setwd("~/Google Drive/Work/UCL/Seth/TestingQuantSeq")
d14QuantDE <- read_csv("DESeqD14QuantSeq.csv")
d14TotalDE <-  read_csv("DESeqD14RNASeqUnsampled.csv")
d14Transcripts <- read_csv("d14Transcripts.csv")[,-1]

KOQuantDE <- read_csv("DESeqKOQuantSeq.csv")
KOTotalDE <- read_csv("DESeqKORNASeqUnsampled.csv")
KOTranscripts <- read_csv("KOTranscripts.csv")[,-1]

Length <- read_csv('GencodeLength.csv')

d14TranscriptLength <- left_join(Length, d14Transcripts[,1:3], by=c("transcriptid"="ensembl_transcript_id")) %>%
  filter(transcript_appris =="principal1")

KOTranscriptLength <- left_join(Length, KOTranscripts[,1:3], by=c("transcriptid"="ensembl_transcript_id")) %>%
  filter(transcript_appris =="principal1")

BinPrincipal1 <- function(DE, Transcripts){
  Table <- left_join(DE, Transcripts, by=c("EnsemblID"="ensembl_gene_id"))
  Table$category <- cut(Table$length, breaks=c(-Inf,100, 1000, 10000,Inf), labels=c("<100", "100-1000", "1000-10000", ">10000"))
  TableFilt <- filter(Table, transcript_appris =="principal1")
  return(TableFilt)
}

d14QuantBreak <- BinPrincipal1(d14QuantDE, d14TranscriptLength)
d14TotalBreak <- BinPrincipal1(d14TotalDE, d14TranscriptLength)

KOQuantBreak <- BinPrincipal1(KOQuantDE, KOTranscriptLength)
KOTotalBreak <- BinPrincipal1(KOTotalDE, KOTranscriptLength)

SummarisePval <- function(Break){
  Count <- Break %>% 
    group_by(category) %>%
    summarise(sig_genes = sum(pvalue<0.05, na.rm = TRUE),
              total_genes = n(),
              prop_sig = sig_genes / total_genes)
  CountedPlot <- ggplot(Count, aes(x=category, y=prop_sig)) + geom_bar(stat="identity") + 
    labs(x="Length of Primary Transcript",  y="Proportion of Significant Genes") + theme(text = element_text(size=20))
  
  return(CountedPlot)
}

SummarisePadj <- function(Break){
  Count <- Break %>% 
    group_by(category) %>%
    summarise(sig_genes = sum(padj<0.05, na.rm = TRUE),
              total_genes = n(),
              prop_sig = sig_genes / total_genes)
  CountedPlot <- ggplot(Count, aes(x=category, y=prop_sig)) + geom_bar(stat="identity") + 
    labs(x="Length of Primary Transcript", y="Proportion of Significant Genes") + theme(text = element_text(size=20))
  return(CountedPlot)
}

d14QuantPvalCount <- SummarisePval(d14QuantBreak) + ggtitle("Proportion of significant genes  by primary \ntranscript length in QuantSeq d14") + theme_bw() + geom_bar(stat="identity", fill="#99CCFF")
d14TotalPvalCount <- SummarisePval(d14TotalBreak) + ggtitle("Proportion of significant genes  by primary \ntranscript length in Total RNA Seq d14") + theme_bw() + geom_bar(stat="identity", fill='#99FF99')

KOQuantPvalCount <- SummarisePval(KOQuantBreak) + ggtitle("Proportion of significant genes  by primary \ntranscript length in QuantSeq KO") + theme_bw() + geom_bar(stat="identity", fill="#0080FF")
KOTotalPvalCount <- SummarisePval(KOTotalBreak) + ggtitle("Proportion of significant genes  by primary \ntranscript length in Total RNA Seq KO") + theme_bw() + geom_bar(stat="identity", fill='#33FF33')

d14QuantPvalCount + d14TotalPvalCount + KOQuantPvalCount + KOTotalPvalCount

d14QuantPadjCount <- SummarisePadj(d14QuantBreak) + ggtitle("Proportion of significant genes  by primary \n         transcript length in QuantSeq d14")
d14TotalPadjCount <- SummarisePadj(d14TotalBreak) + ggtitle("Proportion of significant genes  by primary \n         transcript length in Total RNA Seq d14")

KOQuantPadjCount <- SummarisePadj(KOQuantBreak) + ggtitle("Proportion of significant genes  by primary \n         transcript length in QuantSeq KO")
KOTotalPadjCount <- SummarisePadj(KOTotalBreak) + ggtitle("Proportion of significant genes  by primary \n         transcript length in otal RNA Seq KO")

d14QuantPadjCount + d14TotalPadjCount + KOQuantPadjCount + KOTotalPadjCount
