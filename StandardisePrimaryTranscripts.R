library(tidyverse)
library(ggplot2)
library(biomaRt)

#Choose database
mart <- useMart("ensembl")
#Choose genome
mart <- useDataset("mmusculus_gene_ensembl", mart)

mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "mmusculus_gene_ensembl", 
              host = "www.ensembl.org")
#Set columns you want
attributes <- c("ensembl_transcript_id", "transcript_appris", "ensembl_gene_id", "start_position", "end_position")
#Set column to filter by
filters <- "ensembl_gene_id"

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")

#Join Files Together
d14Join <- full_join(d14Unsampled, d14Quant, by = "EnsemblID") %>%
  filter(!grepl('+', EnsemblID, fixed = TRUE))
KOJoin <- full_join(KOUnsampled, KOQuant, by = "EnsemblID")%>%
  filter(!grepl('+', EnsemblID, fixed = TRUE))

#Get list of gene names
d14GeneID <- unlist(d14Join$EnsemblID)
KOGeneID <- unlist(KOJoin$EnsemblID)

#Get market file for both samples
d14output <- getBM(attributes = attributes, filters = filters, values=d14GeneID, mart=mart)
KOoutput <- getBM(attributes = attributes, filters = filters, values=KOGeneID, mart=mart)

#Calculate length of primary transcript
d14output$length <- d14output$end_position - d14output$start_position
KOoutput$length <- KOoutput$end_position - KOoutput$start_position

#Filter for only primary transcript
d14outfilt <- filter(d14output, transcript_appris == "principal1")
KOoutfilt <- filter(KOoutput, transcript_appris == "principal1")

write_csv(d14outfilt, "d14TranscriptsPrimary.csv")
write_csv(KOoutfilt, "KOTranscriptsPrimary.csv")

write.csv(d14output, "d14Transcripts.csv")
write.csv(KOoutput, "KOTranscripts.csv")

#Load in gencode gtf file
setwd("~/Downloads/")
gencode <- read_tsv("gencode.vM20.basic.annotation.gtf", skip=5, col_names = FALSE)
colnames(gencode) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
GencodeExon <- filter(gencode, feature == "exon")

#split column
GencodeExonSplit <- separate(GencodeExon, attribute, c("gene_id", "transcript_id", "gene_type", "transcript_type", "transcript_name", "exon_number", "exon_id", "level", "transcript_support_level", "tag", "havana_gene", "havana_transcript"), sep=";") %>%
  select(seqname:transcript_id)
GencodeExonSplit$geneid <- rm_between(GencodeExonSplit$gene_id, '"', '"', extract=TRUE)
GencodeExonSplit$transcriptid <- rm_between(GencodeExonSplit$transcript_id, '"', '"', extract=TRUE)
GencodeExonSplit <- select(GencodeExonSplit, seqname:frame, geneid, transcriptid)
GencodeSelected <- GencodeExonSplit

#Remove everything after the . in both gene and transcript IDS. I needed to use ..* rather than just .*
GencodeSelected$geneid <- str_replace(GencodeSelected$geneid ,"\\..*","")
GencodeSelected$transcriptid <- str_replace(GencodeSelected$transcriptid ,"\\..*","")


#Work out the length of each transcript
GencodeSelected$length  <- GencodeSelected$end - GencodeSelected$start 
Length <- GencodeSelected %>% group_by(transcriptid) %>% summarise(length = sum(length), exons = n())

setwd("~/Google Drive/Work/UCL/Seth/TestingQuantSeq")
write.csv(Length, "GencodeLength.csv", row.names = FALSE)