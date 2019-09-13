library(tidyverse)
library(ggplot2)

MakeSignedZscore <- function(Diff){
  signedZ <- Diff %>%
    mutate(ZScore = ifelse(test = Diff$log2FoldChange > 0,
                           yes = qnorm(1 - (Diff$pvalue / 2)),
                           no = qnorm(Diff$pvalue / 2) )) %>%
    na.omit
  # high pass - otherwise very low P is set to Inf
  signedZ$ZScore[signedZ$ZScore > 6] <- 6
  signedZ$ZScore[signedZ$ZScore < -6] <- -6
  signedZ <- select(signedZ, EnsemblID, ZScore)
  return(signedZ)
}

produceZscoreTable <- function(QuantSEQ, TotalSEQ){
  QuantZ <- MakeSignedZscore(QuantSEQ)
  TotalZ <- MakeSignedZscore(TotalSEQ)
  Zscores <- full_join(TotalZ, QuantZ, by = "EnsemblID", suffix = c(".RNASeq", ".Quant")) %>%
    na.omit %>%
    column_to_rownames("EnsemblID")
  return(Zscores)
}

d14Sampled <- read_csv("DESeqD14RNASeqSampled.csv")
KOSampled <- read_csv("DESeqKORNASeqSampled.csv") 
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")

d14SampledZscore <- produceZscoreTable(d14Quant, d14Sampled)
KOSampledZscore <- produceZscoreTable(KOQuant, KOSampled)

d14 <- ggplot(d14SampledZscore, aes(x=ZScore.RNASeq, y=ZScore.Quant)) + geom_point(alpha=.2) + theme_bw() + 
  labs(title="FUS ctrl vs HOM d14 ZScore using downsampled reads", x="RNA-seq", y="QuantSeq") + 
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + 
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) +
  theme(text = element_text(size=20)) 

KO <- ggplot(KOSampledZscore, aes(x=ZScore.RNASeq, y=ZScore.Quant)) +   geom_point(alpha=.2) + theme_bw() + 
  labs(title="FUS ctrl vs HOM KO ZScore using downsampled reads", x="RNA-seq", y="QuantSeq") + 
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + 
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) +
  theme(text = element_text(size=20))

d14 + KO

