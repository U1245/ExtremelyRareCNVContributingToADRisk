load("../MatriceResultat/matrice.Rdata")
load("../CNV_postFiltre.Rdata")
# DELETIONS
dat.del <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DEL_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

dat <- subset(dat.del, p.value..ordinal.regression..using.polr. < 1e-3 & Adjustment == "none" & substr(transcript.name, 1, 2) == "NM")

# load("../MatriceResultat/matrice.Rdata")
# load("../CNV_postFiltre.Rdata")
Rouen <- c("ALZ_AUT.201808.Agilent_V5UTR", "ALZ.201308.Agilent_V5", "ALZ.201312.Agilent_V5", "ALZ.201405.Agilent_V5", "ALZ.201410.Agilent_V5", 
           "ALZ.201601.Agilent_V5", "ALZ.201604.Agilent_V5", "ALZ.201606.Agilent_V5", "ALZ.201702.Agilent_V5UTR", "ALZ.201709.Agilent_V5UTR", 
           "ALZ.201804.Agilent_V5UTR", "ALZ.201907.Agilent_V6UTR", "ALZ.ECASCAD51.Agilent_V6UTR", "ALZ.ECASCAD52.Agilent_V6UTR", 
           "ALZ.ECASCAD53_REDIA1.Agilent_V6UTR", "ALZ.ECASCAD54.Agilent_V6UTR", "ALZ.ECASCAD55.Agilent_V6UTR", "ALZ.ECASCAD56.Agilent_V6UTR", 
           "FREX_Rouen.Agilent_V5UTR") 

transcriptid <- unique(dat$transcript.name)
ddall <- DelMatrice[, c(colnames(DelMatrice)[1:13], transcriptid)]

ddrouen <- DelMatrice[DelMatrice$Cohort == "(FR) ADES-FR" & DelMatrice$datasetid %in% Rouen, c(colnames(DelMatrice)[1:13], transcriptid)]

i <- 1
# ddall[ddall[, transcriptid[i]] == 1,]
extractrouen <- ddrouen[ddrouen[, transcriptid[i]] == 1,]
cnvrouen <- subset(CNV_postFiltre, sample %in% extractrouen$sample & Type == "DEL", select = c(sample, Type, chrom, start, end, reference, ID, Full, Partial))
cnvrouen.extracttranscript <- cnvrouen[grep(transcriptid[i], paste(cnvrouen$Full, cnvrouen$Partial)), ]
cnvrouen.extracttranscript$Gene <- unique(subset(dat, transcript.name == transcriptid[i])$Gene)
cnvrouen.extracttranscript$Transcript <- transcriptid[i]
res <- cnvrouen.extracttranscript
for (i in 2:length(transcriptid)) {
  print(i)
  # ddall[ddall[, transcriptid[i]] == 1,]
  extractrouen <- ddrouen[ddrouen[, transcriptid[i]] == 1,]
  cnvrouen <- subset(CNV_postFiltre, sample %in% extractrouen$sample & Type == "DEL", select = c(sample, Type, chrom, start, end, reference, ID, Full, Partial))
  if (dim(cnvrouen)[1] == 0) {
    cnvrouen.extracttranscript <- data.frame(matrix(rep(NA, 9), ncol = 9))
    colnames(cnvrouen.extracttranscript) <- colnames(cnvrouen)
  } else {
    cnvrouen.extracttranscript <- cnvrouen[grep(transcriptid[i], paste(cnvrouen$Full, cnvrouen$Partial)), ] 
  }
  cnvrouen.extracttranscript$Gene <- unique(subset(dat, transcript.name == transcriptid[i])$Gene)
  cnvrouen.extracttranscript$Transcript <- transcriptid[i]
  # print(cnvrouen.extracttranscript)
  res <- rbind(res, cnvrouen.extracttranscript)
}

res.del <- subset(res, select = c(Gene, Transcript, sample, Type, chrom, start, end, reference, ID))
res.del.ord <- res.del[order(res.del$chrom, res.del$Gene, res.del$start, res.del$end),]

dat.dup <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

dat <- subset(dat.dup, p.value..ordinal.regression..using.polr. < 1e-3 & Adjustment == "none" & substr(transcript.name, 1, 2) == "NM")

# load("../MatriceResultat/matrice.Rdata")

Rouen <- c("ALZ_AUT.201808.Agilent_V5UTR", "ALZ.201308.Agilent_V5", "ALZ.201312.Agilent_V5", "ALZ.201405.Agilent_V5", "ALZ.201410.Agilent_V5", 
           "ALZ.201601.Agilent_V5", "ALZ.201604.Agilent_V5", "ALZ.201606.Agilent_V5", "ALZ.201702.Agilent_V5UTR", "ALZ.201709.Agilent_V5UTR", 
           "ALZ.201804.Agilent_V5UTR", "ALZ.201907.Agilent_V6UTR", "ALZ.ECASCAD51.Agilent_V6UTR", "ALZ.ECASCAD52.Agilent_V6UTR", 
           "ALZ.ECASCAD53_REDIA1.Agilent_V6UTR", "ALZ.ECASCAD54.Agilent_V6UTR", "ALZ.ECASCAD55.Agilent_V6UTR", "ALZ.ECASCAD56.Agilent_V6UTR", 
           "FREX_Rouen.Agilent_V5UTR") 

transcriptid <- unique(dat$transcript.name)
ddall <- DupFullMatrice[, c(colnames(DupFullMatrice)[1:13], transcriptid)]

ddrouen <- DupFullMatrice[DelMatrice$Cohort == "(FR) ADES-FR" & DupFullMatrice$datasetid %in% Rouen, c(colnames(DupFullMatrice)[1:13], transcriptid)]

i <- 1
# ddall[ddall[, transcriptid[i]] == 1,]
extractrouen <- ddrouen[ddrouen[, transcriptid[i]] == 1,]
cnvrouen <- subset(CNV_postFiltre, sample %in% extractrouen$sample & Type == "DUP", select = c(sample, Type, chrom, start, end, reference, ID, Full))
if (dim(cnvrouen)[1] == 0) {
  cnvrouen.extracttranscript <- data.frame(matrix(rep(NA, 8), ncol = 8))
  colnames(cnvrouen.extracttranscript) <- colnames(cnvrouen)
} else {
  cnvrouen.extracttranscript <- cnvrouen[grep(transcriptid[i], cnvrouen$Full), ] 
}
cnvrouen.extracttranscript$Gene <- unique(subset(dat, transcript.name == transcriptid[i])$Gene)
cnvrouen.extracttranscript$Transcript <- transcriptid[i]
res <- cnvrouen.extracttranscript
for (i in 2:length(transcriptid)) {
  print(i)
  # ddall[ddall[, transcriptid[i]] == 1,]
  extractrouen <- ddrouen[ddrouen[, transcriptid[i]] == 1,]
  cnvrouen <- subset(CNV_postFiltre, sample %in% extractrouen$sample & Type == "DUP", select = c(sample, Type, chrom, start, end, reference, ID, Full))
  if (dim(cnvrouen)[1] == 0) {
    cnvrouen.extracttranscript <- data.frame(matrix(rep(NA, 8), ncol = 8))
    colnames(cnvrouen.extracttranscript) <- colnames(cnvrouen)
  } else {
    cnvrouen.extracttranscript <- cnvrouen[grep(transcriptid[i], cnvrouen$Full), ] 
  }
  cnvrouen.extracttranscript$Gene <- unique(subset(dat, transcript.name == transcriptid[i])$Gene)
  cnvrouen.extracttranscript$Transcript <- transcriptid[i]
  # print(cnvrouen.extracttranscript)
  res <- rbind(res, cnvrouen.extracttranscript)
}


res.dup <- subset(res, select = c(Gene, Transcript, sample, Type, chrom, start, end, reference, ID))
res.dup.ord <- res.dup[order(res.dup$chrom, res.dup$Gene, res.dup$start, res.dup$end),]


write.table(res.del.ord, "DATA/list_CNV_DEL_tocheck_Rouen.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(res.dup.ord, "DATA/list_CNV_DUP_tocheck_Rouen.txt", col.names = TRUE, row.names = FALSE, sep = "\t")