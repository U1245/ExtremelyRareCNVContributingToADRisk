# count transcript in dosage union set A et B analysis
# dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)
# sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >=4)
dat <- read.table("RESULTS/Table_for_article_october2023_Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)
sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" &  Ncarriers.CNV.total >=4)
dim(sub)
length(unique(sub$Gene))

# count transcript in set A and B before filtres

filtre.freq.value <- 0.01
filtre.freq.type.recouvrement<- "reciproque"
# filtreDEL <- read.table(paste0("DATA/count_CNV_DEL_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
# filtreDUP <- read.table(paste0("DATA/count_CNV_DUP_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)

filtreDEL <- read.table(paste0("DATA/papier23octobre_count_CNV_DEL_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
filtreDUP <- read.table(paste0("DATA/papier23octobre_count_CNV_DUP_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)

filtreCNV <- merge(filtreDEL, filtreDUP, by = "transcript.name")
filtreCNV$typeNM <- substr(filtreCNV$transcript.name, 1, 2)
# load correspondance genes-transcripts
geneTranscript <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/GeneTranscript.Rdata")
gene1 <- GeneTranscript[!duplicated(GeneTranscript$Transcript),]
gene2 <- GeneTranscript[duplicated(GeneTranscript$Transcript),]
colnames(gene2) <- c("Gene 2", "Transcript")
GeneTranscript2 <- merge(gene1, gene2, all.x = TRUE) # because sometimes 2 genes are associated with one transcript
dat <- merge(filtreCNV, GeneTranscript2, by.x = "transcript.name", by.y = "Transcript", all.x = TRUE, all.y = FALSE)

# count nb transcripts
sub <- subset(dat, typeNM == "NM")
dim(sub)
length(unique(sub$Gene))


# set A
subA <- subset(dat, typeNM == "NM" & exclusion.x == FALSE)
dim(subA)
length(unique(subA$Gene))
subA2 <- subset(dat, typeNM == "NM" & exclusion.x == FALSE & count.cnv.x >= 4)
dim(subA2)
length(unique(subA2$Gene))

# set B
subB <- subset(dat, typeNM == "NM" & exclusion.y == FALSE)
dim(subB)
length(unique(subB$Gene))
subB2 <- subset(dat, typeNM == "NM" & exclusion.y == FALSE & count.cnv.y >= 4)
dim(subB2)
length(unique(subB2$Gene))

# union sets A and B
sub <- subset(dat, typeNM == "NM" & (exclusion.x == FALSE | exclusion.y == FALSE))
dim(sub)
length(unique(sub$Gene))
dat$sumCNV <- dat$count.cnv.x + dat$count.cnv.y
sub2 <- subset(dat, typeNM == "NM" & (exclusion.x == FALSE | exclusion.y == FALSE) & sumCNV >= 4 & count.cnv.x >= 1 & count.cnv.y >= 1)
dim(sub2)
length(unique(sub2$Gene))
