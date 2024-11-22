setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")



Matrices <- import.matrice(type = "CNV", freq = 0.01, type.recouvrement = "reciproque")
### Count transcripts in deletion
genematrix <- Matrices$DelMatrice

# AD
count.cnv <- apply(subset(genematrix, Status %in% c("EOAD", "LOAD"))[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
count.available <- apply(subset(genematrix, Status %in% c("EOAD", "LOAD"))[,-c(1:nbNonGenes)], 2, function(x){sum(!is.na(x))})
res.DEL.AD <- merge(data.frame(transcript.name = names(count.cnv), N.DEL.AD = count.cnv), data.frame(transcript.name = names(count.available), N.DEL.AD.available = count.available))


# CTRL
count.cnv <- apply(subset(genematrix, Status %in% c("Control"))[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
count.available <- apply(subset(genematrix, Status %in% c("Control"))[,-c(1:nbNonGenes)], 2, function(x){sum(!is.na(x))})
res.DEL.CTRL <- merge(data.frame(transcript.name = names(count.cnv), N.DEL.CTRL = count.cnv), data.frame(transcript.name = names(count.available), N.DEL.CTRL.available = count.available))


res.DEL <- merge(res.DEL.AD, res.DEL.CTRL, by = "transcript.name")


### Count transcripts in duplication
genematrix <- Matrices$DupFullMatrice

# AD
count.cnv <- apply(subset(genematrix, Status %in% c("EOAD", "LOAD"))[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
count.available <- apply(subset(genematrix, Status %in% c("EOAD", "LOAD"))[,-c(1:nbNonGenes)], 2, function(x){sum(!is.na(x))})
res.DUP.AD <- merge(data.frame(transcript.name = names(count.cnv), N.DUP.AD = count.cnv), data.frame(transcript.name = names(count.available), N.DUP.AD.available = count.available))


# CTRL
count.cnv <- apply(subset(genematrix, Status %in% c("Control"))[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
count.available <- apply(subset(genematrix, Status %in% c("Control"))[,-c(1:nbNonGenes)], 2, function(x){sum(!is.na(x))})
res.DUP.CTRL <- merge(data.frame(transcript.name = names(count.cnv), N.DUP.CTRL = count.cnv), data.frame(transcript.name = names(count.available), N.DUP.CTRL.available = count.available))


res.DUP <- merge(res.DUP.AD, res.DUP.CTRL, by = "transcript.name")


### merge all data

res <- merge(res.DEL, res.DUP, by = "transcript.name")

res.subset.available <- subset(res, N.DEL.AD.available > 0 & N.DEL.CTRL.available > 0 & N.DUP.AD.available > 0 & N.DUP.CTRL.available > 0)

res.subset.mirror <- subset(res.subset.available, (N.DEL.AD > 0 & N.DEL.CTRL == 0 & N.DUP.AD == 0 & N.DUP.CTRL > 0) | (N.DEL.AD == 0 & N.DEL.CTRL > 0 & N.DUP.AD > 0 & N.DUP.CTRL == 0))
res.subset.mirror$N.total.cnv <- res.subset.mirror$N.DEL.AD + res.subset.mirror$N.DEL.CTRL + res.subset.mirror$N.DUP.AD  + res.subset.mirror$N.DUP.CTRL




# load correspondance genes-transcriptsMatrices$
geneTranscript <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/GeneTranscript.Rdata")
gene1 <- GeneTranscript[!duplicated(GeneTranscript$Transcript),]
gene2 <- GeneTranscript[duplicated(GeneTranscript$Transcript),]
colnames(gene2) <- c("Gene 2", "Transcript")
GeneTranscript2 <- merge(gene1, gene2, all.x = TRUE) # because sometimes 2 genes are associated with one transcript

temp <- merge(res.subset.mirror, GeneTranscript2, by.x = "transcript.name", by.y = "Transcript", all.x = TRUE, all.y = FALSE)
res.mirror.order <- temp[order(temp$N.total.cnv, decreasing = TRUE),]


# save results
write.table(res.mirror.order, "RESULTS/mirror_CNV.txt", sep = "\t", col.names = TRUE, row.names = FALSE)


