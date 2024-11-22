transcript.selection <- function(freq, type.recouvrement) {
  
  load("../gnomadGeneFreq.RData")
  
  Matrices <- import.matrice(type = "CNV", freq = freq, type.recouvrement = type.recouvrement) 
  ### Count transcripts in deletion
  genematrix <- Matrices$DelMatrice
  
  count.cnv <- apply(genematrix[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
  data.count.cnv <- data.frame(transcript.name = names(count.cnv), count.cnv = count.cnv)
  
  ### list of excluded transcript based on frequency
  list.excluded.transcript <- unique(subset(freqGene, type == "DEL" & exclusion == TRUE)$Transcript)
  data.count.cnv$exclusion <- ifelse(data.count.cnv$transcript.name %in% list.excluded.transcript, TRUE, FALSE)
  data.count.cnv.not.excluded <- subset(data.count.cnv, exclusion == FALSE)
  
  write.table(data.count.cnv, 
              paste0("DATA/papier23octobre_count_CNV_DEL_filtre_freq_", freq, "_recouvrement_", type.recouvrement, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  count.cnv1 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv > 0 & data.count.cnv.not.excluded$count.cnv < 3]
  count.cnv3 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv >= 3]
  
  res.del <- data.frame(Dataset = "DelMatrice", N.transcript = length(data.count.cnv.not.excluded$count.cnv), transcript.0.N = sum(data.count.cnv.not.excluded$count.cnv == 0), transcript.0.percent = sum(data.count.cnv.not.excluded$count.cnv == 0)/length(data.count.cnv.not.excluded$count.cnv),
                        transcript.1.2.N = length(count.cnv1), transcript.1.2.percent = length(count.cnv1)/length(data.count.cnv.not.excluded$count.cnv),
                        transcript.3.N = length(count.cnv3), transcript.3.percent = length(count.cnv3)/length(data.count.cnv.not.excluded$count.cnv),
                        transcript.mean = mean(data.count.cnv.not.excluded$count.cnv), transcript.sd = sd(data.count.cnv.not.excluded$count.cnv), transcript.min = min(data.count.cnv.not.excluded$count.cnv), transcript.max = max(data.count.cnv.not.excluded$count.cnv),
                        transcript.median = median(data.count.cnv.not.excluded$count.cnv), transcript.Q1 = summary(data.count.cnv.not.excluded$count.cnv)[[2]], transcript.Q3 = summary(data.count.cnv.not.excluded$count.cnv)[[5]],
                        transcript.3.mean = mean(count.cnv3), transcript.3.sd = sd(count.cnv3), transcript.3.min = min(count.cnv3), transcript.3.max = max(count.cnv3),
                        transcript.3.median = median(count.cnv3), transcript.3.Q1 = summary(count.cnv3)[[2]], transcript.3.Q3 = summary(count.cnv3)[[5]],
                        stringsAsFactors = FALSE)
  
  
  ### Count transcripts in duplication
  genematrix <- Matrices$DupFullMatrice
  
  count.cnv <- apply(genematrix[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
  data.count.cnv <- data.frame(transcript.name = names(count.cnv), count.cnv = count.cnv)
  
  ### list of excluded transcript based on frequency
  list.excluded.transcript <- unique(subset(freqGene, type == "completeDUP" & exclusion == TRUE)$Transcript)
  data.count.cnv$exclusion <- ifelse(data.count.cnv$transcript.name %in% list.excluded.transcript, TRUE, FALSE)
  data.count.cnv.not.excluded <- subset(data.count.cnv, exclusion == FALSE)
  
  write.table(data.count.cnv, 
              paste0("DATA/papier23octobre_count_CNV_DUP_filtre_freq_", freq, "_recouvrement_", type.recouvrement, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  count.cnv1 <- count.cnv[count.cnv > 0 & count.cnv < 3]
  count.cnv3 <- count.cnv[count.cnv >= 3]
  
  
  count.cnv1 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv > 0 & data.count.cnv.not.excluded$count.cnv < 3]
  count.cnv3 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv >= 3]
  
  res.dup <- data.frame(Dataset = "DupFullMatrice", N.transcript = length(data.count.cnv.not.excluded$count.cnv), transcript.0.N = sum(data.count.cnv.not.excluded$count.cnv == 0), transcript.0.percent = sum(data.count.cnv.not.excluded$count.cnv == 0)/length(data.count.cnv.not.excluded$count.cnv),
                        transcript.1.2.N = length(count.cnv1), transcript.1.2.percent = length(count.cnv1)/length(data.count.cnv.not.excluded$count.cnv),
                        transcript.3.N = length(count.cnv3), transcript.3.percent = length(count.cnv3)/length(data.count.cnv.not.excluded$count.cnv),
                        transcript.mean = mean(data.count.cnv.not.excluded$count.cnv), transcript.sd = sd(data.count.cnv.not.excluded$count.cnv), transcript.min = min(data.count.cnv.not.excluded$count.cnv), transcript.max = max(data.count.cnv.not.excluded$count.cnv),
                        transcript.median = median(data.count.cnv.not.excluded$count.cnv), transcript.Q1 = summary(data.count.cnv.not.excluded$count.cnv)[[2]], transcript.Q3 = summary(data.count.cnv.not.excluded$count.cnv)[[5]],
                        transcript.3.mean = mean(count.cnv3), transcript.3.sd = sd(count.cnv3), transcript.3.min = min(count.cnv3), transcript.3.max = max(count.cnv3),
                        transcript.3.median = median(count.cnv3), transcript.3.Q1 = summary(count.cnv3)[[2]], transcript.3.Q3 = summary(count.cnv3)[[5]],
                        stringsAsFactors = FALSE)
  
  ### Count transcripts in deletion + partial duplication
  genematrix <- Matrices$DelMatriceWithPartial
  
  count.cnv <- apply(genematrix[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
  data.count.cnv <- data.frame(transcript.name = names(count.cnv), count.cnv = count.cnv)
  
  ### list of excluded transcript based on frequency (DEL = partialDUP = same filtre than for DEL)
  list.excluded.transcript <- unique(subset(freqGene, type == "DEL" & exclusion == TRUE)$Transcript)
  data.count.cnv$exclusion <- ifelse(data.count.cnv$transcript.name %in% list.excluded.transcript, TRUE, FALSE)
  data.count.cnv.not.excluded <- subset(data.count.cnv, exclusion == FALSE)
  
  write.table(data.count.cnv, 
              paste0("DATA/papier23octobre_count_CNV_DEL_and_DUPpartial_filtre_freq_", freq, "_recouvrement_", type.recouvrement, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  count.cnv1 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv > 0 & data.count.cnv.not.excluded$count.cnv < 3]
  count.cnv3 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv >= 3]
  
  res.delduppartial <- data.frame(Dataset = "DelMatriceWithPartial", N.transcript = length(data.count.cnv.not.excluded$count.cnv), transcript.0.N = sum(data.count.cnv.not.excluded$count.cnv == 0), transcript.0.percent = sum(data.count.cnv.not.excluded$count.cnv == 0)/length(data.count.cnv.not.excluded$count.cnv),
                                  transcript.1.2.N = length(count.cnv1), transcript.1.2.percent = length(count.cnv1)/length(data.count.cnv.not.excluded$count.cnv),
                                  transcript.3.N = length(count.cnv3), transcript.3.percent = length(count.cnv3)/length(data.count.cnv.not.excluded$count.cnv),
                                  transcript.mean = mean(data.count.cnv.not.excluded$count.cnv), transcript.sd = sd(data.count.cnv.not.excluded$count.cnv), transcript.min = min(data.count.cnv.not.excluded$count.cnv), transcript.max = max(data.count.cnv.not.excluded$count.cnv),
                                  transcript.median = median(data.count.cnv.not.excluded$count.cnv), transcript.Q1 = summary(data.count.cnv.not.excluded$count.cnv)[[2]], transcript.Q3 = summary(data.count.cnv.not.excluded$count.cnv)[[5]],
                                  transcript.3.mean = mean(count.cnv3), transcript.3.sd = sd(count.cnv3), transcript.3.min = min(count.cnv3), transcript.3.max = max(count.cnv3),
                                  transcript.3.median = median(count.cnv3), transcript.3.Q1 = summary(count.cnv3)[[2]], transcript.3.Q3 = summary(count.cnv3)[[5]],
                                  stringsAsFactors = FALSE)
  
  
  ### Count transcripts in partial duplication
  genematrix <- Matrices$DupPartialMatrice
  
  count.cnv <- apply(genematrix[,-c(1:nbNonGenes)], 2, function(x){sum(x, na.rm = TRUE)})
  data.count.cnv <- data.frame(transcript.name = names(count.cnv), count.cnv = count.cnv)
  
  ### list of excluded transcript based on frequency (DEL = partialDUP = same filtre than for DEL)
  list.excluded.transcript <- unique(subset(freqGene, type == "DEL" & exclusion == TRUE)$Transcript)
  data.count.cnv$exclusion <- ifelse(data.count.cnv$transcript.name %in% list.excluded.transcript, TRUE, FALSE)
  data.count.cnv.not.excluded <- subset(data.count.cnv, exclusion == FALSE)
  
  write.table(data.count.cnv, 
              paste0("DATA/papier23octobre_count_CNV_DUPpartial_filtre_freq_", freq, "_recouvrement_", type.recouvrement, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  count.cnv1 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv > 0 & data.count.cnv.not.excluded$count.cnv < 3]
  count.cnv3 <- data.count.cnv.not.excluded$count.cnv[data.count.cnv.not.excluded$count.cnv >= 3]
  
  res.duppartial <- data.frame(Dataset = "DupPartial", N.transcript = length(data.count.cnv.not.excluded$count.cnv), transcript.0.N = sum(data.count.cnv.not.excluded$count.cnv == 0), transcript.0.percent = sum(data.count.cnv.not.excluded$count.cnv == 0)/length(data.count.cnv.not.excluded$count.cnv),
                                  transcript.1.2.N = length(count.cnv1), transcript.1.2.percent = length(count.cnv1)/length(data.count.cnv.not.excluded$count.cnv),
                                  transcript.3.N = length(count.cnv3), transcript.3.percent = length(count.cnv3)/length(data.count.cnv.not.excluded$count.cnv),
                                  transcript.mean = mean(data.count.cnv.not.excluded$count.cnv), transcript.sd = sd(data.count.cnv.not.excluded$count.cnv), transcript.min = min(data.count.cnv.not.excluded$count.cnv), transcript.max = max(data.count.cnv.not.excluded$count.cnv),
                                  transcript.median = median(data.count.cnv.not.excluded$count.cnv), transcript.Q1 = summary(data.count.cnv.not.excluded$count.cnv)[[2]], transcript.Q3 = summary(data.count.cnv.not.excluded$count.cnv)[[5]],
                                  transcript.3.mean = mean(count.cnv3), transcript.3.sd = sd(count.cnv3), transcript.3.min = min(count.cnv3), transcript.3.max = max(count.cnv3),
                                  transcript.3.median = median(count.cnv3), transcript.3.Q1 = summary(count.cnv3)[[2]], transcript.3.Q3 = summary(count.cnv3)[[5]],
                                  stringsAsFactors = FALSE)
  
  
  res <- rbind(res.del, res.dup, res.delduppartial, res.duppartial)
  
  write.table(res,  paste0("RESULTS/papier23octobre_Filtre_transcript_description_filtre_freq_", freq, "_recouvrement_", type.recouvrement, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
}

