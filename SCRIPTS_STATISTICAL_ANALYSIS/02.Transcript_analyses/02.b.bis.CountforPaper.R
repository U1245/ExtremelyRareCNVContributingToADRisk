setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")


countforpaper_dosage.analysis_union_setA_B <- function(type.CNV, filtre.freq.value, filtre.freq.type.recouvrement) {
  
  
  Matrices <- import.matrice(type = "dosage", freq = filtre.freq.value, type.recouvrement = filtre.freq.type.recouvrement)
  
  
  # verifier noms matrices
  if (type.CNV == "DEL_DUP") {
    
    filtreDEL <- read.table(paste0("DATA/count_CNV_DEL_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreDUP <- read.table(paste0("DATA/count_CNV_DUP_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreCNV <- merge(filtreDEL, filtreDUP, by = "transcript.name")
    filtreCNV$count.total <- filtreCNV$count.cnv.x + filtreCNV$count.cnv.y
    transcript.list <- as.character(subset(filtreCNV, substr(transcript.name, 1, 2) == "NM")$transcript.name)
    
    genematrix <- Matrices$DosageMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  
  if (type.CNV == "DEL_and_DUPpartial_DUP") {
    
    filtreDEL <- read.table(paste0("DATA/count_CNV_DEL_and_DUPpartial_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreDUP <- read.table(paste0("DATA/count_CNV_DUP_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreCNV <- merge(filtreDEL, filtreDUP, by = "transcript.name")
    filtreCNV$count.total <- filtreCNV$count.cnv.x + filtreCNV$count.cnv.y
    transcript.list <- as.character(subset(filtreCNV, substr(transcript.name, 1, 2) == "NM")$transcript.name)
    
    genematrix <- Matrices$DosageMatriceWithPartial[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  
  B <- length(transcript.list)
  # B <- 20
  cl <- parallel::makeCluster(10)
  doParallel::registerDoParallel(cl)
  
  resultsB <- foreach(iter = 1:B, .combine = "rbind", .packages = c("logistf", "brglm2")) %dopar% {
    # sink("log.txt", append=TRUE)
    # cat(paste("Starting iteration", iter, "\n"))
    cngeneral <- c("transcript.name", "N CTRL", "N LOAD", "N EOAD", "N Total", 
                   "Ncarriers CNV CTRL", "Ncarriers CNV LOAD", "Ncarriers CNV EOAD", "Ncarriers CNV Total", 
                   "Ncarriers DEL CTRL", "Ncarriers DEL LOAD", "Ncarriers DEL EOAD", "Ncarriers DEL Total",
                   "Ncarriers DUP CTRL", "Ncarriers DUP LOAD", "Ncarriers DUP EOAD", "Ncarriers DUP Total")
  
    ### WITHOUT APOE ADJUSTMENT
    
    ### Data for the current transcript (iter)
    
    dat <- genematrix[, c("sample", "Status", transcript.list[iter])]
    dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
    dat$cnv <- dat[,transcript.list[iter]]
    
    
    count <- addmargins(table(dat$cnv, dat$Status.f))
    
    cnv.N <- c(count["Sum",], count["Sum",]-count["2",])
    cnv.rate <- c((count["Sum",] - count["2",])/count["Sum",])
    # count DEL
    countDEL <- count[c(na.omit(as.numeric(rownames(count)) <= 1), FALSE), ]
    if (is.null(dim(countDEL))) {
      cnvDEL.N <- countDEL
    } else {
      cnvDEL.N <- apply(countDEL, 2, sum) 
    }
    cnvDEL.rate <- cnvDEL.N/count["Sum",]
    # count DUP
    countDUP <- count[c(na.omit(as.numeric(rownames(count)) >= 3), FALSE), ]
    if (is.null(dim(countDUP))) {
      cnvDUP.N <- countDUP
    } else {
      cnvDUP.N <- apply(countDUP, 2, sum) 
    }
    cnvDUP.rate <- cnvDUP.N/count["Sum",]
    

    res <- c(cnv.N, cnvDEL.N, cnvDUP.N)
      
    
    res.withoutAdjustment <- data.frame(transcript.name = transcript.list[iter], t(res))
    colnames(res.withoutAdjustment) <- cngeneral
    
   
    
    return(res.withoutAdjustment)
    
  }
  
  stopCluster(cl)
  
  # load correspondance genes-transcripts
  geneTranscript <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/GeneTranscript.Rdata")
  gene1 <- GeneTranscript[!duplicated(GeneTranscript$Transcript),]
  gene2 <- GeneTranscript[duplicated(GeneTranscript$Transcript),]
  colnames(gene2) <- c("Gene 2", "Transcript")
  GeneTranscript2 <- merge(gene1, gene2, all.x = TRUE) # because sometimes 2 genes are associated with one transcript
  
  temp <- merge(resultsB, GeneTranscript2, by.x = "transcript.name", by.y = "Transcript", all.x = TRUE, all.y = FALSE)
  
  # save results
  write.table(temp, paste0("RESULTS/Table_for_article_october2023_Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
}

