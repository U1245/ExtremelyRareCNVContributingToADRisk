setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")


dosage.analysis_union_setA_B <- function(type.CNV, filtre.freq.value, filtre.freq.type.recouvrement) {
  
  
  Matrices <- import.matrice(type = "dosage", freq = filtre.freq.value, type.recouvrement = filtre.freq.type.recouvrement)
  
   
  # verifier noms matrices
  if (type.CNV == "DEL_DUP") {
    
    filtreDEL <- read.table(paste0("DATA/count_CNV_DEL_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreDUP <- read.table(paste0("DATA/count_CNV_DUP_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreCNV <- merge(filtreDEL, filtreDUP, by = "transcript.name")
    filtreCNV$count.total <- filtreCNV$count.cnv.x + filtreCNV$count.cnv.y
    transcript.list <- as.character(subset(filtreCNV, count.cnv.x >= 1 & count.cnv.y >= 1 & count.total >= 3 & (exclusion.x == FALSE | exclusion.y == FALSE))$transcript.name)
   
    genematrix <- Matrices$DosageMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }

  if (type.CNV == "DEL_and_DUPpartial_DUP") {
    
    filtreDEL <- read.table(paste0("DATA/count_CNV_DEL_and_DUPpartial_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreDUP <- read.table(paste0("DATA/count_CNV_DUP_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
    filtreCNV <- merge(filtreDEL, filtreDUP, by = "transcript.name")
    filtreCNV$count.total <- filtreCNV$count.cnv.x + filtreCNV$count.cnv.y
    transcript.list <- as.character(subset(filtreCNV, count.cnv.x >= 1 & count.cnv.y >= 1 & count.total >= 3 & (exclusion.x == FALSE | exclusion.y == FALSE))$transcript.name)
    
    genematrix <- Matrices$DosageMatriceWithPartial[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  
  B <- length(transcript.list)
# B <- 20
  cl <- parallel::makeCluster(30)
  doParallel::registerDoParallel(cl)
  
  resultsB <- foreach(iter = 1:B, .combine = "rbind", .packages = c("logistf", "brglm2")) %dopar% {
    # sink("log.txt", append=TRUE)
    # cat(paste("Starting iteration", iter, "\n"))
    cngeneral <- c("transcript.name", "Adjustment", "N CTRL", "N LOAD", "N EOAD", "N total", 
                   "Ncarriers CNV CTRL", "Ncarriers CNV LOAD", "Ncarriers CNV EOAD", "Ncarriers CNV total", 
                   "Ncarriers DEL CTRL", "Ncarriers DEL LOAD", "Ncarriers DEL EOAD", "Ncarriers DEL total", 
                   "Ncarriers DUP CTRL", "Ncarriers DUP LOAD", "Ncarriers DUP EOAD", "Ncarriers DUP total", 
                   "Rate CNV CTRL", "Rate CNV LOAD", "Rate CNV EOAD", "Rate CNV total", 
                   "Rate DEL CTRL", "Rate DEL LOAD", "Rate DEL EOAD", "Rate DEL total", 
                   "Rate DUP CTRL", "Rate DUP LOAD", "Rate DUP EOAD", "Rate DUP total", 
                   "OR all AD vs CTRL", "OR all AD vs CTRL - CIinf", "OR all AD vs CTRL - CIsup", "p-value all AD vs CTRL", 
                   "OR firth all AD vs CTRL", "OR firth all AD vs CTRL - CIinf", "OR firth all AD vs CTRL - CIsup", "p-value firth all AD vs CTRL", 
                   "OR EOAD vs CTRL", "OR EOAD vs CTRL - CIinf", "OR EOAD vs CTRL - CIsup", "p-value EOAD vs CTRL",
                   "OR LOAD vs CTRL", "OR LOAD vs CTRL - CIinf", "OR LOAD vs CTRL - CIsup", "p-value LOAD vs CTRL",
                   "OR firth EOAD vs CTRL", "OR firth EOAD vs CTRL - CIinf", "OR firth EOAD vs CTRL - CIsup", "p-value firth EOAD vs CTRL", 
                   "OR firth LOAD vs CTRL", "OR firth LOAD vs CTRL - CIinf", "OR firth LOAD vs CTRL - CIsup", "p-value firth LOAD vs CTRL",  
                   "OR clm", "OR clm - CIinf", "OR clm - CIsup", "OR polr", "OR polr - CIinf", "OR polr - CIsup", "p-value (ordinal regression, using clm)", "p-value (ordinal regression, using polr)", "warning clm")
    
    library(ordinal)
    library(MASS)
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
    
    if (0 %in% count["Sum",]) {
      res <- c(cnv.N, cnvDEL.N, cnvDUP.N, cnv.rate, cnvDEL.rate, cnvDUP.rate, rep(NA, 31))
    } else {
      # cat(paste("AD", iter, "\n"))
      # OR AD
      dat.AD <- dat
      dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
      m <- glm(Y ~ cnv, data = dat.AD, family = binomial(link = "logit"))
      # print(summary(m))
      # print(exp(m$coef))
      # print(confint(m))
      # print(summary(m)$coefficients)
      OR.AD <- c(exp(m$coef)["cnv"], exp(confint(m))["cnv",])
      
      pvalue.AD <- summary(m)$coefficients["cnv", "Pr(>|z|)"]
      
      # cat(iter, " LOGISTF AD \n")
      mf <- logistf(Y ~ cnv, data = dat.AD)
      print(summary(mf))
      ORf.AD <- c(exp(mf$coefficients["cnv"]), exp(confint(mf))["cnv",])
      pvaluef.AD <- mf$prob["cnv"]
      
      # OR EOAD
      dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
      if (sum(dat.EOAD$cnv != 2, na.rm = TRUE) == 0) {
        OR.EOAD <- c(NA, NA, NA)
        ORf.EOAD <- c(NA, NA, NA)
        pvalue.EOAD <- NA
        pvaluef.EOAD <- NA
      } else {
        dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
        m <- glm(Y ~ cnv, data = dat.EOAD, family = binomial(link = "logit"))
        OR.EOAD <- c(exp(m$coef)["cnv"], exp(confint(m))["cnv",])
        pvalue.EOAD <- summary(m)$coefficients["cnv", "Pr(>|z|)"]
        # cat(iter, " LOGISTF EOAD \n")
        mf <- logistf(Y ~ cnv, data = dat.EOAD)
        ORf.EOAD <- c(exp(mf$coefficients["cnv"]), exp(confint(mf))["cnv",])
        pvaluef.EOAD <- mf$prob["cnv"]
      }
      
      # OR LOAD
      dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
      if (sum(dat.LOAD$cnv != 2, na.rm = TRUE) == 0) {
        OR.LOAD <- c(NA, NA, NA)
        ORf.LOAD <- c(NA, NA, NA)
        pvalue.LOAD <- NA
        pvaluef.LOAD <- NA
      } else {
        dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
        m <- glm(Y ~ cnv, data = dat.LOAD, family = binomial(link = "logit"))
        OR.LOAD <- c(exp(m$coef)["cnv"], exp(confint(m))["cnv",])
        pvalue.LOAD <- summary(m)$coefficients["cnv", "Pr(>|z|)"]
        # cat(iter, " LOGISTF LOAD \n")
        mf <- logistf(Y ~ cnv, data = dat.LOAD)
        ORf.LOAD <- c(exp(mf$coefficients["cnv"]), exp(confint(mf))["cnv",])
        pvaluef.LOAD <- mf$prob["cnv"]
      }
      # cat(paste("clm", iter, "\n"))
      # Ordinal regression avec clm
      rego <- clm(Status.f ~ cnv, data = na.omit(dat))
      confint.rego <- try(exp(confint(rego))["cnv",], silent = TRUE)
      if (inherits(confint.rego, "try-error")) {
        OR.clm <- c(exp(rego$coefficients["cnv"]), NA, NA)
      } else {
        OR.clm <- c(exp(rego$coefficients["cnv"]), confint.rego)
      }
      
      pvalue.clm <- summary(rego)$coefficients["cnv", "Pr(>|z|)"]
      message.clm <- rego$message
      
      # Ordinal regression avec polr
      rego2 <- polr(Status.f ~ cnv, data = na.omit(dat), Hess = TRUE)
      confint.rego2 <- try(exp(confint(rego2)), silent = TRUE)
      if (inherits(confint.rego2, "try-error")) {
        OR.polr <- c(exp(rego2$coefficients["cnv"]), NA, NA)
      } else {
        OR.polr <- c(exp(rego2$coefficients["cnv"]), confint.rego2)
      }
      
      rego0 <- polr(Status.f ~ 1, data = na.omit(dat))
      statRV <- -2*(logLik(rego0) - logLik(rego2))
      pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
      
      res <- c(cnv.N, cnvDEL.N, cnvDUP.N, cnv.rate, cnvDEL.rate, cnvDUP.rate, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, OR.LOAD, pvalue.LOAD, ORf.EOAD, pvaluef.EOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
      
    }
    res.withoutAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "none", t(res))
    colnames(res.withoutAdjustment) <- cngeneral
    
    ### WITH APOE ADJUSTMENT
    # cat(paste("APOE4", iter, "\n"))
    ### Data for the current transcript
    dat <- subset(genematrix[, c("sample", "Status", "apoe_genotype", transcript.list[iter])], apoe_genotype != "NA")
    dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
    dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
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
    
    if (0 %in% count["Sum",]) {
      res <- c(cnv.N, cnvDEL.N, cnvDUP.N, cnv.rate, cnvDEL.rate, cnvDUP.rate, rep(NA, 31))
    } else { 
      # cat(paste("APOE 4 AD", iter, "\n"))
      # OR AD
      dat.AD <- dat
      dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
      m <- glm(Y ~ cnv + APOE4, data = dat.AD, family = binomial(link = "logit"))
      OR.AD <- c(exp(m$coef)["cnv"], exp(confint(m))["cnv",])
      
      pvalue.AD <- summary(m)$coefficients["cnv", "Pr(>|z|)"]
      
      # cat(iter, " LOGISTF AD \n")
      mf <- logistf(Y ~ cnv + APOE4, data = dat.AD)
      ORf.AD <- c(exp(mf$coefficients["cnv"]), exp(confint(mf))["cnv",])
      pvaluef.AD <- mf$prob["cnv"]
      
      
      # OR EOAD
      dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
      if (sum(dat.EOAD$cnv != 2, na.rm = TRUE) == 0) {
        OR.EOAD <- c(NA, NA, NA)
        ORf.EOAD <- c(NA, NA, NA)
        pvalue.EOAD <- NA
        pvaluef.EOAD <- NA
      } else {
        dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
        m <- glm(Y ~ cnv + APOE4, data = dat.EOAD, family = binomial(link = "logit"))
        OR.EOAD <- c(exp(m$coef)["cnv"], exp(confint(m))["cnv",])
        pvalue.EOAD <- summary(m)$coefficients["cnv", "Pr(>|z|)"]
        
        # cat(iter, " LOGISTF EOAD APOE4 \n")
        mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
        ORf.EOAD <- c(exp(mf$coefficients["cnv"]), exp(confint(mf))["cnv",])
        pvaluef.EOAD <- mf$prob["cnv"]
      }
      
      # OR LOAD
      dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
      if (sum(dat.LOAD$cnv != 2, na.rm = TRUE) == 0) {
        OR.LOAD <- c(NA, NA, NA)
        ORf.LOAD <- c(NA, NA, NA)
        pvalue.LOAD <- NA
        pvaluef.LOAD <- NA
      } else {
        dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
        m <- glm(Y ~ cnv + APOE4, data = dat.LOAD, family = binomial(link = "logit"))
        OR.LOAD <- c(exp(m$coef)["cnv"], exp(confint(m))["cnv",])
        pvalue.LOAD <- summary(m)$coefficients["cnv", "Pr(>|z|)"]
        
        # cat(iter, " LOGISTF LOAD APOE4 \n")
        mf <- logistf(Y ~ cnv + APOE4, data = dat.LOAD)
        ORf.LOAD <- c(exp(mf$coefficients["cnv"]), exp(confint(mf))["cnv",])
        pvaluef.LOAD <- mf$prob["cnv"]
      }
      # cat(paste("APOE4 clm", iter, "\n"))
      # Ordinal regression avec clm
      rego <- clm(Status.f ~ cnv + APOE4, data = na.omit(dat))
      confint.rego <- try(exp(confint(rego))["cnv",], silent = TRUE)
      if (inherits(confint.rego, "try-error")) {
        OR.clm <- c(exp(rego$coefficients["cnv"]), NA, NA)
      } else {
        OR.clm <- c(exp(rego$coefficients["cnv"]), confint.rego)
      }
      
      pvalue.clm <- summary(rego)$coefficients["cnv", "Pr(>|z|)"]
      message.clm <- rego$message
      
      # Ordinal regression avec polr
      rego2 <- polr(Status.f ~ cnv + APOE4, data = na.omit(dat), Hess = TRUE)
      confint.rego2 <- try(exp(confint(rego2))["cnv",], silent = TRUE)
      if (inherits(confint.rego2, "try-error")) {
        OR.polr <- c(exp(rego2$coefficients["cnv"]), NA, NA)
      } else {
        OR.polr <- c(exp(rego2$coefficients["cnv"]), confint.rego2)
      }
      
      rego0 <- polr(Status.f ~ 1 + APOE4, data = na.omit(dat))
      statRV <- -2*(logLik(rego0) - logLik(rego2))
      pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
      
      res <- c(cnv.N, cnvDEL.N, cnvDUP.N, cnv.rate, cnvDEL.rate, cnvDUP.rate, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, OR.LOAD, pvalue.LOAD, ORf.EOAD, pvaluef.EOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
    }
    res.withAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "APOE nb4 (quanti)", t(res))
    colnames(res.withAdjustment) <- cngeneral 
    
    return(rbind(res.withoutAdjustment, res.withAdjustment))
    
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
  # write.table(temp, paste0("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  write.table(temp, paste0("RESULTS/Table_for_article_october2023_Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
}

