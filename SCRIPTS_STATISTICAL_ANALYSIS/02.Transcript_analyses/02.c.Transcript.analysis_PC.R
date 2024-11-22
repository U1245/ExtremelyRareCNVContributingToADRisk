setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")


transcript.analysis.PC <- function(type.CNV, filtre.freq.value, filtre.freq.type.recouvrement) {
  
  Matrices <- import.matrice(type = "CNV", freq = filtre.freq.value, type.recouvrement = filtre.freq.type.recouvrement)
  
  filtreCNV <- read.table(paste0("DATA/count_CNV_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
  transcript.list <- as.character(subset(filtreCNV, count.cnv >= 3 & exclusion == FALSE)$transcript.name)
  
  load("../clinData.Rdata")
  keptsample <- subset(clinData, !is.na(PC1))$sample
  
  if (type.CNV == "DEL") {
    genematrix1 <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DUP") {
    genematrix1 <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DEL_and_DUPpartial") {
    genematrix1 <- Matrices$DelMatriceWithPartial[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DUPpartial") {
    genematrix1 <- Matrices$DupPartialMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  # Attention pour DUP partial seules, ils faut recalculer filtre CNV
  
  genematrix2 <- subset(genematrix1, sample %in% keptsample)
  genematrix <- merge(subset(clinData, sample %in% genematrix2$sample)[, c("sample", paste0("PC", 1:10))], genematrix2, all.x = TRUE, all.y = FALSE, by = "sample")
  
  print(dim(genematrix1))
  print(dim(genematrix2))
  print(dim(genematrix))
  
  B <- length(transcript.list)
  # B <- 9
  cl <- parallel::makeCluster(29) 
  doParallel::registerDoParallel(cl)
  
  resultsB <- foreach(iter = 1:B, .combine = "rbind", .packages = c("logistf", "brglm2")) %dopar% {
    library(ordinal)
    library(MASS)
    cngeneral <- c("transcript.name", "Adjustment", "N CTRL", "N LOAD", "N EOAD", "N total", "Ncarriers CTRL", "Ncarriers LOAD", "Ncarriers EOAD", "Ncarriers total", "Rate CTRL", "Rate LOAD", "Rate EOAD", "Rate total", "OR all AD vs CTRL", "OR all AD vs CTRL - CIinf", "OR all AD vs CTRL - CIsup", "p-value all AD vs CTRL", "OR firth all AD vs CTRL", "OR firth all AD vs CTRL - CIinf", "OR firth all AD vs CTRL - CIsup", "p-value firth all AD vs CTRL",  "OR EOAD vs CTRL", "OR EOAD vs CTRL - CIinf", "OR EOAD vs CTRL - CIsup", "p-value EOAD vs CTRL", "OR LOAD vs CTRL", "OR LOAD vs CTRL - CIinf", "OR LOAD vs CTRL - CIsup", "p-value LOAD vs CTRL", "OR firth EOAD vs CTRL", "OR firth EOAD vs CTRL - CIinf", "OR firth EOAD vs CTRL - CIsup", "p-value firth EOAD vs CTRL", "OR firth LOAD vs CTRL", "OR firth LOAD vs CTRL - CIinf", "OR firth LOAD vs CTRL - CIsup", "p-value firth LOAD vs CTRL",  "OR clm", "OR clm - CIinf", "OR clm - CIsup", "OR polr", "OR polr - CIinf", "OR polr - CIsup", "p-value (ordinal regression, using clm)", "p-value (ordinal regression, using polr)", "warning clm", "OR bracl", "OR bracl - CIinf", "OR bracl - CIsup", "p-value (ordinal regression (adjacent model), using bracl)", "convergence bracl", "type bracl")
    
    ### WITHOUT APOE ADJUSTMENT
    # sink("log.txt", append=TRUE)  
    # cat(paste("Starting iteration", iter, "\n"))  
    ### Data for the current transcript
    dat <- genematrix[, c("sample", "Status", paste0("PC", 1:10), transcript.list[iter])]
    dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
    dat$cnv <- as.factor(dat[,transcript.list[iter]])
    
    ### Counts
    if (dim(table(dat$cnv, dat$Status.f))[1] == 0) {
      res <- rep(NA, 51)
      res.withoutAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "none", t(res), stringsAsFactors = FALSE)
      colnames(res.withoutAdjustment) <-cngeneral
      res <- rep(NA, 51)
      res.withAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "APOE nb4 (quanti)", t(res), stringsAsFactors = FALSE)
      colnames(res.withAdjustment) <- cngeneral
    } else {
      
      count <- addmargins(table(dat$cnv, dat$Status.f))
      
      if (sum(count["0",]==count["Sum",]) == 4) {
        res <- rep(NA, 51)
        res.withoutAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "none", t(res), stringsAsFactors = FALSE)
        colnames(res.withoutAdjustment) <-cngeneral
        res <- rep(NA, 51)
        res.withAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "APOE nb4 (quanti)", t(res), stringsAsFactors = FALSE)
        colnames(res.withAdjustment) <- cngeneral
  
        
      } else {
        cnv.N <- c(count["Sum",], count["1",])
        cnv.rate <- c(count["1",]/count["Sum",])
        
        
        if (0 %in% count["Sum",]) {
          res <- c(cnv.N, cnv.rate, rep(NA, 39))
        } else {
          
          # OR AD
          dat.AD <- dat
          dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
          m <- glm(Y ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.AD, family = binomial(link = "logit"))
          OR.AD <- c(exp(m$coef)["cnv1"], exp(confint(m))["cnv1",])
          
          pvalue.AD <- summary(m)$coefficients["cnv1", "Pr(>|z|)"]
          
          # cat(iter, " LOGISTF AD \n")
          mf <- logistf(Y ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.AD)
          ORf.AD <- c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
          pvaluef.AD <- mf$prob["cnv1"]
          
          
          # OR EOAD
          dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
          if (sum(dat.EOAD$cnv == 1, na.rm = TRUE) == 0) {
            OR.EOAD <- c(NA, NA, NA)
            ORf.EOAD <- c(NA, NA, NA)
            pvalue.EOAD <- NA
            pvaluef.EOAD <- NA
          } else {
            dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
            m <- glm(Y ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.EOAD, family = binomial(link = "logit"))
            OR.EOAD <- c(exp(m$coef)["cnv1"], exp(confint(m))["cnv1",])
            
            pvalue.EOAD <- summary(m)$coefficients["cnv1", "Pr(>|z|)"]
            
            # cat(iter, " LOGISTF EOAD \n")
            mf <- logistf(Y ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.EOAD)
            ORf.EOAD <- c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
            pvaluef.EOAD <- mf$prob["cnv1"]
          }
          
          # OR LOAD
          dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
          if (sum(dat.LOAD$cnv == 1, na.rm = TRUE) == 0) {
            OR.LOAD <- c(NA, NA, NA)
            ORf.LOAD <- c(NA, NA, NA)
            pvalue.LOAD <- NA
            pvaluef.LOAD <- NA
          } else {
            dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
            m <- glm(Y ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.LOAD, family = binomial(link = "logit"))
            OR.LOAD <- c(exp(m$coef)["cnv1"], exp(confint(m))["cnv1",])
            pvalue.LOAD <- summary(m)$coefficients["cnv1", "Pr(>|z|)"]
            
            # cat(iter, " LOGISTF LOAD \n")
            mf <- logistf(Y ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.LOAD)
            ORf.LOAD <- c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
            pvaluef.LOAD <- mf$prob["cnv1"]
          }
          
          # Ordinal regression avec clm
          rego <- clm(Status.f ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = na.omit(dat))
          confint.rego <- try(exp(confint(rego))["cnv1",], silent = TRUE)
          if (inherits(confint.rego, "try-error")) {
            OR.clm <- c(exp(rego$coefficients["cnv1"]), NA, NA)
          } else {
            OR.clm <- c(exp(rego$coefficients["cnv1"]), confint.rego)
          }
          
          pvalue.clm <- summary(rego)$coefficients["cnv1", "Pr(>|z|)"]
          message.clm <- rego$message
          
          # Ordinal regression avec polr
          rego2 <- polr(Status.f ~ cnv + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = na.omit(dat), Hess = TRUE)
          confint.rego2 <- try(exp(confint(rego2)["cnv1",]), silent = TRUE)
          if (inherits(confint.rego2, "try-error")) {
            OR.polr <- c(exp(rego2$coefficients["cnv1"]), NA, NA)
          } else {
            OR.polr <- c(exp(rego2$coefficients["cnv1"]), confint.rego2)
          }
          
          rego0 <- polr(Status.f ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = na.omit(dat))
          statRV <- -2*(logLik(rego0) - logLik(rego2))
          pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
          
          
          # # Ordinal regression avec bracl
          # cat(iter, " bracl \n")
          # regobr <- bracl(Status.f ~ cnv, data = na.omit(dat), parallel = TRUE, type = "MPL_Jeffreys")
          # confint.regobr <- try(exp(-confint(regobr))["cnv1",], silent = TRUE)
          # if (inherits(confint.regobr, "try-error")) {
          #   OR.bracl <- c(exp(-regobr$coefficients["cnv1"]), NA, NA)
          # } else {
          #   OR.bracl <- c(exp(-regobr$coefficients["cnv1"]), confint.regobr)
          # }
          # 
          # 
          # pvalue.bracl <- summary(regobr)$coefficients["cnv1", "Pr(>|z|)"]
          # message.bracl <- c(regobr$converged, regobr$type)
          OR.bracl <- c(NA, NA, NA); pvalue.bracl <- NA; message.bracl <- c(NA, NA)
          
          res <- c(cnv.N, cnv.rate, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD,  OR.EOAD, pvalue.EOAD, OR.LOAD, pvalue.LOAD, ORf.EOAD, pvaluef.EOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm, OR.bracl, pvalue.bracl, message.bracl)
          
        }
        
        res.withoutAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "none", t(res), stringsAsFactors = FALSE)
        colnames(res.withoutAdjustment) <- cngeneral
        
        ### WITH APOE ADJUSTMENT
        
        ### Data for the current transcript
        dat <- subset(genematrix[, c("sample", "Status", "apoe_genotype", paste0("PC", 1:10), transcript.list[iter])], apoe_genotype != "NA")
        dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
        dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
        dat$cnv <- as.factor(dat[,transcript.list[iter]])
        
        ### Counts
        count <- addmargins(table(dat$cnv, dat$Status.f))
        cnv.N <- c(count["Sum",], count["1",])
        cnv.rate <- c(count["1",]/count["Sum",])
        
        
        if (0 %in% count["Sum",]) {
          res <- c(cnv.N, cnv.rate, rep(NA, 39))
        } else {
          
          # OR AD
          dat.AD <- dat
          dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
          m <- glm(Y ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.AD, family = binomial(link = "logit"))
          OR.AD <- c(exp(m$coef)["cnv1"], exp(confint(m))["cnv1",])
          
          pvalue.AD <- summary(m)$coefficients["cnv1", "Pr(>|z|)"]
          
          # cat(iter, " LOGISTF AD \n")
          mf <- logistf(Y ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.AD)
          ORf.AD <- c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
          pvaluef.AD <- mf$prob["cnv1"]
          
          # OR EOAD
          dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
          if (sum(dat.EOAD$cnv == 1, na.rm = TRUE) == 0) {
            OR.EOAD <- c(NA, NA, NA)
            ORf.EOAD <- c(NA, NA, NA)
            pvalue.EOAD <- NA
            pvaluef.EOAD <- NA
          } else {
            dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
            m <- glm(Y ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.EOAD, family = binomial(link = "logit"))
            OR.EOAD <- c(exp(m$coef)["cnv1"], exp(confint(m))["cnv1",])
            pvalue.EOAD <- summary(m)$coefficients["cnv1", "Pr(>|z|)"]
            
            # cat(iter, " LOGISTF EOAD APOE4 \n")
            mf <- logistf(Y ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.EOAD)
            ORf.EOAD <- c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
            pvaluef.EOAD <- mf$prob["cnv1"]
          }
          
          # OR LOAD
          dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
          if (sum(dat.LOAD$cnv == 1, na.rm = TRUE) == 0) {
            OR.LOAD <- c(NA, NA, NA)
            ORf.LOAD <- c(NA, NA, NA)
            pvalue.LOAD <- NA
            pvaluef.LOAD <- NA
          } else {
            dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
            m <- glm(Y ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.LOAD, family = binomial(link = "logit"))
            OR.LOAD <- c(exp(m$coef)["cnv1"], exp(confint(m))["cnv1",])
            pvalue.LOAD <- summary(m)$coefficients["cnv1", "Pr(>|z|)"]
            
            # cat(iter, " LOGISTF LOAD APOE4 \n")
            mf <- logistf(Y ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = dat.LOAD)
            ORf.LOAD <- c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
            pvaluef.LOAD <- mf$prob["cnv1"]
          }
          
          # Ordinal regression avec clm
          rego <- clm(Status.f ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = na.omit(dat))
          confint.rego <- try(exp(confint(rego))["cnv1",], silent = TRUE)
          if (inherits(confint.rego, "try-error")) {
            OR.clm <- c(exp(rego$coefficients["cnv1"]), NA, NA)
          } else {
            OR.clm <- c(exp(rego$coefficients["cnv1"]), confint.rego)
          }
          
          pvalue.clm <- summary(rego)$coefficients["cnv1", "Pr(>|z|)"]
          message.clm <- rego$message
          
          # Ordinal regression avec polr
          rego2 <- polr(Status.f ~ cnv + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = na.omit(dat), Hess = TRUE)
          confint.rego2 <- try(exp(confint(rego2))["cnv1",], silent = TRUE)
          if (inherits(confint.rego2, "try-error")) {
            OR.polr <- c(exp(rego2$coefficients["cnv1"]), NA, NA)
          } else {
            OR.polr <- c(exp(rego2$coefficients["cnv1"]), confint.rego2)
          }
          
          
          rego0 <- polr(Status.f ~ 1 + APOE4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = na.omit(dat))
          statRV <- -2*(logLik(rego0) - logLik(rego2))
          pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
          
          
          # # Ordinal regression avec bracl
          # cat(iter, " BRACL APOE4 \n")
          # regobr <- bracl(Status.f ~ cnv + APOE4, data = na.omit(dat), parallel = TRUE)
          # confint.regobr <- try(exp(-confint(regobr))["cnv1",], silent = TRUE)
          # if (inherits(confint.regobr, "try-error")) {
          #   OR.bracl <- c(exp(-regobr$coefficients["cnv1"]), NA, NA)
          # } else {
          #   OR.bracl <- c(exp(-regobr$coefficients["cnv1"]), confint.regobr)
          # }
          # 
          # pvalue.bracl <- summary(regobr)$coefficients["cnv1", "Pr(>|z|)"]
          # message.bracl <- c(regobr$converged, regobr$type)
          OR.bracl <- c(NA, NA, NA); pvalue.bracl <- NA; message.bracl <- c(NA, NA)
          
          res <- c(cnv.N, cnv.rate, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, OR.LOAD, pvalue.LOAD, ORf.EOAD, pvaluef.EOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm, OR.bracl, pvalue.bracl, message.bracl)
          
        }
        res.withAdjustment <- data.frame(transcript.name = transcript.list[iter], Adjustment = "APOE nb4 (quanti)", t(res), stringsAsFactors = FALSE)
        colnames(res.withAdjustment) <- cngeneral
      }
    }
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
  write.table(temp, paste0("RESULTS/Analysis_by_transcript_PC_adjusted_ordinal_regression_and_subset_analyses_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
}

