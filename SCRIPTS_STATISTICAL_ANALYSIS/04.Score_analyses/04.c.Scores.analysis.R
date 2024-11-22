setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

# type.CNV = "DEL_and_DUPpartial"; filtre.freq.value = 0.01; filtre.freq.type.recouvrement = "reciproque"; score = "pLI"; p.value_threshold = 1e-5

score_genes.analysis <- function(type.CNV, filtre.freq.value, filtre.freq.type.recouvrement, score, p.value_threshold) {
  
  Matrices <- import.matrice(type = "CNV", freq = filtre.freq.value, type.recouvrement = filtre.freq.type.recouvrement)
  
  filtreCNV <- read.table(paste0("DATA/count_CNV_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
  transcript.list1 <- as.character(subset(filtreCNV, count.cnv >= 1 & exclusion == FALSE)$transcript.name)
  
  transcript.list <- transcript.list1[substr(transcript.list1, 1, 2) == "NM"]
 
  differentialMissingness <- read.table(paste0("DATA/Differential_missingness_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
  
  if (type.CNV == "DEL") {
    genematrix <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DUP") {
    genematrix <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DEL_and_DUPpartial") {
    genematrix <- Matrices$DelMatriceWithPartial[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  
  
  geneTranscript <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/GeneTranscript.Rdata")
  
  
  
  if (score == "pLI") {
    
    pLI1 <- read.table("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/Script_R_stat/Gene_annotation/scores/pLI/supplementary_dataset_11_full_constraint_metrics.tsv", sep = "\t", header = TRUE)
    pLI <- subset(pLI1, canonical == "true", select = c("gene", "pLI", "oe_lof_upper"))
    pLI$oe_lof_upper_modif <- 1 / pLI$oe_lof_upper
    
    # create ensemble of lists of genes
    ensemble.list.Gene <- list(pLI[, c("gene", "pLI")], pLI[, c("gene", "oe_lof_upper_modif")])
    name.list.Gene <- c("pLI", "oe_loff_upper_modif")
    
  }
  
  if (score == "cell") {
    # from "A cross-disorder dosage sensitivity map of the human genome", Cell 2022 accessible via BiblioInserm
    pLI1 <- read.table("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/STATISTICAL_ANALYSES/DATA/SCORES_AND_LIST_OF_GENES/1-s2.0-S0092867422007887-mmc7.csv", sep = "\t", dec = ",", header = TRUE)
    pLI <- subset(pLI1, select = c("Gene", "pHaplo", "pTriplo"))
    pLI$gene <- pLI$Gene
  
    # create ensemble of lists of genes
    ensemble.list.Gene <- list(pLI[, c("gene", "pHaplo")], pLI[, c("gene", "pTriplo")])
    name.list.Gene <- c("pHaplo", "pTriplo")
    
  }
  
  if (score == "DS") {
    DS <- read.table("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/Script_R_stat/Gene_annotation/scores/Hawrylycz et al 2015 DS score/NIHMS731485-supplement-4.csv", header = TRUE, sep = ",")
    DS$Pearson.modif <- DS$Pearson - min(DS$Pearson)
    # DS$Pearson.modif.log <- log(DS$Pearson.modif+1)
    DS$gene <- DS$Gene
    # create ensemble of lists of genes
    ensemble.list.Gene <- list(DS[, c("gene", "Pearson.modif")]) #, DS[, c("gene", "Pearson.modif.log")])
    name.list.Gene <- c("Pearson.modif") #, "Pearson.modif.log")
    
  }
  
  if (score == "GTEx") {
    GTEx <- read.table("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/Script_R_stat/Gene_annotation/scores/GTEx/GTEx_modif1.txt", sep = "\t", header = TRUE)
    
    GTEx$gene <- GTEx$Description
    
    ensemble.list.Gene <- list(GTEx[, c("gene", "Brain...Amygdala")], GTEx[, c("gene", "Brain...Anterior.cingulate.cortex..BA24.")], GTEx[, c("gene", "Brain...Caudate..basal.ganglia.")],
                               GTEx[, c("gene", "Brain...Cerebellar.Hemisphere")], GTEx[, c("gene", "Brain...Cerebellum")], GTEx[, c("gene", "Brain...Cortex")], GTEx[, c("gene", "Brain...Frontal.Cortex..BA9.")], 
                               GTEx[, c("gene", "Brain...Hippocampus")], GTEx[, c("gene", "Brain...Hypothalamus")], GTEx[, c("gene", "Brain...Nucleus.accumbens..basal.ganglia.")], 
                               GTEx[, c("gene", "Brain...Putamen..basal.ganglia.")], GTEx[, c("gene", "Brain...Spinal.cord..cervical.c.1.")], GTEx[, c("gene", "Brain...Substantia.nigra")],
                               GTEx[, c("gene", "Whole.Blood")])
    name.list.Gene <- c("Brain...Amygdala", "Brain...Anterior.cingulate.cortex..BA24.", "Brain...Caudate..basal.ganglia.", 
                        "Brain...Cerebellar.Hemisphere",  "Brain...Cerebellum", "Brain...Cortex", "Brain...Frontal.Cortex..BA9.", "Brain...Hippocampus", "Brain...Hypothalamus",
                        "Brain...Nucleus.accumbens..basal.ganglia.", "Brain...Putamen..basal.ganglia.", "Brain...Spinal.cord..cervical.c.1.", "Brain...Substantia.nigra",
                        "Whole.Blood")
    
  }
  
  
  
  B <- length(ensemble.list.Gene)
  
  # cl <- parallel::makeCluster(30)
  # doParallel::registerDoParallel(cl)
  
  resultsB <- foreach(iterB = 1:B, .combine = "rbind") %do% {
    
    cngeneral <- c("list.of.gene.name", "Adjustment", "N CTRL", "N LOAD", "N EOAD", "N total", 
                   "Ncarriers CTRL", "Ncarriers LOAD", "Ncarriers EOAD", "Ncarriers total", 
                   "median CTRL", "Q1 CTRL", "Q3 CTRL", "median LOAD", "Q1 LOAD", "Q3 LOAD", "median EOAD", "Q1 EOAD", "Q3 EOAD", 
                   "OR allAD vs CTRL", "OR allAD vs CTRL - CIinf", "OR allAD vs CTRL - CIsup", "pvalue allAD vs CTRL",
                   "OR firth allAD vs CTRL", "OR firth allAD vs CTRL - CIinf", "OR firth allAD vs CTRL - CIsup", "pvalue firth allAD vs CTRL",
                   "OR EOAD vs CTRL", "OR EOAD vs CTRL - CIinf", "OR EOAD vs CTRL - CIsup", "pvalue EOAD vs CTRL",
                   "OR firth EOAD vs CTRL", "OR firth EOAD vs CTRL - CIinf", "OR firth EOAD vs CTRL - CIsup", "pvalue firth EOAD vs CTRL",
                   "OR LOAD vs CTRL", "OR LOAD vs CTRL - CIinf", "OR LOAD vs CTRL - CIsup", "pvalue LOAD vs CTRL",
                   "OR firth LOAD vs CTRL", "OR firth LOAD vs CTRL - CIinf", "OR firth LOAD vs CTRL - CIsup", "pvalue firth LOAD vs CTRL",
                   "OR clm", "OR clm - CIinf", "OR clm - CIsup", "OR polr", "OR polr - CIinf", "OR polr - CIsup", "p-value (ordinal regression, using clm)", "p-value (ordinal regression, using polr)", "warning clm")
    
    cat(iterB, "/", B, "\n")
    
    library(ordinal)
    library(MASS)
    
    ### WITHOUT APOE ADJUSTMENT

    transcript.DMexcluded.ordinal.regression <- as.character(subset(differentialMissingness, Selection == "none" & (is.na(DM.for.EOAD.vs.CTRL.p.value..chi2.) | DM.for.EOAD.vs.CTRL.p.value..chi2. < p.value_threshold))$transcript.name)
    transcript.list.analysis <- transcript.list[!(transcript.list %in% transcript.DMexcluded.ordinal.regression)]
    GeneTranscript.list.analysis <- subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]]$gene & Transcript %in% transcript.list.analysis)
    list.geneSCORE <- sort(unique(GeneTranscript.list.analysis$Gene))
    
    # if several transcripts are related to the same gene, we summarize the information at gene level
    summary_by_gene <- foreach(iter = 1:length(list.geneSCORE), .combine = "cbind") %do% {
      dat_current_gene <- genematrix[, subset(GeneTranscript.list.analysis, Gene == list.geneSCORE[iter])$Transcript]
      if ( is.null(dim(dat_current_gene)) ) {
        return(dat_current_gene)
      } else {
        dat_current_gene[, list.geneSCORE[iter]] <- apply(dat_current_gene, 1, function(x) { as.numeric(sum(x) > 1)})
        return(dat_current_gene[, list.geneSCORE[iter]])
      }
    }
    
    # compute the score
    current_score_dat <- ensemble.list.Gene[[iterB]][ensemble.list.Gene[[iterB]]$gene %in% list.geneSCORE, ]
    current_score_dat.noduplica <- current_score_dat[!duplicated(current_score_dat$gene),]
    current_score <- current_score_dat.noduplica[order(current_score_dat.noduplica$gene), ] 
    
    sum.listGene <- apply(summary_by_gene, 1, function(x) {sum(x*current_score[, 2], na.rm = TRUE)})
    dat <- cbind(genematrix[, c("sample", "Status")], sum.listGene)
    dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
    dat$status.listGene <- factor(ifelse(dat$sum.listGene > 0, 1, 0), levels = c(0, 1))
    
    
    ### Counts
    count <- addmargins(table(dat$status.listGene, dat$Status.f))
    cnv.N <- c(count["Sum",], count["1",])
    
    cnv.median.q1.q3 <- c(summary(subset(dat, Status == "Control")$sum.listGene)[c(3, 2, 5)], summary(subset(dat, Status == "LOAD")$sum.listGene)[c(3, 2, 5)], summary(subset(dat, Status == "EOAD")$sum.listGene)[c(3, 2, 5)])
    
    ### Analysis of the score as binary predictor 
    
    # OR AD
    dat.AD <- dat
    dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
    m <- glm(Y ~ status.listGene, data = dat.AD, family = binomial(link = "logit"))
    OR.AD <- c(exp(m$coef)["status.listGene1"], exp(confint(m))["status.listGene1",])
    
    pvalue.AD <- summary(m)$coefficients["status.listGene1", "Pr(>|z|)"]
    
    # cat(iter, " LOGISTF AD \n")
    mf <- logistf(Y ~ status.listGene, data = dat.AD)
    ORf.AD <- c(exp(mf$coefficients["status.listGene1"]), exp(confint(mf))["status.listGene1",])
    pvaluef.AD <- mf$prob["status.listGene1"]
    
    
    
    # OR EOAD
    dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
    if (sum(dat.EOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.EOAD <- c(NA, NA, NA)
      ORf.EOAD <- c(NA, NA, NA)
      pvalue.EOAD <- NA
      pvaluef.EOAD <- NA
    } else {
      dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
      m <- glm(Y ~ status.listGene, data = dat.EOAD, family = binomial(link = "logit"))
      OR.EOAD <- c(exp(m$coef)["status.listGene1"], exp(confint(m))["status.listGene1",])
      pvalue.EOAD <- summary(m)$coefficients["status.listGene1", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ status.listGene, data = dat.EOAD)
      ORf.EOAD <- c(exp(mf$coefficients["status.listGene1"]), exp(confint(mf))["status.listGene1",])
      pvaluef.EOAD <- mf$prob["status.listGene1"]
    }
    
    # OR LOAD
    dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
    if (sum(dat.LOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.LOAD <- c(NA, NA, NA)
      ORf.LOAD <- c(NA, NA, NA)
      pvalue.LOAD <- NA
      pvaluef.LOAD <- NA
    } else {
      dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
      m <- glm(Y ~ status.listGene, data = dat.LOAD, family = binomial(link = "logit"))
      OR.LOAD <- c(exp(m$coef)["status.listGene1"], exp(confint(m))["status.listGene1",])
      pvalue.LOAD <- summary(m)$coefficients["status.listGene1", "Pr(>|z|)"]
      mf <- logistf(Y ~ status.listGene, data = dat.LOAD)
      ORf.LOAD <- c(exp(mf$coefficients["status.listGene1"]), exp(confint(mf))["status.listGene1",])
      pvaluef.LOAD <- mf$prob["status.listGene1"]
    }
    
    # Ordinal regression with clm
    rego <- clm(Status.f ~ status.listGene, data = na.omit(dat))
    confint.rego <- try(exp(confint(rego))["status.listGene1",], silent = TRUE)
    if (inherits(confint.rego, "try-error")) {
      OR.clm <- c(exp(rego$coefficients["status.listGene1"]), NA, NA)
    } else {
      OR.clm <- c(exp(rego$coefficients["status.listGene1"]), confint.rego)
    }
    pvalue.clm <- summary(rego)$coefficients["status.listGene1", "Pr(>|z|)"]
    message.clm <- rego$message
    
    # Ordinal regression with polr
    rego2 <- polr(Status.f ~ status.listGene, data = na.omit(dat), Hess = TRUE)
    confint.rego2 <- try(exp(confint(rego2)), silent = TRUE)
    if (inherits(confint.rego2, "try-error")) {
      OR.polr <- c(exp(rego2$coefficients["status.listGene1"]), NA, NA)
    } else {
      OR.polr <- c(exp(rego2$coefficients["status.listGene1"]), confint.rego2)
    }
    rego0 <- polr(Status.f ~ 1, data = na.omit(dat))
    statRV <- -2*(logLik(rego0) - logLik(rego2))
    pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
    
    res <- c(cnv.N, cnv.median.q1.q3, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
    
    res.binary.withoutAdjustment <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (yes/no)"), Adjustment = "none", t(res), stringsAsFactors = FALSE)
    colnames(res.binary.withoutAdjustment) <- cngeneral
    
    ### Analysis of the score as continuous predictor 
    
    # OR AD
    dat.AD <- dat
    dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
    m <- glm(Y ~ sum.listGene, data = dat.AD, family = binomial(link = "logit"))
    OR.AD <- c(exp(m$coef)["sum.listGene"], exp(confint(m))["sum.listGene",])
    
    pvalue.AD <- summary(m)$coefficients["sum.listGene", "Pr(>|z|)"]
    
    # cat(iter, " LOGISTF AD \n")
    mf <- logistf(Y ~ sum.listGene, data = dat.AD)
    ORf.AD <- c(exp(mf$coefficients["sum.listGene"]), exp(confint(mf))["sum.listGene",])
    pvaluef.AD <- mf$prob["sum.listGene"]
    
    
    # OR EOAD
    dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
    if (sum(dat.EOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.EOAD <- c(NA, NA, NA)
      ORf.EOAD <- c(NA, NA, NA)
      pvalue.EOAD <- NA
      pvaluef.EOAD <- NA
    } else {
      dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
      m <- glm(Y ~ sum.listGene, data = dat.EOAD, family = binomial(link = "logit"))
      OR.EOAD <- c(exp(m$coef)["sum.listGene"], exp(confint(m))["sum.listGene",])
      pvalue.EOAD <- summary(m)$coefficients["sum.listGene", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ sum.listGene, data = dat.EOAD)
      ORf.EOAD <- c(exp(mf$coefficients["sum.listGene"]), exp(confint(mf))["sum.listGene",])
      pvaluef.EOAD <- mf$prob["sum.listGene"]
    }
    
    # OR LOAD
    dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
    if (sum(dat.LOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.LOAD <- c(NA, NA, NA)
      ORf.LOAD <- c(NA, NA, NA)
      pvalue.LOAD <- NA
      pvaluef.LOAD <- NA
    } else {
      dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
      m <- glm(Y ~ sum.listGene, data = dat.LOAD, family = binomial(link = "logit"))
      OR.LOAD <- c(exp(m$coef)["sum.listGene"], exp(confint(m))["sum.listGene",])
      pvalue.LOAD <- summary(m)$coefficients["sum.listGene", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ sum.listGene, data = dat.LOAD)
      ORf.LOAD <- c(exp(mf$coefficients["sum.listGene"]), exp(confint(mf))["sum.listGene",])
      pvaluef.LOAD <- mf$prob["sum.listGene"]
    }

    # Ordinal regression with clm
    rego <- clm(Status.f ~ sum.listGene, data = na.omit(dat))
    confint.rego <- try(exp(confint(rego))["sum.listGene",], silent = TRUE)
    if (inherits(confint.rego, "try-error")) {
      OR.clm <- c(exp(rego$coefficients["sum.listGene"]), NA, NA)
    } else {
      OR.clm <- c(exp(rego$coefficients["sum.listGene"]), confint.rego)
    }
    pvalue.clm <- summary(rego)$coefficients["sum.listGene", "Pr(>|z|)"]
    message.clm <- rego$message
    
    # Ordinal regression with polr
    rego2 <- polr(Status.f ~ sum.listGene, data = na.omit(dat), Hess = TRUE)
    confint.rego2 <- try(exp(confint(rego2)), silent = TRUE)
    if (inherits(confint.rego2, "try-error")) {
      OR.polr <- c(exp(rego2$coefficients["sum.listGene"]), NA, NA)
    } else {
      OR.polr <- c(exp(rego2$coefficients["sum.listGene"]), confint.rego2)
    }
    rego0 <- polr(Status.f ~ 1, data = na.omit(dat))
    statRV <- -2*(logLik(rego0) - logLik(rego2))
    pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
    
    res <- c(cnv.N, cnv.median.q1.q3, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
    

    res.sum.withoutAdjustment <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (sum)"), Adjustment = "none", t(res), stringsAsFactors = FALSE)
    colnames(res.sum.withoutAdjustment) <- cngeneral
    
    ### WITH APOE ADJUSTMENT

    genematrixAPOE4 <- subset(genematrix, apoe_genotype != "NA")
    
    ### Data for the current transcript
    
    transcript.DMexcluded.ordinal.regression <- as.character(subset(differentialMissingness, Selection == "APOE available" & (is.na(DM.for.ordinal.regression.p.value..chi2.) | DM.for.ordinal.regression.p.value..chi2. < p.value_threshold))$transcript.name)
    transcript.list.analysis <- transcript.list[!(transcript.list %in% transcript.DMexcluded.ordinal.regression)]
    GeneTranscript.list.analysis <- subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]]$gene & Transcript %in% transcript.list.analysis)
    list.geneSCORE <- sort(unique(GeneTranscript.list.analysis$Gene))
    
    # if several transcripts are related to the same gene, we summarize the information at gene level
    summary_by_gene <- foreach(iter = 1:length(list.geneSCORE), .combine = "cbind") %do% {
      dat_current_gene <- genematrixAPOE4[, subset(GeneTranscript.list.analysis, Gene == list.geneSCORE[iter])$Transcript]
      if ( is.null(dim(dat_current_gene)) ) {
        return(dat_current_gene)
      } else {
        dat_current_gene[, list.geneSCORE[iter]] <- apply(dat_current_gene, 1, function(x) { as.numeric(sum(x) > 1)})
        return(dat_current_gene[, list.geneSCORE[iter]])
      }
    }
    
    # compute the score
    current_score_dat <- ensemble.list.Gene[[iterB]][ensemble.list.Gene[[iterB]]$gene %in% list.geneSCORE, ]
    current_score_dat.noduplica <- current_score_dat[!duplicated(current_score_dat$gene),]
    current_score <- current_score_dat.noduplica[order(current_score_dat.noduplica$gene), ] 
    
    sum.listGene <- apply(summary_by_gene, 1, function(x) {sum(x*current_score[, 2], na.rm = TRUE)})
    
    dat <- cbind(genematrixAPOE4[, c("sample", "Status", "apoe_genotype")], sum.listGene)
    dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
    dat$status.listGene <- factor(ifelse(dat$sum.listGene > 0, 1, 0), levels = c(0, 1))
    
    
    dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
    
    ### Counts
    count <- addmargins(table(dat$status.listGene, dat$Status.f))
    cnv.N <- c(count["Sum",], count["1",])
    
    cnv.median.q1.q3 <- c(summary(subset(dat, Status == "Control")$sum.listGene)[c(3, 2, 5)], summary(subset(dat, Status == "LOAD")$sum.listGene)[c(3, 2, 5)], summary(subset(dat, Status == "EOAD")$sum.listGene)[c(3, 2, 5)])
    
    
    ### Analysis of the score as binary predictor 
    
    # OR AD
    dat.AD <- dat
    dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
    m <- glm(Y ~ status.listGene + APOE4, data = dat.AD, family = binomial(link = "logit"))
    OR.AD <- c(exp(m$coef)["status.listGene1"], exp(confint(m))["status.listGene1",])
    
    pvalue.AD <- summary(m)$coefficients["status.listGene1", "Pr(>|z|)"]
    
    # cat(iter, " LOGISTF AD \n")
    mf <- logistf(Y ~ status.listGene + APOE4, data = dat.AD)
    ORf.AD <- c(exp(mf$coefficients["status.listGene1"]), exp(confint(mf))["status.listGene1",])
    pvaluef.AD <- mf$prob["status.listGene1"]
    
    # OR EOAD
    dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
    if (sum(dat.EOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.EOAD <- c(NA, NA, NA)
      ORf.EOAD <- c(NA, NA, NA)
      pvalue.EOAD <- NA
      pvaluef.EOAD <- NA
    } else {
      dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
      m <- glm(Y ~ status.listGene + APOE4, data = dat.EOAD, family = binomial(link = "logit"))
      OR.EOAD <- c(exp(m$coef)["status.listGene1"], exp(confint(m))["status.listGene1",])
      pvalue.EOAD <- summary(m)$coefficients["status.listGene1", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ status.listGene + APOE4, data = dat.EOAD)
      ORf.EOAD <- c(exp(mf$coefficients["status.listGene1"]), exp(confint(mf))["status.listGene1",])
      pvaluef.EOAD <- mf$prob["status.listGene1"]
    }
    
    # OR LOAD
    dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
    if (sum(dat.LOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.LOAD <- c(NA, NA, NA)
      ORf.LOAD <- c(NA, NA, NA)
      pvalue.LOAD <- NA
      pvaluef.LOAD <- NA
    } else {
      dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
      m <- glm(Y ~ status.listGene + APOE4, data = dat.LOAD, family = binomial(link = "logit"))
      OR.LOAD <- c(exp(m$coef)["status.listGene1"], exp(confint(m))["status.listGene1",])
      pvalue.LOAD <- summary(m)$coefficients["status.listGene1", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ status.listGene + APOE4, data = dat.LOAD)
      ORf.LOAD <- c(exp(mf$coefficients["status.listGene1"]), exp(confint(mf))["status.listGene1",])
      pvaluef.LOAD <- mf$prob["status.listGene1"]
    }
    
    
    
    # Ordinal regression with clm
    rego <- clm(Status.f ~ status.listGene + APOE4, data = na.omit(dat))
    confint.rego <- try(exp(confint(rego))["status.listGene1",], silent = TRUE)
    if (inherits(confint.rego, "try-error")) {
      OR.clm <- c(exp(rego$coefficients["status.listGene1"]), NA, NA)
    } else {
      OR.clm <- c(exp(rego$coefficients["status.listGene1"]), confint.rego)
    }
    pvalue.clm <- summary(rego)$coefficients["status.listGene1", "Pr(>|z|)"]
    message.clm <- rego$message
    
    
    # Ordinal regression with polr
    rego2 <- polr(Status.f ~ status.listGene + APOE4, data = na.omit(dat), Hess = TRUE)
    confint.rego2 <- try(exp(confint(rego2))["status.listGene1",], silent = TRUE)
    if (inherits(confint.rego2, "try-error")) {
      OR.polr <- c(exp(rego2$coefficients["status.listGene1"]), NA, NA)
    } else {
      OR.polr <- c(exp(rego2$coefficients["status.listGene1"]), confint.rego2)
    }
    rego0 <- polr(Status.f ~ 1 + APOE4, data = na.omit(dat))
    statRV <- -2*(logLik(rego0) - logLik(rego2))
    pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
    
    res <- c(cnv.N, cnv.median.q1.q3, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
    
    res.binary.withAdjustment <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (yes/no)"), Adjustment = "APOE nb4 (quanti)", t(res), stringsAsFactors = FALSE)
    colnames(res.binary.withAdjustment) <- cngeneral
    
    ### Analysis of the score as continuous predictor 
    
    # OR AD
    dat.AD <- dat
    dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
    m <- glm(Y ~ sum.listGene + APOE4, data = dat.AD, family = binomial(link = "logit"))
    OR.AD <- c(exp(m$coef)["sum.listGene"], exp(confint(m))["sum.listGene",])
    
    pvalue.AD <- summary(m)$coefficients["sum.listGene", "Pr(>|z|)"]
    
    # cat(iter, " LOGISTF AD \n")
    mf <- logistf(Y ~ sum.listGene + APOE4, data = dat.AD)
    ORf.AD <- c(exp(mf$coefficients["sum.listGene"]), exp(confint(mf))["sum.listGene",])
    pvaluef.AD <- mf$prob["sum.listGene"]
    
    
    # OR EOAD
    dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
    if (sum(dat.EOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.EOAD <- c(NA, NA, NA)
      ORf.EOAD <- c(NA, NA, NA)
      pvalue.EOAD <- NA
      pvaluef.EOAD <- NA
    } else {
      dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
      m <- glm(Y ~ sum.listGene + APOE4, data = dat.EOAD, family = binomial(link = "logit"))
      OR.EOAD <- c(exp(m$coef)["sum.listGene"], exp(confint(m))["sum.listGene",])
      pvalue.EOAD <- summary(m)$coefficients["sum.listGene", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ sum.listGene + APOE4, data = dat.EOAD)
      ORf.EOAD <- c(exp(mf$coefficients["sum.listGene"]), exp(confint(mf))["sum.listGene",])
      pvaluef.EOAD <- mf$prob["sum.listGene"]
    }
    
    # OR LOAD
    dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
    if (sum(dat.LOAD$status.listGene == 1, na.rm = TRUE) == 0) {
      OR.LOAD <- c(NA, NA, NA)
      ORf.LOAD <- c(NA, NA, NA)
      pvalue.LOAD <- NA
      pvaluef.LOAD <- NA
    } else {
      dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
      m <- glm(Y ~ sum.listGene + APOE4, data = dat.LOAD, family = binomial(link = "logit"))
      OR.LOAD <- c(exp(m$coef)["sum.listGene"], exp(confint(m))["sum.listGene",])
      pvalue.LOAD <- summary(m)$coefficients["sum.listGene", "Pr(>|z|)"]
      
      mf <- logistf(Y ~ sum.listGene + APOE4, data = dat.LOAD)
      ORf.LOAD <- c(exp(mf$coefficients["sum.listGene"]), exp(confint(mf))["sum.listGene",])
      pvaluef.LOAD <- mf$prob["sum.listGene"]
    }
    
    # Ordinal regression with clm
    rego <- clm(Status.f ~ sum.listGene + APOE4, data = na.omit(dat))
    confint.rego <- try(exp(confint(rego))["sum.listGene",], silent = TRUE)
    if (inherits(confint.rego, "try-error")) {
      OR.clm <- c(exp(rego$coefficients["sum.listGene"]), NA, NA)
    } else {
      OR.clm <- c(exp(rego$coefficients["sum.listGene"]), confint.rego)
    }
    pvalue.clm <- summary(rego)$coefficients["sum.listGene", "Pr(>|z|)"]
    message.clm <- rego$message
    
    # Ordinal regression with polr
    rego2 <- polr(Status.f ~ sum.listGene + APOE4, data = na.omit(dat), Hess = TRUE)
    confint.rego2 <- try(exp(confint(rego2))["sum.listGene",], silent = TRUE)
    if (inherits(confint.rego2, "try-error")) {
      OR.polr <- c(exp(rego2$coefficients["sum.listGene"]), NA, NA)
    } else {
      OR.polr <- c(exp(rego2$coefficients["sum.listGene"]), confint.rego2)
    }
    rego0 <- polr(Status.f ~ 1 + APOE4, data = na.omit(dat))
    statRV <- -2*(logLik(rego0) - logLik(rego2))
    pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
    
    res <- c(cnv.N, cnv.median.q1.q3, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
    
    
    res.sum.withAdjustment <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (sum)"), Adjustment = "APOE nb4 (quanti)", t(res), stringsAsFactors = FALSE)
    colnames(res.sum.withAdjustment) <- cngeneral
    
    
    res <- rbind(res.binary.withoutAdjustment, res.sum.withoutAdjustment, res.binary.withAdjustment, res.sum.withAdjustment)

    
    pdf(paste0("TEMP/TEST_Boxplot_for_Analysis_by_score_of_gene_", score, "_", name.list.Gene[iterB], "_ordinal_regression_and_subset_analyses_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, "_differential_missingness_threshold_", p.value_threshold, ".pdf"))
   
    twoQ3 <- 2*summary(dat$sum.listGene)[5]
    p <- ggplot(dat, aes(y = sum.listGene, fill = as.factor(APOE4), color = as.factor(APOE4))) + geom_boxplot()
    print(p)
    p <- ggplot(dat, aes(y = sum.listGene, fill = as.factor(APOE4))) + geom_boxplot() + coord_cartesian(ylim = c(0, twoQ3))
    print(p)
    
    p <- ggplot(dat, aes(y = sum.listGene, x = Status, fill = as.factor(APOE4), color = as.factor(APOE4))) + geom_boxplot()
    print(p)
    p <- ggplot(dat, aes(y = sum.listGene, x = Status, fill = as.factor(APOE4))) + geom_boxplot() + coord_cartesian(ylim = c(0, twoQ3))
    print(p)

    p <- ggplot(dat, aes(y = sum.listGene, x = Status)) + geom_boxplot()
    print(p)
    p <- ggplot(dat, aes(y = sum.listGene, x = Status)) + geom_boxplot() + coord_cartesian(ylim = c(0, twoQ3))
    print(p)
    
    dev.off()
    
    return(res)
    
    
  }

  # save results
  write.table(resultsB, paste0("RESULTS/TEST_Analysis_by_score_of_gene_", score, "_ordinal_regression_and_subset_analyses_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, "_differential_missingness_threshold_", p.value_threshold, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

  
  
}
  