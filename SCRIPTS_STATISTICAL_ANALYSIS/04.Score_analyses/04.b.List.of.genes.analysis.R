setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

list_of_genes.analysis <- function(type.CNV, filtre.freq.value, filtre.freq.type.recouvrement, list_of_genes, p.value_threshold) {
  
  Matrices <- import.matrice(type = "CNV", freq = filtre.freq.value, type.recouvrement = filtre.freq.type.recouvrement)
  
  filtreCNV <- read.table(paste0("DATA/count_CNV_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
  transcript.list1 <- as.character(subset(filtreCNV, count.cnv >= 1 & exclusion == FALSE)$transcript.name)
  transcript.list <- transcript.list1[substr(transcript.list1, 1, 2) == "NM"]
  transcript.list.in.matrix <- as.character(subset(filtreCNV, exclusion == FALSE & substr(transcript.name, 1, 2) == "NM")$transcript.name)
  transcript.list.in.matrix.before.freq.filtre <- subset(filtreCNV, substr(transcript.name, 1, 2) == "NM")$transcript.name
  
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
  
  
  
  if (list_of_genes == "ABeta") {
    aBeta <- read.table("DATA/SCORES_AND_LIST_OF_GENES/List of genes Abeta network v4 full_14OCT2022.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    # create ensemble of lists of genes
    ensemble.list.Gene <- list()
    ensemble.list.Gene[[1]] <- aBeta$Gene
    # create names of lists of genes
    name.list.Gene <- c("A-Beta")
    # Increment the ensemble of lists of genes
    k <- 2
    for (iter in unique(aBeta$Category)) {
      ensemble.list.Gene[[k]] <- aBeta[aBeta$Category == iter, "Gene"]
      name.list.Gene[k] <- iter
      k <- k + 1
      for (iter2 in unique(aBeta[aBeta$Category == iter, "SubCategory"])) {
        ensemble.list.Gene[[k]] <- aBeta[aBeta$Category == iter & aBeta$SubCategory == iter2, "Gene"]
        name.list.Gene[k] <- paste(iter, " -> Sub: ", iter2)
        k <- k + 1
      }
    }
  }
  
  if (list_of_genes == "GWAS_EADB") {
    
    EADB <- read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    ensemble.list.Gene <- list()
    ensemble.list.Gene[[1]] <- unique(EADB$V1)
    name.list.Gene <- c("GWAS - EADB")
    
  }
  
  if (list_of_genes == "All") {
    
    ensemble.list.Gene <- list()
    ensemble.list.Gene[[1]] <- unique(GeneTranscript$Gene)
    name.list.Gene <- c("All genes in GeneTranscript file")
    
  }
  
  # ensemble.list.Gene <- aBeta$Gene[1:100]
  
  B <- length(ensemble.list.Gene)
  
  # cl <- parallel::makeCluster(30)
  # doParallel::registerDoParallel(cl)
  
  resultsB <- foreach(iterB = 1:B, .combine = "rbind") %do% {
  # for(iterB in 1:B) {
    cat(iterB, "/", B, "\n")
    cngeneral <- c("list.of.gene.name", "Adjustment", "N CTRL", "N LOAD", "N EOAD", "N total", 
                   "Ncarriers CTRL", "Ncarriers LOAD", "Ncarriers EOAD", "Ncarriers total", 
                   "Rate CTRL", "Rate LOAD", "Rate EOAD", "Rate total", 
                   "N genes in the list", "N genes with related transcript in the matrix", 
                   "N genes with related analysable transcript in the matrix", "N genes with related transcript in the matrix after filter for differential missingness", 
                   "N genes encompassing at least one CNV among analysable transcripts", "N analysable transcripts encompassing at least one CNV", 
                   "min of N genes affected per sample", "max of N genes affected per sample", "min of N genes being NA per sample", "max of N genes being NA per sample", 
                   "OR allAD vs CTRL", "OR allAD vs CTRL - CIinf", "OR allAD vs CTRL - CIsup", "p-value allAD vs CTRL",
                   "OR firth allAD vs CTRL", "OR firth allAD vs CTRL - CIinf", "OR firth allAD vs CTRL - CIsup", "p-value firth allAD vs CTRL",
                   "OR EOAD vs CTRL", "OR EOAD vs CTRL - CIinf", "OR EOAD vs CTRL - CIsup", "pvalue EOAD vs CTRL",
                   "OR firth EOAD vs CTRL", "OR firth EOAD vs CTRL - CIinf", "OR firth EOAD vs CTRL - CIsup", "pvalue firth EOAD vs CTRL",
                   "OR LOAD vs CTRL", "OR LOAD vs CTRL - CIinf", "OR LOAD vs CTRL - CIsup", "pvalue LOAD vs CTRL",
                   "OR firth LOAD vs CTRL", "OR firth LOAD vs CTRL - CIinf", "OR firth LOAD vs CTRL - CIsup", "pvalue firth LOAD vs CTRL",
                   "OR clm", "OR clm - CIinf", "OR clm - CIsup", "OR polr", "OR polr - CIinf", "OR polr - CIsup", "p-value (ordinal regression, using clm)", "p-value (ordinal regression, using polr)", "warning clm")
    
 
    
    library(ordinal)
    library(MASS)
    
    ### WITHOUT APOE ADJUSTMENT
    
    transcript.DMexcluded.ordinal.regression <- as.character(subset(differentialMissingness, Selection == "none" & (is.na(DM.for.EOAD.vs.CTRL.p.value..chi2.) | DM.for.EOAD.vs.CTRL.p.value..chi2. < p.value_threshold))$transcript.name)
    transcript.list.analysis <- transcript.list[!(transcript.list %in% transcript.DMexcluded.ordinal.regression)]
    GeneTranscript.list.analysis <- subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.analysis)
    list.geneSCORE <- unique(GeneTranscript.list.analysis$Gene)
    
    if (length(list.geneSCORE) == 0) { # if none gene correspond to transcript passing filters
      
      count.genes.original.matrix <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.in.matrix.before.freq.filtre)$Gene))
      count.genes.analysable.before.DM.exclusion <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.in.matrix)$Gene))
      transcript.list.afterDMexclusion <- transcript.list.in.matrix[!(transcript.list.in.matrix %in% transcript.DMexcluded.ordinal.regression)]
      count.genes.analysable.after.DMexclusion <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.afterDMexclusion)$Gene))
      count.genes.in.list.with1 <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list)$Gene))
      res1 <- c(rep(NA, 12),  length(ensemble.list.Gene[[iterB]]), count.genes.original.matrix, count.genes.analysable.before.DM.exclusion,
                count.genes.analysable.after.DMexclusion, count.genes.in.list.with1, rep(NA, 38))
      res <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (no gene passing transcript filter)"), Adjustment = "", t(res1), stringsAsFactors = FALSE)
      colnames(res) <- cngeneral
    } else {
      
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
      if (is.null(dim(summary_by_gene))) {
        summary_by_gene2 = data.frame(V1 = summary_by_gene)
        colnames(summary_by_gene2) <- list.geneSCORE
        summary_by_gene <- summary_by_gene2
      } else {
        colnames(summary_by_gene) <- list.geneSCORE
      }
      #write.table(cbind(genematrix[, c("sample", "Status", "apoe_genotype")], summary_by_gene), paste0("DATA/Summary_by_list_of_gene_", list_of_genes, "_", name.list.Gene[iter], "_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, "_differential_missingness_threshold_", p.value_threshold, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
      # compute the score
      if (is.null(dim(summary_by_gene))) {
        sum.listGene <- summary_by_gene[,1]
      } else {
        sum.listGene <- apply(summary_by_gene, 1, function(x) {sum(x, na.rm = TRUE)})
      } 
      dat <- cbind(genematrix[, c("sample", "Status")], sum.listGene)
      dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
      dat$status.listGene <- factor(ifelse(dat$sum.listGene > 0, 1, 0), levels = c(0, 1))
      
      ### Counts
      count <- addmargins(table(dat$status.listGene, dat$Status.f))
      cnv.N <- c(count["Sum",], count["1",])
      cnv.rate <- c(count["1",]/count["Sum",])
      
      ### Nb genes included in the analysis
      count.gene.transcript.included <- c(length(list.geneSCORE), length(unique(GeneTranscript.list.analysis$Transcript)))
      if (is.null(dim(summary_by_gene))) {
        count.gene.affected.by.individual <- summary_by_gene
        count.gene.NA.by.individual <- as.numeric(is.na(summary_by_gene))
      } else {
        count.gene.affected.by.individual <- apply(summary_by_gene, 1, function(x) {sum(x == 1, na.rm = TRUE)})
        count.gene.NA.by.individual <- apply(summary_by_gene, 1, function(x) {sum(is.na(x))})
      } 
      
      count.gene.affected.by.individual.min <- min(count.gene.affected.by.individual)
      count.gene.affected.by.individual.max <- max(count.gene.affected.by.individual) 
      
      count.gene.NA.by.individual.min <- min(count.gene.NA.by.individual)
      count.gene.NA.by.individual.max <- max(count.gene.NA.by.individual) 
      
      # nb genes analysable
      count.genes.original.matrix <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.in.matrix.before.freq.filtre)$Gene))
      count.genes.analysable.before.DM.exclusion <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.in.matrix)$Gene))
      transcript.list.afterDMexclusion <- transcript.list.in.matrix[!(transcript.list.in.matrix %in% transcript.DMexcluded.ordinal.regression)]
      count.genes.analysable.after.DMexclusion <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.afterDMexclusion)$Gene))
      
      count.ordinal.regression <- c(cnv.N, cnv.rate, length(ensemble.list.Gene[[iterB]]), count.genes.original.matrix, count.genes.analysable.before.DM.exclusion,
                                    count.genes.analysable.after.DMexclusion, count.gene.transcript.included, 
                                    count.gene.affected.by.individual.min, count.gene.affected.by.individual.max,
                                    count.gene.NA.by.individual.min, count.gene.NA.by.individual.max)
      
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
      
      res <- c(count.ordinal.regression, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
      
      
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
      if (sum(dat.EOAD$sum.listGene == 1, na.rm = TRUE) == 0) {
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
      if (sum(dat.LOAD$sum.listGene == 1, na.rm = TRUE) == 0) {
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
      
      res <- c(count.ordinal.regression, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
      
      
      res.sum.withoutAdjustment <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (sum)"), Adjustment = "none", t(res), stringsAsFactors = FALSE)
      colnames(res.sum.withoutAdjustment) <- cngeneral
      
      ### WITH APOE ADJUSTMENT
      
      genematrixAPOE4 <- subset(genematrix, apoe_genotype != "NA")
      
      ### Data for the current transcript
      
      transcript.DMexcluded.ordinal.regression <- as.character(subset(differentialMissingness, Selection == "APOE available" & is.na(DM.for.ordinal.regression.p.value..chi2.) | DM.for.ordinal.regression.p.value..chi2. < p.value_threshold)$transcript.name)
      transcript.list.analysis <- transcript.list[!(transcript.list %in% transcript.DMexcluded.ordinal.regression)]
      GeneTranscript.list.analysis <- subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.analysis)
      list.geneSCORE <- unique(GeneTranscript.list.analysis$Gene)
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
      if (is.null(dim(summary_by_gene))) {
        sum.listGene <- summary_by_gene
      } else {
        sum.listGene <- apply(summary_by_gene, 1, function(x) {sum(x, na.rm = TRUE)})
      } 
      dat <- cbind(genematrixAPOE4[, c("sample", "Status", "apoe_genotype")], sum.listGene)
      dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
      dat$status.listGene <- factor(ifelse(dat$sum.listGene > 0, 1, 0), levels = c(0, 1))
      dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
      
      ### Counts
      count <- addmargins(table(dat$status.listGene, dat$Status.f))
      cnv.N <- c(count["Sum",], count["1",])
      cnv.rate <- c(count["1",]/count["Sum",])
      
      ### Nb genes included in the analysis
      count.gene.transcript.included <- c(length(list.geneSCORE), length(unique(GeneTranscript.list.analysis$Transcript)))
      
      if (is.null(dim(summary_by_gene))) {
        count.gene.affected.by.individual <- summary_by_gene
        count.gene.NA.by.individual <- as.numeric(is.na(summary_by_gene))
      } else {
        count.gene.affected.by.individual <- apply(summary_by_gene, 1, function(x) {sum(x == 1, na.rm = TRUE)})
        count.gene.NA.by.individual <- apply(summary_by_gene, 1, function(x) {sum(is.na(x))})
      } 
      
      count.gene.affected.by.individual.min <- min(count.gene.affected.by.individual)
      count.gene.affected.by.individual.max <- max(count.gene.affected.by.individual) 
      
      count.gene.NA.by.individual.min <- min(count.gene.NA.by.individual)
      count.gene.NA.by.individual.max <- max(count.gene.NA.by.individual) 
      
      # nb genes analysable
      count.genes.original.matrix <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.in.matrix.before.freq.filtre)$Gene))
      count.genes.analysable.before.DM.exclusion <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.in.matrix)$Gene))
      transcript.list.afterDMexclusion <- transcript.list.in.matrix[!(transcript.list.in.matrix %in% transcript.DMexcluded.ordinal.regression)]
      count.genes.analysable.after.DMexclusion <- length(unique(subset(GeneTranscript, Gene %in% ensemble.list.Gene[[iterB]] & Transcript %in% transcript.list.afterDMexclusion)$Gene))
      
      count.ordinal.regression <- c(cnv.N, cnv.rate, length(ensemble.list.Gene[[iterB]]), count.genes.original.matrix, count.genes.analysable.before.DM.exclusion,
                                    count.genes.analysable.after.DMexclusion, count.gene.transcript.included, 
                                    count.gene.affected.by.individual.min, count.gene.affected.by.individual.max,
                                    count.gene.NA.by.individual.min, count.gene.NA.by.individual.max)
      
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
      
      res <- c(count.ordinal.regression, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
      
      
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
      if (sum(dat.EOAD$sum.listGene == 1, na.rm = TRUE) == 0) {
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
      if (sum(dat.LOAD$sum.listGene == 1, na.rm = TRUE) == 0) {
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
      
      res <- c(count.ordinal.regression, OR.AD, pvalue.AD, ORf.AD, pvaluef.AD, OR.EOAD, pvalue.EOAD, ORf.EOAD, pvaluef.EOAD, OR.LOAD, pvalue.LOAD, ORf.LOAD, pvaluef.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
      
      
      res.sum.withAdjustment <- data.frame(list.of.gene.name = paste0(name.list.Gene[[iterB]], " (sum)"), Adjustment = "APOE nb4 (quanti)", t(res), stringsAsFactors = FALSE)
      colnames(res.sum.withAdjustment) <- cngeneral
      
      res <- rbind(res.binary.withoutAdjustment, res.sum.withoutAdjustment, res.binary.withAdjustment, res.sum.withAdjustment)
      
    }
    return(res)
    # print(res)
  }
  
  
  # save results
  write.table(resultsB, paste0("RESULTS/Analysis_by_list_of_gene_", list_of_genes, "_ordinal_regression_and_subset_analyses_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, "_differential_missingness_threshold_", p.value_threshold, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
} 
