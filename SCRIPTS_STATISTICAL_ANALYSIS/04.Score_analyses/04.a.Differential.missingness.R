setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

differential.missingness <- function(type.CNV, filtre.freq.value, filtre.freq.type.recouvrement, nbNonGenes. = nbNonGenes) {
  
  
  Matrices <- import.matrice(type = "CNV", freq = filtre.freq.value, type.recouvrement = filtre.freq.type.recouvrement)
  
  filtreCNV <- read.table(paste0("DATA/count_CNV_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
  transcript.list <- as.character(subset(filtreCNV, count.cnv >= 1 & exclusion == FALSE)$transcript.name)
  
  if (type.CNV == "DEL") {
    genematrix1 <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DUP") {
    genematrix1 <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }
  if (type.CNV == "DEL_and_DUPpartial") {
    genematrix1 <- Matrices$DelMatriceWithPartial[, c("sample", "Status", "apoe_genotype", transcript.list)]
  }


  sum.na <- apply(genematrix1[, (nbNonGenes. + 1):(dim(genematrix1)[2])], 2, function(x){sum(is.na(x))})
  # transcript.withoutNA <- names(sum.na[sum.na == 0])
  transcript.withNA <- names(sum.na[sum.na != 0])
  genematrix <- genematrix1[ , c("sample", "Status", "apoe_genotype", transcript.withNA)]
  
  B <- length(transcript.withNA)
  # cl <- parallel::makeCluster(20)
  # doParallel::registerDoParallel(cl)

  
  # relative missingness
  resultsB <- foreach(iter = 1:B, .combine = "rbind") %do% {
    cat(iter, "/", B, "/n")
    ### WHEN ANALYSES WERE MADE WITHOUT ADJUSTMENT
    
    
    # count available individuals
    count <- tapply(genematrix[, transcript.withNA[iter]], genematrix$Status, length)
    count.total <- c(count["Control"], count["LOAD"], count["EOAD"], sum(count))
    
    # count individuals with NA
    sum.na <- tapply(genematrix[, transcript.withNA[iter]], genematrix$Status, function(x){sum(is.na(x))})
    count.na <- c(sum.na["Control"], sum.na["LOAD"], sum.na["EOAD"], sum(sum.na))
    
    # rate individuals with NA
    rate.na <- count.na/count.total
    
    # count individuals with 1
    sum.1 <- tapply(genematrix[, transcript.withNA[iter]], genematrix$Status, function(x){sum((!is.na(x) & x==1))})
    count.1 <- c(sum.1["Control"], sum.1["LOAD"], sum.1["EOAD"], sum(sum.1))
    
    
    # matrice for chi2 and fisher exact test
    m <- matrix(c(sum.na[c("Control","LOAD", "EOAD")], count[c("Control","LOAD", "EOAD")]), byrow = FALSE, ncol = 2)
    
    # Test of missingness in case of ordinal regression
    p.chisq <- chisq.test(m)$p.value
    valid.chisq <- as.numeric(sum(chisq.test(m)$expected < 5) == 0)
    p.fisher <- fisher.test(m, simulate.p.value = TRUE, B = 50000)$p.value
    
    # Test of missingness in case of EOAD vs CTRL
    m.EOAD <- m[c(1,3),]
    p.chisq.EOADvsCTRL <- chisq.test(m.EOAD)$p.value
    valid.chisq.EOADvsCTRL <- as.numeric(sum(chisq.test(m.EOAD)$expected < 5) == 0)
    p.fisher.EOADvsCTRL <- fisher.test(m.EOAD)$p.value
    
    # Test of missingness in case of LOAD vs CTRL
    m.LOAD <- m[c(1,2),]
    p.chisq.LOADvsCTRL <- chisq.test(m.LOAD)$p.value
    valid.chisq.LOADvsCTRL <- as.numeric(sum(chisq.test(m.LOAD)$expected < 5) == 0)
    p.fisher.LOADvsCTRL <- fisher.test(m.LOAD)$p.value
    
    # results
    res <- c(count.total, count.na, rate.na, count.1, p.chisq, p.fisher, valid.chisq, p.chisq.EOADvsCTRL, p.fisher.EOADvsCTRL, valid.chisq.EOADvsCTRL, p.chisq.LOADvsCTRL, p.fisher.LOADvsCTRL, valid.chisq.LOADvsCTRL)
    res.withoutAdjustment <- data.frame(transcript.name = transcript.withNA[iter], Adjustment = "none", t(res))
    colnames(res.withoutAdjustment) <- c("transcript.name", "Selection", "N CTRL", "N LOAD", "N EOAD", "N total", "N NA CTRL", "N NA LOAD", "N NA EOAD", "N NA total", "Rate NA CTRL", "Rate NA LOAD", "Rate NA EOAD", "Rate NA total",
                                         "N 1 CTRL", "N 1 LOAD", "N 1 EOAD", "N 1 total",
                                         "DM for ordinal regression p-value (chi2)", "DM for ordinal regression p-value (fisher exact)", "DM for ordinal regression validity of chi2",
                                         "DM for EOAD vs CTRL p-value (chi2)", "DM for EOAD vs CTRL p-value (fisher exact)", "DM for EOAD vs CTRL validity of chi2",
                                         "DM for LOAD vs CTRL p-value (chi2)", "DM for LOAD vs CTRL p-value (fisher exact)", "DM for LOAD vs CTRL validity of chi2")
    

    ### WHEN ANALYSES WERE MADE WITH ADJUSTMENT ON APOE4
    
    genematrixAPOE4 <- subset(genematrix, apoe_genotype != "NA")
    
    # count available individuals
    count <- tapply(genematrixAPOE4[, transcript.withNA[iter]], genematrixAPOE4$Status, length)
    count.total <- c(count["Control"], count["LOAD"], count["EOAD"], sum(count))
    
    # count individuals with NA
    sum.na <- tapply(genematrixAPOE4[, transcript.withNA[iter]], genematrixAPOE4$Status, function(x){sum(is.na(x))})
    count.na <- c(sum.na["Control"], sum.na["LOAD"], sum.na["EOAD"], sum(sum.na))
    
    # count individuals with 1
    sum.1 <- tapply(genematrix[, transcript.withNA[iter]], genematrix$Status, function(x){sum((!is.na(x) & x==1))})
    count.1 <- c(sum.1["Control"], sum.1["LOAD"], sum.1["EOAD"], sum(sum.1))
    
    
    # rate individuals with NA
    rate.na <- count.na/count.total
    
    if (sum(sum.na) == 0) {
      res <- c(count.total, count.na, rate.na, rep(1, 9))
    } else {
      
      # matrice for chi2 and fisher exact test
      m <- matrix(c(sum.na[c("Control", "LOAD", "EOAD")], count[c("Control", "LOAD", "EOAD")]), byrow = FALSE, ncol = 2)
      
      # Test of missingness in case of ordinal regression
      p.chisq <- chisq.test(m)$p.value
      valid.chisq <- as.numeric(sum(chisq.test(m)$expected < 5) == 0)
      p.fisher <- fisher.test(m, simulate.p.value = TRUE, B = 50000)$p.value
      
      # Test of missingness in case of EOAD vs CTRL
      m.EOAD <- m[c(1,3),]
      p.chisq.EOADvsCTRL <- chisq.test(m.EOAD)$p.value
      valid.chisq.EOADvsCTRL <- as.numeric(sum(chisq.test(m.EOAD)$expected < 5) == 0)
      p.fisher.EOADvsCTRL <- fisher.test(m.EOAD)$p.valu
      
      # Test of missingness in case of LOAD vs CTRL
      m.LOAD <- m[c(1,2),]
      p.chisq.LOADvsCTRL <- chisq.test(m.LOAD)$p.value
      valid.chisq.LOADvsCTRL <- as.numeric(sum(chisq.test(m.LOAD)$expected < 5) == 0)
      p.fisher.LOADvsCTRL <- fisher.test(m.LOAD)$p.value
      
      # results
      res <- c(count.total, count.na, rate.na, count.1, p.chisq, p.fisher, valid.chisq, p.chisq.EOADvsCTRL, p.fisher.EOADvsCTRL, valid.chisq.EOADvsCTRL, p.chisq.LOADvsCTRL, p.fisher.LOADvsCTRL, valid.chisq.LOADvsCTRL)
    }
    
    res.withAdjustment <- data.frame(transcript.name = transcript.withNA[iter], Adjustment = "APOE available", t(res))
    colnames(res.withAdjustment) <- c("transcript.name", "Selection", "N CTRL", "N LOAD", "N EOAD", "N total", "N NA CTRL", "N NA LOAD", "N NA EOAD", "N NA total", "Rate NA CTRL", "Rate NA LOAD", "Rate NA EOAD", "Rate NA total",
                                      "N 1 CTRL", "N 1 LOAD", "N 1 EOAD", "N 1 total",
                                      "DM for ordinal regression p-value (chi2)", "DM for ordinal regression p-value (fisher exact)", "DM for ordinal regression validity of chi2",
                                      "DM for EOAD vs CTRL p-value (chi2)", "DM for EOAD vs CTRL p-value (fisher exact)", "DM for EOAD vs CTRL validity of chi2",
                                      "DM for LOAD vs CTRL p-value (chi2)", "DM for LOAD vs CTRL p-value (fisher exact)", "DM for LOAD vs CTRL validity of chi2")
    
    return(rbind(res.withoutAdjustment, res.withAdjustment))
    
  }
  
  # stopCluster(cl)
  
  # save results
  write.table(resultsB, paste0("DATA/Differential_missingness_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

}

