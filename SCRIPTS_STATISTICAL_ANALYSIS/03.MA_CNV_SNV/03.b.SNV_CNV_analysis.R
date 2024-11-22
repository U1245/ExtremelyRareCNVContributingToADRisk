setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")


library("readr")

load("../MatriceResultat/matrice.Rdata")
filtre.freq.value = 0.01; filtre.freq.type.recouvrement = "reciproque"; LOAD90 = FALSE

### Extract gene of interest from diffrent sources
# CNV analysis
res1.del <- read.table(paste0("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_", ifelse(isTRUE(LOAD90), "LOADinf90_", ""), "and_subset_analyses_", "DOSAGE_DEL_DUP", "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
list.gene.del_CNV <- c("PRAMEF26", "ADI1", "MTUS1", "MBL2", "FADS6", "PLIN4", "KLC3", "ERCC2", "ZNF74", "SCARF2", "KLHL22", "MED15", "PI4KA", "SERPIND1", "SNAP29", "CRKL", "LZTR1", "THAP7", "P2RX6", "SLC7A4", "LRRC74B")

# Holstege et al exome
list.gene.SNV <- c("SORL1", "TREM2", "ABCA7", "ATP8B4", "ABCA1", "ADAM10", "CBX3", "B3GNT4", "SRC", "MOT1", "SLC16A1", "PRSS3", "TYROBP")
# GWAS Holstege et al
list.gene.GWAS <- as.character(read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt")$V1)
# Combine all genes
list.gene.del <- unique(c(list.gene.del_CNV, list.gene.SNV, list.gene.GWAS))

# load other useful information
load("../GeneTranscript.Rdata")
load("../clinData.Rdata")
B  <- length(list.gene.del)
resultsB <- foreach(iter = 1:B, .combine = "rbind") %do% {
  gene.name <- list.gene.del[iter]
  # for each gene, make the analysis
  # for (gene.name in list.gene.del) {
  cat("-----------------\n", gene.name, "\n----------------\n")
  transcript.list1 <- unique(subset(GeneTranscript, Gene == gene.name)$Transcript)
  transcript.list <- transcript.list1[substr(transcript.list1, 1, 2) == "NM"]
  genematrix <- DelMatrice[, c("sample", "datasetid", "Status", transcript.list)]
  genematrix.vcfiid <- merge(genematrix, clinData, by = c("sample", "Status", "datasetid"), all.x = TRUE, all.y = FALSE)
  genematrix.vcfiid$status_vcf <- ifelse(genematrix.vcfiid$Status == "Control", "control", ifelse(genematrix.vcfiid$Status == "LOAD", "case>=65", ifelse(genematrix.vcfiid$Status == "EOAD", "case<65", NA)))
  genematrix.vcfiid$vcfID_ok <- ifelse(genematrix.vcfiid$vcfID == ".", genematrix.vcfiid$sample, genematrix.vcfiid$vcfID)
  # Pour chaque gene liste des transcripts avec <= 25% de missingness
  transcript.list.miss <- c()
  for (i in transcript.list) {
    if (sum(is.na(genematrix[, i]))/(dim(genematrix)[1]) <= 0.25) {
      transcript.list.miss <- c(transcript.list.miss, i)
    }
  }
  genematrix.vcfiid.right.columns <- subset(genematrix.vcfiid, select = c("vcfID_ok", "Status", "datasetid", "status_vcf", transcript.list.miss))
  genematrix.vcfiid.right.columns.no.NA <- na.omit(genematrix.vcfiid.right.columns)
  # count by Gene
  if (length(transcript.list.miss) == 1) {
    genematrix.vcfiid.right.columns.no.NA$cnv <- apply(matrix(genematrix.vcfiid.right.columns.no.NA[, transcript.list.miss], ncol = 1), 1, function(x){as.numeric(sum(x)>0)})
  } else {
    genematrix.vcfiid.right.columns.no.NA$cnv <- apply(genematrix.vcfiid.right.columns.no.NA[, transcript.list.miss], 1, function(x){as.numeric(sum(x)>0)})
  }
  
  # table(genematrix.vcfiid.right.columns.no.NA$Status, genematrix.vcfiid.right.columns.no.NA$cnv)
  
  
  # lire individus exclus
  system(paste0('grep "list removed sample" ../LOF_DATA/LOG/log', gene.name, '.txt | sed -e "s/;/\\n/g" | sed -e "s/.*=//" > ../LOF_DATA/LOG/list.removed.sample.', gene.name, '.txt'))
  my_file <- try(readLines(paste0("../LOF_DATA/LOG/list.removed.sample.", gene.name, ".txt")))
  if (length(my_file) == 0) {
    res <- rep(NA, 61)
    excl <- "no SNV"
  } else {
    if (sum(my_file != "")>0) {
      excluded.individual.snv <- read.table(paste0("../LOF_DATA/LOG/list.removed.sample.", gene.name, ".txt"), stringsAsFactors = FALSE)$V1
      datfinal <-  subset(genematrix.vcfiid.right.columns.no.NA, !(vcfID_ok %in% excluded.individual.snv))
    } else {
      datfinal <-  genematrix.vcfiid.right.columns.no.NA
    }
    # SNV LoF
    snv <- read.table(paste0("../LOF_DATA/RESULT/", gene.name, "_variants.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # position
    snv$position <- sapply(strsplit(snv$variant_id, split = ":"), "[", 2)
    
    pp <- try(readLines(paste0("../LOF_DATA/RESULT/varToExclude.", gene.name, ".tsv")))
    if (length(pp) == 0) {
      position.excluded <- c()
    } else {
      position.excluded <- unique(as.character(read.table(paste0("../LOF_DATA/RESULT/varToExclude.", gene.name, ".tsv"), stringsAsFactors = FALSE)$V2))
    }
    # extract only snv of interest
    snv.interest <- subset(snv, (is.na(exclude_reason) | exclude_reason == "") & !(position %in% position.excluded))
    
    # extract carriers
    idcases = unlist(strsplit(gsub("\\[", "", gsub("]", "", snv.interest$cases)), split = " "))
    idctrls = unlist(strsplit(gsub("\\[", "", gsub("]", "", snv.interest$controls)), split = " "))
    
    datfinal$snv <- ifelse(datfinal$vcfID_ok %in% unique(c(idcases, idctrls)), 1, 0)
    datfinal$LOF <- as.factor(ifelse(datfinal$cnv+datfinal$snv>0, 1, 0))

    
    datfinal$Status.f <- factor(datfinal$Status, c("Control", "LOAD", "EOAD"))
    dat <- datfinal 
    count <- addmargins(table(dat$LOF, dat$Status.f))
    
    if (dim(count)[1] < 3) {
      res <- rep(NA, 61)
      excl <- "no LoF"
    } else {
      
      # OR AD
      dat.AD <- dat
      dat.AD$Y <- ifelse(dat.AD$Status.f == "Control", 0, 1)
      m <- glm(Y ~ LOF, data = dat.AD, family = binomial(link = "logit"))
      OR.AD <- c(exp(m$coef)["LOF1"], exp(confint(m))["LOF1",], summary(m)$coefficients["LOF1", "Pr(>|z|)"])

      mf <- logistf(Y ~ LOF, data = dat.AD)
      ORf.AD <- c(exp(mf$coefficients["LOF1"]), exp(confint(mf))["LOF1",], mf$prob["LOF1"])

      
      # OR EOAD
      dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
      if (sum(dat.EOAD$LOF == 1, na.rm = TRUE) == 0) {
        OR.EOAD <- c(NA, NA, NA, NA)
        ORf.EOAD <- c(NA, NA, NA, NA)
      } else {
        dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
        m <- glm(Y ~ LOF, data = dat.EOAD, family = binomial(link = "logit"))
        OR.EOAD <- c(exp(m$coef)["LOF1"], exp(confint(m))["LOF1",], summary(m)$coefficients["LOF1", "Pr(>|z|)"])

        mf <- logistf(Y ~ LOF, data = dat.EOAD)
        ORf.EOAD <- c(exp(mf$coefficients["LOF1"]), exp(confint(mf))["LOF1",], mf$prob["LOF1"])
      }
      
      # OR LOAD
      dat.LOAD <- subset(dat, Status.f %in% c("LOAD", "Control"))
      if (sum(dat.LOAD$LOF == 1, na.rm = TRUE) == 0) {
        OR.LOAD <- c(NA, NA, NA, NA)
        ORf.LOAD <- c(NA, NA, NA, NA)
      } else {
        dat.LOAD$Y <- ifelse(dat.LOAD$Status.f == "LOAD", 1, 0)
        m <- glm(Y ~ LOF, data = dat.LOAD, family = binomial(link = "logit"))
        OR.LOAD <- c(exp(m$coef)["LOF1"], exp(confint(m))["LOF1",], summary(m)$coefficients["LOF1", "Pr(>|z|)"])
        
        mf <- logistf(Y ~ LOF, data = dat.LOAD)
        ORf.LOAD <- c(exp(mf$coefficients["LOF1"]), exp(confint(mf))["LOF1",], mf$prob["LOF1"])
      }
      
      # Ordinal regression avec clm
      rego <- clm(Status.f ~ LOF, data = na.omit(dat))
      confint.rego <- try(exp(confint(rego))["LOF1",], silent = TRUE)
      if (inherits(confint.rego, "try-error")) {
        OR.clm <- c(exp(rego$coefficients["LOF1"]), NA, NA)
      } else {
        OR.clm <- c(exp(rego$coefficients["LOF1"]), confint.rego)
      }
      
      pvalue.clm <- summary(rego)$coefficients["LOF1", "Pr(>|z|)"]
      message.clm <- rego$message
      
      # Ordinal regression avec polr
      rego2 <- polr(Status.f ~ LOF, data = na.omit(dat), Hess = TRUE)
      confint.rego2 <- try(exp(confint(rego2)), silent = TRUE)
      if (inherits(confint.rego2, "try-error")) {
        OR.polr <- c(exp(rego2$coefficients["LOF1"]), NA, NA)
      } else {
        OR.polr <- c(exp(rego2$coefficients["LOF1"]), confint.rego2)
      }
      
      rego0 <- polr(Status.f ~ 1, data = na.omit(dat))
      statRV <- -2*(logLik(rego0) - logLik(rego2))
      pvalue.polr <- 1 - pchisq(statRV, df = length(coef(rego2)) - length(coef(rego0)))
      
      
      count <- addmargins(table(dat$LOF, dat$Status.f))
      
      LOF.N <- c(count["Sum",], count["1",])
      LOF.rate <- c(count["1",]/count["Sum",])
      
      count.snv <- addmargins(table(dat$snv, dat$Status.f))
      count.cnv <- addmargins(table(dat$cnv, dat$Status.f))
      
      
      if (dim(count.snv)[1] < 3) {
        res <- rep(NA, 61)
        excl <- "no SNV"
      } else if (dim(count.cnv)[1] < 3) {
        res <- rep(NA, 61)
        excl <- "no CNV"
      } else {
        snv.N <- c(count.snv["1",])
        snv.rate <- c(count.snv["1",]/count.snv["Sum",])
        cnv.N <- c(count.cnv["1",])
        cnv.rate <- c(count.cnv["1",]/count.cnv["Sum",])
        
        res <- c(LOF.N, LOF.rate, snv.N, snv.rate, cnv.N, cnv.rate, OR.AD, ORf.AD, OR.EOAD, ORf.EOAD, OR.LOAD, ORf.LOAD, OR.clm, OR.polr, pvalue.clm, pvalue.polr, message.clm)
        excl <- ""
      }
      
    }
  }
  res.withoutAdjustment <- data.frame(Gene = gene.name, Transcripts = paste(transcript.list.miss, collapse = ";"), exc = excl, Adjustment = "none", t(res), stringsAsFactors = FALSE)
  colnames(res.withoutAdjustment) <- c("Gene", "transcript.name", "reason exclusion", "Adjustment", "N CTRL", "N LOAD", "N EOAD", "N total", "Ncarriers CTRL", 
                                       "Ncarriers LOAD", "Ncarriers EOAD", "Ncarriers total", "Rate CTRL", "Rate LOAD", "Rate EOAD", "Rate total", 
                                       "Ncarriers snv CTRL", "Ncarriers snv LOAD", "Ncarriers snv EOAD", "Ncarriers snv total", "Rate snv CTRL", "Rate snv LOAD", "Rate snv EOAD", "Rate snv total",
                                       "Ncarriers cnv CTRL", "Ncarriers cnv LOAD", "Ncarriers cnv EOAD", "Ncarriers cnv total", "Rate cnv CTRL", "Rate cnv LOAD", "Rate cnv EOAD", "Rate cnv total",
                                       "OR allAD vs CTRL", "OR allAD vs CTRL - CIinf", "OR allAD vs CTRL - CIsupOR", "p-value allAD vs CTRL",
                                       "OR firth allAD vs CTRL", "OR firth allAD vs CTRL - CIinf", "OR firth allAD vs CTRL - CIsupOR", "p-value firth allAD vs CTRL",
                                       "OR EOAD vs CTRL", "OR EOAD vs CTRL - CIinf", "OR EOAD vs CTRL - CIsupOR", "p-value EOAD vs CTRL",
                                       "OR firth EOAD vs CTRL", "OR firth EOAD vs CTRL - CIinf", "OR firth EOAD vs CTRL - CIsupOR", "p-value firth EOAD vs CTRL",
                                       "OR LOAD vs CTRL", "OR LOAD vs CTRL - CIinf", "OR LOAD vs CTRL - CIsup", "p-value LOAD vs CTRL",
                                       "OR firth LOAD vs CTRL", "OR firth LOAD vs CTRL - CIinf", "OR firth LOAD vs CTRL - CIsup", "p-value firth LOAD vs CTRL",
                                       "OR clm", "OR clm - CIinf", "OR clm - CIsup", "OR polr", "OR polr - CIinf", "OR polr - CIsup", "p-value (ordinal regression, using clm)", "p-value (ordinal regression, using polr)", "warning clm")
  
  return(res.withoutAdjustment)
}

write.table(resultsB, "RESULTS/Analysis_DEL_SNV_LOF_ordinal_regression_and_subset_analyses_filtre_freq_0.01_recouvrement_reciproque_add_TYROBP.txt", col.names = TRUE, row.names = FALSE, sep = "\t")



