setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

Matrices <- import.matrice(type = "CNV", freq = 0.01, type.recouvrement = "reciproque")

# DEL
type.CNV <- "DEL"
filtre.freq.value = 0.01
filtre.freq.type.recouvrement = "reciproque"
filtreCNV <- read.table(paste0("DATA/count_CNV_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)

# adi1
transcripts.list <- c("NM_001306077.2", "NM_018269.4")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set A. just 1 DEL
# analysis :
dat <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", "NM_001306077.2")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_001306077.2)
dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])

# PLIN4
transcripts.list <- c("NM_018269.4", "NM_001393888.1", "NM_001393889.1", "NM_001393890.1", "NM_001393891.1")
subset(filtreCNV, transcript.name %in% transcripts.list)
# analysis :
dat <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", "NM_018269.4")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_018269.4)
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])

# MUTS1
transcripts.list <- c("NM_001001924.3", "NM_001363057.2", "NM_001363058.2", "NM_001363059.2", "NM_001363061.2")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => not in set A

# GSDMD
transcripts.list <- c("NM_001166237.1")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => not in set A

# KLC3
transcripts.list <- c("NM_177417.3")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set A. just 2 DEL
# analysis :
dat <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", "NM_177417.3")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_177417.3)
dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])




# DUP
type.CNV <- "DUP"
filtre.freq.value = 0.01
filtre.freq.type.recouvrement = "reciproque"
filtreCNV <- read.table(paste0("DATA/count_CNV_", type.CNV, "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)


# PLIN4
transcripts.list <- c("NM_001393888.1", "NM_001393889.1", "NM_001393890.1", "NM_001393891.1")
subset(filtreCNV, transcript.name %in% transcripts.list)
# analysis :
dat <- Matrices$DelMatrice[, c("sample", "Status", "apoe_genotype", "NM_001393888.1")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_001393888.1)
dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])


# PRAMEF26
transcripts.list <- c("NM_001306072.3")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set B; just 1 DUP
# analysis :
dat <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", "NM_001306072.3")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_001306072.3)
dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])



# MUTS1
transcripts.list <- c("NM_001001924.3", "NM_001363057.2", "NM_001363058.2", "NM_001363059.2", "NM_001363061.2")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set B, just 1 DUP in LOAD
# analysis :
dat <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", "NM_001001924.3")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_001001924.3)
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])


# GSDMD
transcripts.list <- c("NM_001166237.1")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set B, just 2 DUP
# analysis :
dat <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", "NM_001166237.1")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_001166237.1)
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])


# MBL2
transcripts.list <- c("NM_000242.3", "NM_001378373.1", "NM_001378374.1")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set B, just 1 DUP
# analysis :
dat <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", "NM_000242.3")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_000242.3)
dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])

# FADS6
transcripts.list <- c("NM_178128.6")
subset(filtreCNV, transcript.name %in% transcripts.list)
# => in set B, just 1 DUP
# analysis :
dat <- Matrices$DupFullMatrice[, c("sample", "Status", "apoe_genotype", "NM_178128.6")]
dat$Status.f <- factor(dat$Status, levels = c("Control", "LOAD", "EOAD"))
dat$cnv <- as.factor(dat$NM_178128.6)
dat$APOE4 <- ifelse(dat$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(dat$apoe_genotype %in% c("24", "34"), 1, ifelse(dat$apoe_genotype == "44", 2, NA)))
# EOAD vs CTRL
dat.EOAD <- subset(dat, Status.f %in% c("EOAD", "Control"))
dat.EOAD$Y <- ifelse(dat.EOAD$Status.f == "EOAD", 1, 0)
mf <- logistf(Y ~ cnv, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat.EOAD)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# All AD vs CTRL
dat$Y <- ifelse(dat$Status.f == "Control", 0, 1)
mf <- logistf(Y ~ cnv, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])
# apoe
mf <- logistf(Y ~ cnv + APOE4, data = dat)
summary(mf)
c(exp(mf$coefficients["cnv1"]), exp(confint(mf))["cnv1",])


