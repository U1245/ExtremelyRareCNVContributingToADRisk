setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

Matrices <- import.matrice(type = "dosage", freq = 0.01, type.recouvrement = "reciproque")

cn <- colnames(Matrices$DosageMatrice)[1:13]
dat <- Matrices$DosageMatrice[,c(cn, "NM_178128.6")]

res <- table(dat$Status, dat$NM_178128.6)
res

# modèle de base

dat.EOAD <- subset(dat, Status!="LOAD")
dat.EOAD$Y <- ifelse(dat.EOAD$Status == "EOAD", 1, 0)
mf <- logistf(Y ~ NM_178128.6, data = dat.EOAD)
summary(mf)
ORf.EOAD <- c(exp(mf$coefficients["NM_178128.6"]), exp(confint(mf))["NM_178128.6",])
pvaluef.EOAD <- mf$prob["NM_178128.6"]
ORf.EOAD
pvaluef.EOAD

# avec suppression du cnv faux positif

dat[dat$sample=="EXT-DUP-1506-001","NM_178128.6"] <- 2

res2 <- table(dat$Status, dat$NM_178128.6)
res2

dat.EOAD <- subset(dat, Status!="LOAD")
dat.EOAD$Y <- ifelse(dat.EOAD$Status == "EOAD", 1, 0)
mf <- logistf(Y ~ NM_178128.6, data = dat.EOAD)
summary(mf)
ORf.EOAD <- c(exp(mf$coefficients["NM_178128.6"]), exp(confint(mf))["NM_178128.6",])
pvaluef.EOAD <- mf$prob["NM_178128.6"]
ORf.EOAD
pvaluef.EOAD

