setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

dat.info <- read.table("DATA/SNV_LOF_Table.txt", header = TRUE, sep = "\t")

# Matrices <- import.matrice(type = "dosage", freq = 0.01, type.recouvrement = "reciproque")
# matrice_info <- Matrices$DosageMatrice[, c("sample", "Status", "apoe_genotype", "Age", "gender")]
# list.eoad <- subset(matrice_info, Status == "EOAD")$sample
# list.load <- subset(matrice_info, Status == "LOAD")$sample
# list.ctrl <- subset(matrice_info, Status == "Control")$sample

load("../clinData.Rdata")
matrice_info <- subset(clinData, exclusion == 0, select = c("vcfID", "Status"))
list.eoad <- subset(matrice_info, Status == "EOAD")$vcfID
list.load <- subset(matrice_info, Status == "LOAD")$vcfID
list.ctrl <- subset(matrice_info, Status == "Control")$vcfID

dat.info$eoad <- NA
dat.info$load <- NA

last <- dim(dat.info)[1]
for (i in 1:last) {
  oo <- as.character(dat.info[i, "cases"])
  if (substr(oo, nchar(oo), nchar(oo))==")") { next }
  oo2 <- gsub("[", "", oo, fixed = TRUE)
  oo3 <- gsub("]", "", oo2, fixed = TRUE)
  oo4 <- strsplit(oo3, " ")
  dat.info[i, "load"] <- sum(oo4[[1]] %in% list.load)
  dat.info[i, "eoad"] <- length(oo4[[1]]) - sum(oo4[[1]] %in% list.load)
}

for (i in 1:last) {
  oo <- as.character(dat.info[i, "controls"])
  if (substr(oo, nchar(oo), nchar(oo))==")") { next }
  oo2 <- gsub("[", "", oo, fixed = TRUE)
  oo3 <- gsub("]", "", oo2, fixed = TRUE)
  oo4 <- strsplit(oo3, " ")
  dat.info[i, "control_n"] <- sum(oo4[[1]] %in% list.ctrl)
}

dat.info2 <- subset(dat.info, exclude == "in analysis")

list.gene <- read.table("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/LOF_DATA/outputGenelist.txt")

dat.info3 <- subset(dat.info2, Gene %in% list.gene$V1)

write.table(dat.info3, "DATA/SNV_LOF_Table_simplified.txt", col.names = TRUE, row.names = FALSE, sep = "\t")


