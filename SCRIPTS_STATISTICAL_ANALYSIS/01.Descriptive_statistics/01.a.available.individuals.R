# Load packages and dataset

setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")
Matrices <- import.matrice(type = "CNV", freq = 0.01, type.recouvrement = "reciproque") 
# describe individuals according to Status, Age, gender, APOE

# DELETION
res <- data.frame(Dataset = as.character(), Group = as.character(), N = as.double(), age.available = as.double(), age.na = as.double(), age.mean = as.double(), age.sd = as.double(), 
                  age.min = as.double(), age.max = as.double(), age.median = as.double(), age.Q1 = as.double(), age.Q3 = as.double(), 
                  gender.available = as.double(), gender.na = as.double(), gender.1.N = as.double(), gender.1.percent = as.double(), gender.2.N = as.double(), gender.2.percent = as.double(),
                  APOE.available = as.double(), APOE.na = as.double(), APOE.22.N = as.double(), APOE.22.percent = as.double(), APOE.23.N = as.double(), APOE.23.percent = as.double(),
                  APOE.24.N = as.double(), APOE.24.percent = as.double(), APOE.33.N = as.double(), APOE.33.percent = as.double(), APOE.34.N = as.double(), APOE.34.percent = as.double(),
                  APOE.44.N = as.double(), APOE.44.percent = as.double(), APOE4.non.carrier.N = as.double(), APOE4.non.carrier.percent = as.double(),
                  APOE4.heterozygous.N = as.double(), APOE4.heterozygous.percent = as.double(), APOE4.homozygous.N = as.double(), APOE4.homozygous.percent = as.double(), stringsAsFactors = FALSE)

DelMatrice2 <- Matrices$DelMatrice[,1:nbNonGenes]
DelMatrice2$APOE4 <- ifelse(DelMatrice2$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(DelMatrice2$apoe_genotype %in% c("24", "34"), 1, ifelse(DelMatrice2$apoe_genotype %in% c("44"), 2, NA)))
DelMatrice2$APOENA <- ifelse(DelMatrice2$apoe_genotype == "NA", NA, DelMatrice2$apoe_genotype)
DelMatrice2$genderNA <- DelMatrice2$gender
table(DelMatrice2$APOENA, DelMatrice2$apoe_genotype)
table(DelMatrice2$APOE4, DelMatrice2$apoe_genotype)
table(DelMatrice2$APOE4, DelMatrice2$APOENA)
table(DelMatrice2$genderNA, DelMatrice2$gender)


# Deletion All
k <- 1; genematrix <- DelMatrice2
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "all individuals"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]

res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion EOAD
k <- 2; genematrix <- subset(DelMatrice2, Status == "EOAD")
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "EOAD"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion LOAD
k <- 3; genematrix <- subset(DelMatrice2, Status == "LOAD")
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "LOAD"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion LOAD < 90
k <- 4; genematrix <- subset(DelMatrice2, Status == "LOAD" & !is.na(Age) & as.numeric(Age) < 90)
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "LOAD < 90"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion LOAD < 90
k <- 5; genematrix <- subset(DelMatrice2, Status == "LOAD" & !is.na(Age) & as.numeric(Age) >= 90)
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "LOAD >= 90"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion LOAD age NA
k <- 6; genematrix <- subset(DelMatrice2, Status == "LOAD" & is.na(Age))
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "LOAD with age NA"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion CTRL
k <- 7; genematrix <- subset(DelMatrice2, Status == "Control")
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "CTRL"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion All AD
k <- 8; genematrix <- subset(DelMatrice2, Status != "Control")
res[k, "Dataset"] <- "DelMatrice"
res[k, "Group"] <- "All AD"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

res.DEL <- res

# DUPLICATION

res <- data.frame(Dataset = as.character(), Group = as.character(), N = as.double(), age.available = as.double(), age.na = as.double(), age.mean = as.double(), age.sd = as.double(), 
                  age.min = as.double(), age.max = as.double(), age.median = as.double(), age.Q1 = as.double(), age.Q3 = as.double(), 
                  gender.available = as.double(), gender.na = as.double(), gender.1.N = as.double(), gender.1.percent = as.double(), gender.2.N = as.double(), gender.2.percent = as.double(),
                  APOE.available = as.double(), APOE.na = as.double(), APOE.22.N = as.double(), APOE.22.percent = as.double(), APOE.23.N = as.double(), APOE.23.percent = as.double(),
                  APOE.24.N = as.double(), APOE.24.percent = as.double(), APOE.33.N = as.double(), APOE.33.percent = as.double(), APOE.34.N = as.double(), APOE.34.percent = as.double(),
                  APOE.44.N = as.double(), APOE.44.percent = as.double(), APOE4.non.carrier.N = as.double(), APOE4.non.carrier.percent = as.double(),
                  APOE4.heterozygous.N = as.double(), APOE4.heterozygous.percent = as.double(), APOE4.homozygous.N = as.double(), APOE4.homozygous.percent = as.double(), stringsAsFactors = FALSE)

DupFullMatrice2 <- Matrices$DupFullMatrice[,1:nbNonGenes]
DupFullMatrice2$APOE4 <- ifelse(DupFullMatrice2$apoe_genotype %in% c("22", "23", "33"), 0, ifelse(DupFullMatrice2$apoe_genotype %in% c("24", "34"), 1, ifelse(DupFullMatrice2$apoe_genotype %in% c("44"), 2, NA)))
DupFullMatrice2$APOENA <- ifelse(DupFullMatrice2$apoe_genotype == "NA", NA, DupFullMatrice2$apoe_genotype)
DupFullMatrice2$genderNA <- DupFullMatrice2$gender
table(DupFullMatrice2$APOE4, DupFullMatrice2$apoe_genotype)
table(DupFullMatrice2$APOE4, DupFullMatrice2$APOENA)
table(DupFullMatrice2$genderNA, DupFullMatrice2$gender)

# Duplication All
k <- 1; genematrix <- DupFullMatrice2
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "all individuals"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Duplication EOAD
k <- 2; genematrix <- subset(DupFullMatrice2, Status == "EOAD")
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "EOAD"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Duplication LOAD
k <- 3; genematrix <- subset(DupFullMatrice2, Status == "LOAD")
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "LOAD"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Duplication LOAD < 90
k <- 4; genematrix <- subset(DupFullMatrice2, Status == "LOAD" & !is.na(Age) & as.numeric(Age) < 90)
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "LOAD < 90"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Duplication LOAD < 90
k <- 5; genematrix <- subset(DupFullMatrice2, Status == "LOAD" & !is.na(Age) & as.numeric(Age) >= 90)
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "LOAD >= 90"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Duplication LOAD age NA
k <- 6; genematrix <- subset(DupFullMatrice2, Status == "LOAD" & is.na(Age))
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "LOAD with age NA"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Duplication CTRL
k <- 7; genematrix <- subset(DupFullMatrice2, Status == "Control")
res[k, "Dataset"] <- "DupFullMatrice"
res[k, "Group"] <- "CTRL"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))

# Deletion All AD
k <- 8; genematrix <- subset(DupFullMatrice2, Status != "Control")
res[k, "Dataset"] <- "DupMatrice"
res[k, "Group"] <- "All AD"
res[k, "N"] <- dim(genematrix)[1]
res[k, "age.available"] <- length(na.omit(genematrix$Age))
res[k, "age.na"] <- sum(is.na(genematrix$Age))
res[k, "age.mean"] <- mean(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.sd"] <- sd(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.min"] <- min(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.max"] <- max(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.median"] <- median(as.numeric(genematrix$Age), na.rm = TRUE)
res[k, "age.Q1"] <- summary(as.numeric(genematrix$Age))[2]
res[k, "age.Q3"] <- summary(as.numeric(genematrix$Age))[5]
res[k, "gender.available"] <- length(na.omit(genematrix$genderNA))
res[k, "gender.na"] <- sum(is.na(genematrix$genderNA))
res[k, "gender.1.N"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)
res[k, "gender.1.percent"] <- sum(genematrix$genderNA == "1", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "gender.2.N"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)
res[k, "gender.2.percent"] <- sum(genematrix$genderNA == "2", na.rm = TRUE)/length(na.omit(genematrix$genderNA))
res[k, "APOE.available"] <- length(na.omit(genematrix$APOENA))
res[k, "APOE.na"] <- sum(is.na(genematrix$APOENA))
res[k, "APOE.22.N"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)
res[k, "APOE.22.percent"] <- sum(genematrix$APOENA == "22", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.23.N"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)
res[k, "APOE.23.percent"] <- sum(genematrix$APOENA == "23", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.24.N"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)
res[k, "APOE.24.percent"] <- sum(genematrix$APOENA == "24", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.33.N"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)
res[k, "APOE.33.percent"] <- sum(genematrix$APOENA == "33", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.34.N"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)
res[k, "APOE.34.percent"] <- sum(genematrix$APOENA == "34", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE.44.N"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)
res[k, "APOE.44.percent"] <- sum(genematrix$APOENA == "44", na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.non.carrier.N"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)
res[k, "APOE4.non.carrier.percent"] <- sum(genematrix$APOE4 == 0, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.heterozygous.N"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)
res[k, "APOE4.heterozygous.percent"] <- sum(genematrix$APOE4 == 1, na.rm = TRUE)/length(na.omit(genematrix$APOENA))
res[k, "APOE4.homozygous.N"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)
res[k, "APOE4.homozygous.percent"] <- sum(genematrix$APOE4 == 2, na.rm = TRUE)/length(na.omit(genematrix$APOENA))


res.DUP <- res

# save results

res.ALL <- rbind(res.DEL, res.DUP)
write.table(res.ALL, "RESULTS/Descriptive_statistics.txt", sep = "\t", col.names = TRUE, row.names = FALSE)


