

############## UNION set A et B
dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

### Si on ne garde que les CNVs vus au moins 4 fois et 1 test par Gene x N


sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >=4)
dim(sub)
length(unique(sub$Gene))


sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.EOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)

sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)

dim(sub.uniqueGene)

summary(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL)

sub.uniqueGene$pvalue.fdr <- p.adjust(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL, method = "fdr")

head(sub.uniqueGene[order(sub.uniqueGene$pvalue.fdr), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)
head(sub.uniqueGene[order(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)




############## DELETIONS only
dat <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DEL_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

### Si on ne garde que les CNVs vus au moins 4 fois et 1 test par Gene x N


sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.total >=4)
dim(sub)
length(unique(sub$Gene)) #1043


sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CTRL, sub$Ncarriers.EOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)

sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)

dim(sub.uniqueGene) #1100

summary(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL)

sub.uniqueGene$pvalue.fdr <- p.adjust(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL, method = "fdr")

head(sub.uniqueGene[order(sub.uniqueGene$pvalue.fdr), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)
head(sub.uniqueGene[order(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)





############## DUPLICATION only
dat <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

### Si on ne garde que les CNVs vus au moins 4 fois et 1 test par Gene x N


sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.total >=4)
dim(sub)
length(unique(sub$Gene)) #1252


sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CTRL, sub$Ncarriers.EOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)

sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)

dim(sub.uniqueGene) #1331

summary(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL)

sub.uniqueGene$pvalue.fdr <- p.adjust(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL, method = "fdr")

head(sub.uniqueGene[order(sub.uniqueGene$pvalue.fdr), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)
head(sub.uniqueGene[order(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)




############## UNION set A et B adjusted for PC
dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_PC_adjusted_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

### Si on ne garde que les CNVs vus au moins 4 fois et 1 test par Gene x N


sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >=4)
dim(sub)
length(unique(sub$Gene))


sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.EOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)

sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)

dim(sub.uniqueGene)

summary(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL)

sub.uniqueGene$pvalue.fdr <- p.adjust(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL, method = "fdr")

head(sub.uniqueGene[order(sub.uniqueGene$pvalue.fdr), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)
head(sub.uniqueGene[order(sub.uniqueGene$p.value.firth.EOAD.vs.CTRL), c("Gene", "p.value.firth.EOAD.vs.CTRL", "pvalue.fdr")], 25)



################## TESTS SEPTEMBRE 2024 : FDR LOAD ves CTRL et All AD vs CTRL

############## UNION set A et B LOAD vs CTRL
dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

### Si on ne garde que les CNVs vus au moins 4 fois et 1 test par Gene x N


sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >=4)
dim(sub)
length(unique(sub$Gene))


sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.LOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.LOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)

sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)

dim(sub.uniqueGene)

summary(sub.uniqueGene$p.value.firth.LOAD.vs.CTRL)

sub.uniqueGene$pvalue.fdr <- p.adjust(sub.uniqueGene$p.value.firth.LOAD.vs.CTRL, method = "fdr")

head(sub.uniqueGene[order(sub.uniqueGene$pvalue.fdr), c("Gene", "p.value.firth.LOAD.vs.CTRL", "pvalue.fdr")], 25)
head(sub.uniqueGene[order(sub.uniqueGene$p.value.firth.LOAD.vs.CTRL), c("Gene", "p.value.firth.LOAD.vs.CTRL", "pvalue.fdr")], 25)

# exploration des genes qui sortent MUC5B et SLC25A24 


############## UNION set A et B all AD vs CTRL
dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

### Si on ne garde que les CNVs vus au moins 4 fois et 1 test par Gene x N


sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >=4)
dim(sub)
length(unique(sub$Gene))


sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$N.LOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.EOAD, sub$Ncarriers.CNV.LOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)

sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)

dim(sub.uniqueGene)

summary(sub.uniqueGene$p.value.firth.all.AD.vs.CTRL)

sub.uniqueGene$pvalue.fdr <- p.adjust(sub.uniqueGene$p.value.firth.all.AD.vs.CTRL, method = "fdr")

head(sub.uniqueGene[order(sub.uniqueGene$pvalue.fdr), c("Gene", "p.value.firth.all.AD.vs.CTRL", "pvalue.fdr")], 25)
head(sub.uniqueGene[order(sub.uniqueGene$p.value.firth.all.AD.vs.CTRL), c("Gene", "p.value.firth.all.AD.vs.CTRL", "pvalue.fdr")], 25)

# rien ne sort en all AD vs CTRL
