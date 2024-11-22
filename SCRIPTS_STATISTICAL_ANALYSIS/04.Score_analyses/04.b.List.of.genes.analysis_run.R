

source("SCRIPTS/04.Score_analyses/04.b.List.of.genes.analysis.R")

### ABeta

# recouvrement réciproque, freq 1%
# DEL
list_of_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "ABeta", p.value_threshold = 1e-5) 
# complete DUP
list_of_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "ABeta", p.value_threshold = 1e-5) 
# DEL+partial DUP
list_of_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "ABeta", p.value_threshold = 1e-5) 



### GWAS - EADB

# recouvrement réciproque, freq 1%
# DEL
list_of_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "GWAS_EADB", p.value_threshold = 1e-5) 
# complete DUP
list_of_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "GWAS_EADB", p.value_threshold = 1e-5) 
# DEL+partial DUP
list_of_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "GWAS_EADB", p.value_threshold = 1e-5) 



### All genes

# recouvrement réciproque, freq 1%
# DEL
list_of_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "All", p.value_threshold = 1e-5) 
# complete DUP
list_of_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "All", p.value_threshold = 1e-5) 
# DEL+partial DUP
list_of_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "All", p.value_threshold = 1e-5) 



### Interprétation des résultats
list.EADB <- unique(read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt")$V1)
load("../GeneTranscript.Rdata")
dat.EADB <- rbind(subset(GeneTranscript, Gene %in% list.EADB & substr(Transcript, 1, 2) == "NM"), data.frame(Gene = list.EADB[!(list.EADB %in% GeneTranscript$Gene)], Transcript = NA))
info.count <- read.table("DATA/count_CNV_DEL_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE)
dat.EADB.count <- merge(dat.EADB, info.count, by.x = "Transcript", by.y = "transcript.name", all.x = TRUE, all.y = FALSE)
info.DM <- subset(read.table("DATA/Differential_missingness_DEL_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE), Selection == "none")
dat.EADB.count.DM <- merge(dat.EADB.count, subset(info.DM, select = c(transcript.name, DM.for.EOAD.vs.CTRL.p.value..chi2.)), by.x = "Transcript", by.y = "transcript.name", all.x = TRUE, all.y = FALSE)
load("../gnomadGeneFreq.RData")
dat.EADB.count.DM.freq <- merge(dat.EADB.count.DM, subset(freqGene, type == "DEL", select = c(Transcript, cumulDGV_AF, cumulEUR_AF)), by = "Transcript", all.x = TRUE, all.y = FALSE)
dim(dat.EADB.count.DM.freq)          
head(dat.EADB.count.DM.freq)
length(unique(dat.EADB.count.DM.freq$Gene))
length(list.EADB)
write.table(dat.EADB.count.DM.freq, "DATA/details_about_EADB_GWAS_list2.txt", col.names = TRUE, row.names = FALSE, sep ="\t")


### Interprétation des résultats
list.EADB <- unique(read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt")$V1)
load("../GeneTranscript.Rdata")
dat.EADB <- rbind(subset(GeneTranscript, Gene %in% list.EADB & substr(Transcript, 1, 2) == "NM"), data.frame(Gene = list.EADB[!(list.EADB %in% GeneTranscript$Gene)], Transcript = NA))
info.count <- read.table("DATA/count_CNV_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE)
dat.EADB.count <- merge(dat.EADB, info.count, by.x = "Transcript", by.y = "transcript.name", all.x = TRUE, all.y = FALSE)
info.DM <- subset(read.table("DATA/Differential_missingness_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE), Selection == "none")
dat.EADB.count.DM <- merge(dat.EADB.count, subset(info.DM, select = c(transcript.name, DM.for.EOAD.vs.CTRL.p.value..chi2.)), by.x = "Transcript", by.y = "transcript.name", all.x = TRUE, all.y = FALSE)
load("../gnomadGeneFreq.RData")
dat.EADB.count.DM.freq <- merge(dat.EADB.count.DM, subset(freqGene, type == "DEL", select = c(Transcript, cumulDGV_AF, cumulEUR_AF)), by = "Transcript", all.x = TRUE, all.y = FALSE)
dim(dat.EADB.count.DM.freq)          
head(dat.EADB.count.DM.freq)
length(unique(dat.EADB.count.DM.freq$Gene))
length(list.EADB)
write.table(dat.EADB.count.DM.freq, "DATA/details_about_EADB_GWAS_list2.txt", col.names = TRUE, row.names = FALSE, sep ="\t")


### quels sont les genes deletes ?
EADB <- read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
genes.EADB <- EADB$V1






