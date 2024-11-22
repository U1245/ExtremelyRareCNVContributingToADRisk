setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

source("SCRIPTS/02.Transcript_analyses/02.a.Transcript.selection.R")
# recouvrement réciproque, freq 1%
transcript.selection(freq = 0.01, type.recouvrement = "reciproque")

# # recouvrement réciproque, freq 10%
# transcript.selection(freq = 0.1, type.recouvrement = "reciproque")
# 
# 
# # recouvrement unilatéral, freq 1%
# transcript.selection(freq = 0.01, type.recouvrement = "unilateral")
# 
# # recouvrement unilatéral, freq 10%
# transcript.selection(freq = 0.1, type.recouvrement = "unilateral")


## count exclusion for dosage analysis
info.del <- read.table("DATA/count_CNV_DEL_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)
info.dup <- read.table("DATA/count_CNV_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", sep = "\t", header = TRUE)

info <- merge(info.del, info.dup, by = "transcript.name")
length(info$transcript.name) #145,522
# nb not excluded for frequency
length(subset(info, exclusion.x == FALSE & exclusion.y == FALSE)$transcript.name) 
# 140,047
length(subset(info, exclusion.x == FALSE & exclusion.y == FALSE & (count.cnv.x+count.cnv.y>=3))$transcript.name)
# 14,312
length(subset(info, exclusion.x == FALSE & exclusion.y == FALSE & (count.cnv.x+count.cnv.y==0))$transcript.name)
# 94,915
length(subset(info, exclusion.x == FALSE & exclusion.y == FALSE & ((count.cnv.x+count.cnv.y==1) | (count.cnv.x+count.cnv.y==2)) )$transcript.name)
# 30,820
length(subset(info, exclusion.x == FALSE & exclusion.y == FALSE & (count.cnv.x+count.cnv.y>=3) & count.cnv.x >=1 & count.cnv.y >=1)$transcript.name)
# 7,775
length(subset(info, exclusion.x == FALSE & exclusion.y == FALSE & (count.cnv.x+count.cnv.y>=3) & (count.cnv.x ==0 | count.cnv.y ==0))$transcript.name)
# 6,537