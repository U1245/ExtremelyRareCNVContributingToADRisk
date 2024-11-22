

source("SCRIPTS/02.Transcript_analyses/02.e.Transcript_union_setA_B_dosage.analysis.R")
# recouvrement réciproque, freq 1%
# DEL/DUP
dosage.analysis_union_setA_B(type.CNV = "DEL_DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 

### table sup du papier

dat.papier1 <- read.table("RESULTS/Table_for_article_october2023_Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")
dat.papier2 <- subset(dat.papier1, Adjustment == "none" & substr(transcript.name, 1, 2) == "NM" & Ncarriers.CNV.total >= 4)


