

source("SCRIPTS/04.Score_analyses/04.b.List.of.genes.analysis_PC.R")


### GWAS - EADB
# recouvrement réciproque, freq 1%
# DEL
list_of_genes.analysis_PC(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", list_of_genes = "GWAS_EADB", p.value_threshold = 1e-5) 





