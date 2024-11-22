

source("SCRIPTS/04.Score_analyses/04.a.Differential.missingness.R")
# recouvrement réciproque, freq 1%
# DEL
differential.missingness(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 
# complete DUP
differential.missingness(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 
# DEL+partial DUP
differential.missingness(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 


