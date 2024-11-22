
source("SCRIPTS/02.Transcript_analyses/02.c.Transcript.analysis_PC.R")
# recouvrement réciproque, freq 1%
# DEL
transcript.analysis.PC(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 
# complete DUP
transcript.analysis.PC(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 
