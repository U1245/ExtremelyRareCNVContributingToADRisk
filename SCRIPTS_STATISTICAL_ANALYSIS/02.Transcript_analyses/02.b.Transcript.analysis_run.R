

source("SCRIPTS/02.Transcript_analyses/02.b.Transcript.analysis.R")
# recouvrement réciproque, freq 1%
# DEL
transcript.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 
# complete DUP
transcript.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 
# # DEL+partial DUP
transcript.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque")
# # partial DUP
# transcript.analysis(type.CNV = "DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 

# # recouvrement réciproque, freq 10%
# # DEL
# transcript.analysis(type.CNV = "DEL", filtre.freq.value = 0.1, filtre.freq.type.recouvrement = "reciproque")
# # complete DUP
# transcript.analysis(type.CNV = "DUP", filtre.freq.value = 0.1, filtre.freq.type.recouvrement = "reciproque")
# # DEL+partial DUP
# transcript.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.1, filtre.freq.type.recouvrement = "reciproque")

# 
# # recouvrement unilatéral, freq 1%
# # DEL
# transcript.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "unilateral") 
# # complete DUP
# transcript.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "unilateral") 
# # DEL+partial DUP
# transcript.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "unilateral") 
# 

# # recouvrement unilatéral, freq 10%
# # DEL
# transcript.analysis(type.CNV = "DEL", filtre.freq.value = 0.1, filtre.freq.type.recouvrement = "unilateral") 
# # complete DUP
# transcript.analysis(type.CNV = "DUP", filtre.freq.value = 0.1, filtre.freq.type.recouvrement = "unilateral") 
# # DEL+partial DUP
# transcript.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.1, filtre.freq.type.recouvrement = "unilateral") 




source("SCRIPTS/02.Transcript_analyses/02.b.bis.CountforPaper.R")
# recouvrement réciproque, freq 1%
countforpaper_dosage.analysis_union_setA_B(type.CNV = "DEL_DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque") 




