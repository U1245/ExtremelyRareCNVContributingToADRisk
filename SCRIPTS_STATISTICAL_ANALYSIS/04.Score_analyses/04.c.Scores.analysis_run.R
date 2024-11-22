
source("SCRIPTS/04.Score_analyses/04.c.Scores.analysis.R")

### pLI
# recouvrement réciproque, freq 1%
# DEL
score_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "pLI", p.value_threshold = 1e-5)
# complete DUP
score_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "pLI", p.value_threshold = 1e-5)
# DEL+partial DUP
score_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "pLI", p.value_threshold = 1e-5)

### pHaplo/pTriplo
# recouvrement réciproque, freq 1%
# DEL
score_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "cell", p.value_threshold = 1e-5)
# complete DUP
score_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "cell", p.value_threshold = 1e-5)
# DEL+partial DUP
score_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "cell", p.value_threshold = 1e-5)

### DS
# recouvrement réciproque, freq 1%
# DEL
score_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "DS", p.value_threshold = 1e-5)
# complete DUP
score_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "DS", p.value_threshold = 1e-5)
# DEL+partial DUP
score_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "DS", p.value_threshold = 1e-5)

### GTEx
# recouvrement réciproque, freq 1%
# DEL
score_genes.analysis(type.CNV = "DEL", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "GTEx", p.value_threshold = 1e-5)
# complete DUP
score_genes.analysis(type.CNV = "DUP", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "GTEx", p.value_threshold = 1e-5)
# DEL+partial DUP
score_genes.analysis(type.CNV = "DEL_and_DUPpartial", filtre.freq.value = 0.01, filtre.freq.type.recouvrement = "reciproque", score = "GTEx", p.value_threshold = 1e-5)


