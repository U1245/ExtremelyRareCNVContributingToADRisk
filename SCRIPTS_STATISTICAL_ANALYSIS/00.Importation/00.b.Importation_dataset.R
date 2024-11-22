############################################# #
## Load dataset
############################################# #

# 2 dossiers de matrices
# MatriceResultat et MatriceResultatFreqModifie
# qui contiennent chacun 2 types de Rdata (matrices "standard" vs matrice "dosage") pour un filtre de frequence à 1% et à 10%
# MatriceResultat => le filtre de fréquence est appliqué bilatéraliement
# MatriceResultatFreqModifie => le filtre de fréquence est appliqué unilatéralement
# matrice "standard" : DupFullMatrice"     "DupPartialMatrice"     "DelMatrice"     "DelMatriceWithPartial"
# matrice "dosage" : "DosageMatriceWithPartial"     "DosageMatrice"

import.matrice <- function(type, freq, type.recouvrement) {
  
  if (type.recouvrement == "reciproque") {
    if (type == "CNV") {
      if (freq == 0.01) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultat/matrice.Rdata")
      }
      if (freq == 0.1) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultat/matrice10pct.Rdata")  
      }
    }
    if (type == "dosage") {
      if (freq == 0.01) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultat/matriceDosage.Rdata")
      }
      if (freq == 0.1) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultat/matriceDosage10pct.Rdata")  
      }
    }
  }
  
  if (type.recouvrement == "unilateral") {
    if (type == "CNV") {
      if (freq == 0.01) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultatFreqModifie/matrice.Rdata")
      }
      if (freq == 0.1) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultatFreqModifie/matrice10pct.Rdata")  
      }
    }
    if (type == "dosage") {
      if (freq == 0.01) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultatFreqModifie/matriceDosage.Rdata")
      }
      if (freq == 0.1) {
        Object <- load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultatFreqModifie/matriceDosage10pct.Rdata")  
      }
    }
  }
  return(mget(Object))
}

# nb columns that are not transcripts
nbNonGenes <- 13

