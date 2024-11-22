
library("readr")

load("../MatriceResultat/matrice.Rdata")
filtre.freq.value = 0.01; filtre.freq.type.recouvrement = "reciproque"; LOAD90 = FALSE

### Extract gene of interest from diffrent sources
# CNV analysis
res1.del <- read.table(paste0("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_", ifelse(isTRUE(LOAD90), "LOADinf90_", ""), "and_subset_analyses_", "DOSAGE_DEL_DUP", "_filtre_freq_", filtre.freq.value, "_recouvrement_", filtre.freq.type.recouvrement, ".txt"), sep = "\t", header = TRUE)
list.gene.del_CNV <- c("PRAMEF26", "ADI1", "MTUS1", "MBL2", "FADS6", "PLIN4", "KLC3", "ERCC2", "ZNF74", "SCARF2", "KLHL22", "MED15", "PI4KA", "SERPIND1", "SNAP29", "CRKL", "LZTR1", "THAP7", "P2RX6", "SLC7A4", "LRRC74B")
# Holstege et al exome
list.gene.SNV <- c("SORL1", "TREM2", "ABCA7", "ATP8B4", "ABCA1", "ADAM10", "CBX3", "B3GNT4", "SRC", "MOT1", "SLC16A1", "PRSS3", "TYROBP")
# GWAS Holstege et al
list.gene.GWAS <- as.character(read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt")$V1)
# Combine all genes
list.gene.del <- unique(c(list.gene.del_CNV, list.gene.SNV, list.gene.GWAS))

### Table of description by transcript
desc <- data.frame(Gene = as.character(), Transcript = as.character(), missingness = as.double(), nEOAD = as.double(), perEOAD = as.double(), 
                   nLOAD = as.double(), perLOAD = as.double(), nCTRL = as.double(), perCTRL = as.double(), 
                   nGeneEOAD = as.double(), nGeneLOAD = as.double(), nGeneCTRL = as.double(),
                   nGenedelEOAD = as.double(), perGenedelEOAD = as.double(), nGenedelLOAD = as.double(), perGenedelLOAD = as.double(), 
                   nGenedelCTRL = as.double(), perGenedelCTRL = as.double(), nECASCAD = as.double(), stringsAsFactors = FALSE)
k <- 1

## Sont inclus dans l'analyse que les transcripts avec <= 25% de missingness et étant étiquetés "NM"
for (gene.name in list.gene.del) {
  cat("------------------------------\n", gene.name, "\n", "------------------------------\n")
  load("../GeneTranscript.Rdata")
  load("../clinData.Rdata")
  transcript.list1 <- unique(subset(GeneTranscript, Gene == gene.name)$Transcript)
  transcript.list <- transcript.list1[substr(transcript.list1, 1, 2) == "NM"]
  if (length(transcript.list) > 0) {
    genematrix <- DelMatrice[,c("sample", "datasetid", "Status", transcript.list)]
    genematrix.vcfiid <- merge(genematrix, clinData, by = c("sample", "Status", "datasetid"), all.x = TRUE, all.y = FALSE)
    genematrix.vcfiid$status_vcf <- ifelse(genematrix.vcfiid$Status == "Control", "control", ifelse(genematrix.vcfiid$Status == "LOAD", "case>=65", ifelse(genematrix.vcfiid$Status == "EOAD", "case<65", NA)))
    
    table(genematrix.vcfiid$Status, genematrix.vcfiid$status_vcf)
    
    genematrix.vcfiid$vcfID_ok <- ifelse(genematrix.vcfiid$vcfID == ".", genematrix.vcfiid$sample, genematrix.vcfiid$vcfID)
    
    # Pour chaque gene liste des transcripts avec <= 25% de missingness
    transcript.list.miss <- c()
    for (i in transcript.list) {
      if (sum(is.na(genematrix[, i]))/(dim(genematrix)[1]) <= 0.25) {
        transcript.list.miss <- c(transcript.list.miss, i)
      }
    }
    
    genematrix.vcfiid.right.columns <- subset(genematrix.vcfiid, select = c("vcfID_ok", "Status", "datasetid", "status_vcf", transcript.list.miss))
    genematrix.vcfiid.right.columns.no.NA <- na.omit(genematrix.vcfiid.right.columns)
    
    # count by Gene
    if (length(transcript.list.miss) == 1) {
      genematrix.vcfiid.right.columns.no.NA$cnv <- apply(matrix(genematrix.vcfiid.right.columns.no.NA[, transcript.list.miss], ncol = 1), 1, sum)
    } else {
      genematrix.vcfiid.right.columns.no.NA$cnv <- apply(genematrix.vcfiid.right.columns.no.NA[, transcript.list.miss], 1, sum)
    }
    
    # selectionner colonnes à fournir à Olivier
    genematrix.vcfiid.right.columns.bis <- subset(genematrix.vcfiid.right.columns.no.NA, select = c("vcfID_ok", "status_vcf"))
    colnames(genematrix.vcfiid.right.columns.bis) <- c("iid", "status")
    
    # list of samples for Olivier
    write_tsv(genematrix.vcfiid.right.columns.bis, paste0("../LOF_DATA/SAMPLES/samples_", gene.name, ".tsv"))
    
    # Count pour description => à revoir car le N peut changer en fonctions de la fonction qu'applique Olivier pour obtenir les SNV
    for (i in transcript.list) {
      desc[k, ] <- c(gene.name, i, sum(is.na(genematrix[, i]))/(dim(genematrix)[1]), 
                     # EOAD
                     sum(genematrix[genematrix$Status == "EOAD", i] == 1), sum(genematrix[genematrix$Status == "EOAD", i] == 1)/sum(!is.na(genematrix[genematrix$Status == "EOAD", i])),
                     # LOAD
                     sum(genematrix[genematrix$Status == "LOAD", i] == 1), sum(genematrix[genematrix$Status == "LOAD", i] == 1)/sum(!is.na(genematrix[genematrix$Status == "EOAD", i])),
                     # CTRL
                     sum(genematrix[genematrix$Status == "Control", i] == 1), sum(genematrix[genematrix$Status == "Control", i] == 1)/sum(!is.na(genematrix[genematrix$Status == "Control", i])),
                     # N
                     dim(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "EOAD", ])[1], 
                     dim(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "LOAD", ])[1], 
                     dim(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "Control", ])[1],
                     # EOAD
                     sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "EOAD", "cnv"] >= 1), sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "EOAD", "cnv"] >= 1)/sum(!is.na(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "EOAD", i])),
                     # LOAD
                     sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "LOAD", "cnv"] >= 1), sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "LOAD", "cnv"] >= 1)/sum(!is.na(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "EOAD", i])),
                     # CTRL
                     sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "Control", "cnv"] >= 1), sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "Control", "cnv"] >= 1)/sum(!is.na(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$Status == "Control", i])),
                     # nECASCAD
                     sum(genematrix.vcfiid.right.columns.no.NA[genematrix.vcfiid.right.columns.no.NA$datasetid %in% c("ALZ.ECASCAD57.Agilent_V6UTR",
                                                                                                                      "ALZ.ECASCAD58.Agilent_V6UTR",
                                                                                                                      "ALZ.ECASCAD59.Agilent_V6UTR",
                                                                                                                      "ALZ.ECASCAD60_61.Agilent_V6UTR"), i] == 1)
      )
      k <- k + 1
    }
  } else {
    print("Gene not found")
  }
}
# write.table(desc, "DATA/count_deletion_by_gene_for_MA_CNV_SNV.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(desc, "DATA/count_deletion_by_gene_for_MA_CNV_SNV_addTYROBP.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
