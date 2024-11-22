setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")


library("readr")
# list of genes
list.gene.del_CNV <- c("PRAMEF26", "ADI1", "MTUS1", "MBL2", "FADS6", "PLIN4", "KLC3", "ERCC2", "ZNF74", "SCARF2", "KLHL22", "MED15", "PI4KA", "SERPIND1", "SNAP29", "CRKL", "LZTR1", "THAP7", "P2RX6", "SLC7A4", "LRRC74B")
list.gene.SNV <- c("SORL1", "TREM2", "ABCA7", "ATP8B4", "ABCA1", "ADAM10", "CBX3", "B3GNT4", "SRC", "MOT1", "SLC16A1", "PRSS3")
list.gene.GWAS <- as.character(read.table("DATA/SCORES_AND_LIST_OF_GENES/list_genes_GWAS_Bellenguez.txt")$V1)
list.gene.del <- unique(c(list.gene.del_CNV, list.gene.SNV, list.gene.GWAS))

# load other useful information
load("../GeneTranscript.Rdata")
load("../clinData.Rdata")

B  <- length(list.gene.del)

resultsB <- foreach(iter = 1:B, .combine = "rbind") %do% {
  gene.name <- list.gene.del[iter]
  
  system(paste0('grep "n_variant" ../LOF_DATA/LOG/log', gene.name, '.txt | sed -e "s/;/\\n/g" | sed -e "s/.*=//" > ../LOF_DATA/LOG/n_variant_', gene.name, '.txt'))
  my_file <- try(readLines(paste0("../LOF_DATA/LOG/n_variant_", gene.name, ".txt")))
  
  if (length(my_file) != 0) {
    snv <- read.table(paste0("../LOF_DATA/RESULT/", gene.name, "_variants.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # position
    snv$position <- sapply(strsplit(snv$variant_id, split = ":"), "[", 2)
    
    pp <- try(readLines(paste0("../LOF_DATA/RESULT/varToExclude.", gene.name, ".tsv")))
    if (length(pp) == 0) {
      position.excluded <- c()
    } else {
      position.excluded <- unique(as.character(read.table(paste0("../LOF_DATA/RESULT/varToExclude.", gene.name, ".tsv"), stringsAsFactors = FALSE)$V2))
    }
    # extract only snv of interest
    snv.interest <- subset(snv, (is.na(exclude_reason) | exclude_reason == "") & !(position %in% position.excluded))
    if (dim(snv.interest)[1] == 0) {
      snv.interest[1,]=NA  # ad a temporary new row of NA values
      snv.interest[,'exclude'] = NA # adding new column, called for example 'new_column'
      snv.interest = snv.interest[0,] 
    } else {
      snv.interest$exclude <- "in analysis"
    }
      snv.excluded <- subset(snv, (!is.na(exclude_reason) & !exclude_reason == "") |  position %in% position.excluded)
    if (dim(snv.excluded)[1] == 0) {
      snv.excluded[1,]=NA  # ad a temporary new row of NA values
      snv.excluded[,'exclude'] = NA # adding new column, called for example 'new_column'
      snv.excluded = snv.excluded[0,] 
    } else {
      snv.excluded$exclude <- "excluded"
      snv.excluded$exclude_reason <- ifelse(snv.excluded$position %in% position.excluded, paste0( snv.excluded$exclude_reason, "; position excluded by Olivier"), snv.excluded$exclude_reason)
    }
    dat_res <- rbind(snv.interest, snv.excluded)
    dat_res$Gene <- gene.name
    return(dat_res)
  } else {
    snv <- read.table(paste0("../LOF_DATA/RESULT/ADI1_variants.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    dat_res <- snv[1,]
    dat_res[1,] <- rep(NA, 18)
    dat_res$exclude <- NA
    dat_res$Gene <- NA
    dat_res[1,"Gene"] <- paste0(gene.name, " (No SNV)")
    return(dat_res)
  }
}


results_simplified <- subset(resultsB, select = c(Gene, variant_id, position, exclude, exclude_reason, cases, controls, annotation, annotation_impact, lof_info))
write.table(results_simplified, "DATA/SNV_LOF_Table.txt", col.names = TRUE, row.names = FALSE, sep = "\t")



