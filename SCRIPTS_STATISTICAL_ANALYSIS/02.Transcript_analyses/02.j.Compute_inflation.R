
estlambda <- function(data, plot=FALSE, proportion=1.0,
                      method="regression", filter=TRUE, df=1,... ) {
  data <- data[which(!is.na(data))]
  if (proportion>1.0 || proportion<=0)
    stop("proportion argument should be greater then zero and less than or equal to one")
  
  ntp <- round( proportion * length(data) )
  if ( ntp<1 ) stop("no valid measurements")
  if ( ntp==1 ) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate=1.0, se=999.99))
  }
  if ( ntp<10 ) warning(paste("number of points is too small:", ntp))
  if ( min(data)<0 ) stop("data argument has values <0")
  if ( max(data)<=1 ) {
    #		lt16 <- (data < 1.e-16)
    #		if (any(lt16)) {
    #			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
    #			data[lt16] <- 1.e-16
    #		}
    data <- qchisq(data, 1, lower.tail=FALSE)
  }
  if (filter)
  {
    data[which(abs(data)<1e-8)] <- NA
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  #	s <- summary(lm(data~offset(ppoi)))$coeff
  #       bug fix thanks to Franz Quehenberger
  
  out <- list()
  if (method=="regression") {
    s <- summary( lm(data~0+ppoi) )$coeff
    out$estimate <- s[1,1]
    out$se <- s[1,2]
  } else if (method=="median") {
    out$estimate <- median(data, na.rm=TRUE)/qchisq(0.5, df)
    out$se <- NA
  } else if (method=="KS") {
    limits <- c(0.5, 100)
    out$estimate <- estLambdaKS(data, limits=limits, df=df)
    if ( abs(out$estimate-limits[1])<1e-4 || abs(out$estimate-limits[2])<1e-4 )
      warning("using method='KS' lambda too close to limits, use other method")
    out$se <- NA
  } else {
    stop("'method' should be either 'regression' or 'median'!")
  }
  
  if (plot) {
    # lim <- c(0, max(data, ppoi,na.rm=TRUE))
    # #		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
    # oldmargins <- par()$mar
    # par(mar=oldmargins + 0.2)
    # plot(ppoi, data,
    #      xlab=expression("Expected " ~ chi^2),
    #      ylab=expression("Observed " ~ chi^2),
    #      ...)
    # abline(a=0, b=1)
    # abline(a=0, b=out$estimate, col="red")
    # legend("topleft", paste0("lambda = ", round(out$estimate, 2)), bty="n")
    # par(mar=oldmargins)
    # graphe personnalisé
    lambda.value = round(out$estimate, 2)
    p.qqplot <- ggplot(mapping=aes(x=ppoi, y=data)) + geom_point() + theme_minimal() +
      labs(x = expression("Expected " ~ chi^2), y = expression("Observed " ~ chi^2)) + 
      theme(axis.title=element_text(size=15)) +
      annotate(geom="text", x=15, y=3, size = 6, label=bquote(paste(lambda, " = ", .(lambda.value)))) +
      coord_cartesian(xlim = c(0,20), ylim = c(0,20)) + geom_abline(intercept = 0, slope = 1)  + geom_abline(intercept = 0, slope = out$estimate, color = "red")
    print(p.qqplot)
  }
  
  out
}




# DOSAGE 

# par(mfrow=c(1,2))

# sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")
# sub$duplicatedGene <- duplicated(sub$Gene)
# sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)
# p.EOADvsCTRL.firth <- sub.uniqueGene$p.value.firth.EOAD.vs.CTRL
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Dosage analysis : EOAD vs CTRL firth (>= 3 CNV carriers, 1 test per gene), median method")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : EOAD vs CTRL firth (>= 3 CNV carriers, 1 test per gene), regression method")
# 

# dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")
# 
# sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >= 4)
# sub$duplicatedGene <- duplicated(sub$Gene)
# sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)
# p.EOADvsCTRL.firth <- sub.uniqueGene$p.value.firth.EOAD.vs.CTRL
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Dosage analysis : EOAD vs CTRL firth (>= 4 CNV carriers, 1 test per gene), median method")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : EOAD vs CTRL firth (>= 4 CNV carriers, 1 test per gene), regression method")
# 
# sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")
# sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.EOAD)
# sub$duplicatedGene <- duplicated(sub$infoduplicated)
# sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)
# p.EOADvsCTRL.firth <- sub.uniqueGene$p.value.firth.EOAD.vs.CTRL
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Dosage analysis : EOAD vs CTRL firth (>= 3 CNV carriers, 1 test per gene x N), median method")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : EOAD vs CTRL firth (>= 3 CNV carriers, 1 test per gene x N), regression method")

# jpeg("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_main_analysis_union_setA_B_DOSAGE_DEL_DUP_recouvrement_reciproque.jpeg", width = 7, height = 7)

jpeg("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_main_analysis_union_setA_B_DOSAGE_DEL_DUP_recouvrement_reciproque.jpeg", width = 700, height = 700)

dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")

sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >= 4)
sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.EOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)
sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)
p.EOADvsCTRL.firth <- sub.uniqueGene$p.value.firth.EOAD.vs.CTRL
estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Without PC adjustment")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : EOAD vs CTRL firth (>= 4 CNV carriers, 1 test per gene x N), regression method")

dev.off()
# jpeg("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_PC_adjusted_main_analysis_union_setA_B_DOSAGE_DEL_DUP_recouvrement_reciproque.jpeg", width = 7, height = 7)

jpeg("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_PC_adjusted_main_analysis_union_setA_B_DOSAGE_DEL_DUP_recouvrement_reciproque.jpeg", width = 700, height = 700)

dat <- read.table("RESULTS/Analysis_by_transcript_union_setA_B_PC_adjusted_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")

sub <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none" & Ncarriers.CNV.total >= 4)
sub$infoduplicated <- paste0(sub$Gene, sub$N.CTRL, sub$N.EOAD, sub$Ncarriers.CNV.CTRL, sub$Ncarriers.CNV.EOAD)
sub$duplicatedGene <- duplicated(sub$infoduplicated)
sub.uniqueGene <- subset(sub, duplicatedGene == FALSE)
p.EOADvsCTRL.firth <- sub.uniqueGene$p.value.firth.EOAD.vs.CTRL
estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("With PC adjustment")

dev.off()


 
# ############################### OLD ############################################
# 
# # DUPLICATION 
# dat <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")
# p.allADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.all.AD.vs.CTRL
# p.EOADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.EOAD.vs.CTRL
# p.LOADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.LOAD.vs.CTRL
# p.ordinal.regression <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value..ordinal.regression..using.polr.
# p.ordinal.regression.clm <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value..ordinal.regression..using.clm.
# 
# pdf("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_ordinal_regression_and_subset_analyses_DUP_recouvrement_reciproque.pdf", width = 18, height = 9)
# par(mfrow=c(1,2))
# estlambda(data = p.allADvsCTRL.firth, plot = TRUE, method = "median")
# title("Duplication analysis : all AD vs CTRL firth, median method")
# estlambda(data = p.allADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Duplication analysis : all AD vs CTRL firth, regression method")
# 
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Duplication analysis : EOAD vs CTRL firth, median method")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Duplication analysis : EOAD vs CTRL firth, regression method")
# 
# estlambda(data = p.LOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Duplication analysis : LOAD vs CTRL firth, median method")
# estlambda(data = p.LOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Duplication analysis : LOAD vs CTRL firth, regression method")
# 
# estlambda(data = p.ordinal.regression, plot = TRUE, method = "median")
# title("Duplication analysis : ordinal regression, median method")
# estlambda(data = p.ordinal.regression, plot = TRUE, method = "regression")
# title("Duplication analysis : ordinal regression, regression method")
# 
# estlambda(data = p.ordinal.regression.clm, plot = TRUE, method = "median")
# title("Duplication analysis (PC) : clm ordinal regression, median method")
# estlambda(data = p.ordinal.regression.clm, plot = TRUE, method = "regression")
# title("Duplication analysis (PC) : clm ordinal regression, regression method")
# dev.off()
# 
# 
# 
# # DELETION 
# dat <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DEL_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")
# p.allADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.all.AD.vs.CTRL
# p.EOADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.EOAD.vs.CTRL
# p.LOADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.LOAD.vs.CTRL
# p.ordinal.regression <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value..ordinal.regression..using.polr.
# 
# pdf("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_ordinal_regression_and_subset_analyses_DEL_recouvrement_reciproque.pdf", width = 18, height = 9)
# par(mfrow=c(1,2))
# estlambda(data = p.allADvsCTRL.firth, plot = TRUE, method = "median")
# title("Deletion analysis : all AD vs CTRL firth, median method")
# estlambda(data = p.allADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Deletion analysis : all AD vs CTRL firth, regression method")
# 
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Deletion analysis : EOAD vs CTRL firth, median method")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Deletion analysis : EOAD vs CTRL firth, regression method")
# 
# estlambda(data = p.LOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Deletion analysis : LOAD vs CTRL firth, median method")
# estlambda(data = p.LOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Deletion analysis : LOAD vs CTRL firth, regression method")
# 
# estlambda(data = p.ordinal.regression, plot = TRUE, method = "median")
# title("Deletion analysis : ordinal regression, median method")
# estlambda(data = p.ordinal.regression, plot = TRUE, method = "regression")
# title("Deletion analysis : ordinal regression, regression method")
# dev.off()
# 
# 
# 
# 
# # DOSAGE 
# dat <- read.table("RESULTS/Analysis_by_transcript_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_filtre_freq_0.01_recouvrement_reciproque.txt", header = TRUE, sep = "\t")
# p.allADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.all.AD.vs.CTRL
# p.EOADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.EOAD.vs.CTRL
# p.LOADvsCTRL.firth <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value.firth.LOAD.vs.CTRL
# p.ordinal.regression <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value..ordinal.regression..using.polr.
# p.ordinal.regression.clm <- subset(dat, substr(transcript.name, 1, 2) == "NM" & Adjustment == "none")$p.value..ordinal.regression..using.clm.
# 
# pdf("RESULTS/Figures/inflation_NM_filtre_freq_0.01_Analysis_by_transcript_ordinal_regression_and_subset_analyses_DOSAGE_DEL_DUP_recouvrement_reciproque.pdf", width = 18, height = 9)
# par(mfrow=c(1,2))
# estlambda(data = p.allADvsCTRL.firth, plot = TRUE, method = "median")
# title("Dosage analysis : all AD vs CTRL firth, median method")
# estlambda(data = p.allADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : all AD vs CTRL firth, regression method")
# 
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Dosage analysis : EOAD vs CTRL firth, median method")
# estlambda(data = p.EOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : EOAD vs CTRL firth, regression method")
# 
# estlambda(data = p.LOADvsCTRL.firth, plot = TRUE, method = "median")
# title("Dosage analysis : LOAD vs CTRL firth, median method")
# estlambda(data = p.LOADvsCTRL.firth, plot = TRUE, method = "regression")
# title("Dosage analysis : LOAD vs CTRL firth, regression method")
# 
# estlambda(data = p.ordinal.regression, plot = TRUE, method = "median")
# title("Dosage analysis : ordinal regression, median method")
# estlambda(data = p.ordinal.regression, plot = TRUE, method = "regression")
# title("Dosage analysis : ordinal regression, regression method")
# 
# estlambda(data = p.ordinal.regression.clm, plot = TRUE, method = "median")
# title("Dosage analysis (PC) : clm ordinal regression, median method")
# estlambda(data = p.ordinal.regression.clm, plot = TRUE, method = "regression")
# title("Dosage analysis (PC) : clm ordinal regression, regression method")
# dev.off()
# 
# 
# 
