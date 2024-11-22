
# Load dosage matrice
load("/storage/store-04/Save/Neuro/ADES_ADSP_CNV/AnalyseFinal/MatriceResultat/matriceDosage.Rdata")

# Load PCA information
load("../clinData.Rdata")

# Data with PCA
keptsample <- subset(clinData, !is.na(PC1) & exclusion == 0)
dim(keptsample)

# Add carriers information for genes of interest to PCA info data
dosage.info <- subset(DosageMatrice, select = c(sample, NM_000242.3, NM_178128.6, NM_001316767.2, NM_153334.7, NM_000185.4))
# MBL2 = NM_000242.3, FADS6 = NM_178128.6, HJV = NM_001316767.2
# SCARF2 = NM_153334.7, SERPIND1 = NM_000185.4

dat <- merge(keptsample, dosage.info, by = "sample")
dim(dat)

# Project carriers on PCA
library(ggplot2)
library(gridExtra)
p.MBL2 <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_000242.3 != 2), aes(x=PC1, y=PC2, color=as.factor(NM_000242.3), shape = Status)) +
  scale_color_manual(values = c("orange", "red", "blue")) + theme_minimal() +
  labs(color = "MBL2 dosage")

p.MBL2bis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_000242.3 != 2), aes(x=PC3, y=PC4, color=as.factor(NM_000242.3), shape = Status)) +
  scale_color_manual(values = c("orange", "red", "blue")) + theme_minimal() +
  labs(color = "MBL2 dosage")

p.FADS6 <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_178128.6 != 2), aes(x=PC1, y=PC2, color=as.factor(NM_178128.6), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "FADS6 dosage")

p.FADS6bis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_178128.6 != 2), aes(x=PC3, y=PC4, color=as.factor(NM_178128.6), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "FADS6 dosage")

p.HJV <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_001316767.2 != 2), aes(x=PC1, y=PC2, color=as.factor(NM_001316767.2), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "HJV dosage")

p.HJVbis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_001316767.2 != 2), aes(x=PC3, y=PC4, color=as.factor(NM_001316767.2), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "HJV dosage")


p.SCARF2 <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_153334.7 != 2), aes(x=PC1, y=PC2, color=as.factor(NM_153334.7), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "SCARF2 dosage")

p.SCARF2bis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_153334.7 != 2), aes(x=PC3, y=PC4, color=as.factor(NM_153334.7), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "SCARF2 dosage")

p.SERPIND1 <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_153334.7 != 2), aes(x=PC1, y=PC2, color=as.factor(NM_153334.7), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "SERPIND1 dosage")

p.SERPIND1bis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_000185.4 != 2), aes(x=PC3, y=PC4, color=as.factor(NM_000185.4), shape = Status)) +
  scale_color_manual(values = c("red", "blue")) + theme_minimal() +
  labs(color = "SERPIND1 dosage")

p.info.datasetid <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status, color = datasetid)) + geom_point(alpha=0.2) +
 theme_minimal() + theme(legend.position = "none")
p.info.datasetidbis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status, color = datasetid)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")

p.info.site <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")
p.info.sitebis <- ggplot(dat, aes(x=PC3, y=PC4, shape = Status, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")

p.info.eoad <- ggplot(subset(dat, Status == "EOAD"), aes(x=PC1, y=PC2, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none") + xlim(-25,75) + ylim(-30,30)
p.info.eoadbis <- ggplot(subset(dat, Status == "EOAD"), aes(x=PC3, y=PC4, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")

p.info.load <- ggplot(subset(dat, Status == "LOAD"), aes(x=PC1, y=PC2, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none") + xlim(-25,75) + ylim(-30,30)
p.info.loadbis <- ggplot(subset(dat, Status == "LOAD"), aes(x=PC3, y=PC4, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")

p.info.ctrl <- ggplot(subset(dat, Status == "Control"), aes(x=PC1, y=PC2, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none") + xlim(-25,75) + ylim(-30,30)
p.info.ctrlbis <- ggplot(subset(dat, Status == "Control"), aes(x=PC3, y=PC4, color = site_nice)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")

p.info.status <- ggplot(dat, aes(x=PC1, y=PC2, color = Status)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")
p.info.statusbis <- ggplot(dat, aes(x=PC3, y=PC4, color = Status)) + geom_point(alpha=0.2) +
  theme_minimal() + theme(legend.position = "none")

pdf("RESULTS/Figures/PC_analysis.pdf", height = 6, width = 12)
grid.arrange(p.MBL2, p.MBL2bis, ncol = 2, nrow = 1)
grid.arrange(p.FADS6, p.FADS6bis, ncol = 2, nrow = 1)
grid.arrange(p.HJV, p.HJVbis, ncol = 2, nrow = 1)
grid.arrange(p.SCARF2, p.SCARF2bis, ncol = 2, nrow = 1)
grid.arrange(p.info.datasetid, p.info.datasetidbis, ncol = 2, nrow = 1)
grid.arrange(p.info.site, p.info.sitebis, ncol = 2, nrow = 1)
grid.arrange(p.info.eoad, p.info.eoadbis, ncol = 2, nrow = 1)
grid.arrange(p.info.load, p.info.loadbis, ncol = 2, nrow = 1)
grid.arrange(p.info.ctrl, p.info.ctrlbis, ncol = 2, nrow = 1)
grid.arrange(p.info.status, p.info.statusbis, ncol = 2, nrow = 1)
# grid.arrange(p.SERPIND1, p.SERPIND1bis, ncol = 2, nrow = 1)
dev.off()

dat$MBL2 <- as.factor(dat$NM_000242.3)
p.MBL2 <- ggplot(dat, aes(x=PC1, y=PC2, shape = Status)) + geom_point(color = "gray") +
  geom_point(data=subset(dat, NM_000242.3 != 2), aes(x=PC1, y=PC2, color=MBL2, shape = Status)) +
  scale_color_manual(values = c("0"="orange", "1"="red", "2"="gray", "3"="blue"), labels = c("0 (homozygote deletion)", "1 (heterozygote deletion)", "2 copies", "3 (duplication)")) + theme_minimal() +
  labs(color = "MBL2 dosage") 
 

pdf("RESULTS/Figures/PC_analysis_MBL2.pdf", height = 6, width = 7)
print(p.MBL2)
# grid.arrange(p.SERPIND1, p.SERPIND1bis, ncol = 2, nrow = 1)
dev.off()

jpeg("RESULTS/Figures/PC_analysis_MBL2.jpeg", height = 600, width = 700)
print(p.MBL2)
# grid.arrange(p.SERPIND1, p.SERPIND1bis, ncol = 2, nrow = 1)
dev.off()

# proportion dans le petit nuage en bas à droite de early versus late
n1 <- subset(dat, PC1>43.75 & PC2< -5)
table(n1$Status)
# 35/(35+656)
# idem dans le gros nuage
n2 <- subset(dat, (PC1<37.5 & PC2>5) | (PC1<12.5))
table(n2$Status)
# 4024/(8767+4024)
