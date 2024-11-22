setwd("SCRIPTS/00.Importation/")
files.sources <- list.files()
sapply(files.sources, source)
setwd("../../")

Matrices <- import.matrice(type = "dosage", freq = 0.01, type.recouvrement = "reciproque")

dat <- Matrices$DosageMatrice

res <- table(dat$apoe_genotype, dat$NM_001130867.2)


res <- table(dat$apoe_genotype, dat$NM_001130867.2, dat$Status)

# CTRL N=9466 sans NA sur APOE 
# ,  = Control
#       1    2    3
# 22    0   78    0
# 23    0 1464    2
# 24    0  182    1
# 33    0 5986    5
# 34    0 1637   10
# 44    0  100    1
# NA    0   92    1

#       nc    dup
# 0 e4  7528    7
# 1 e4  1819   11
# 2 e4   100    1

m <- matrix(c(7528,1819,100,7,11,1), ncol = 2)
chisq.test(m)$expected
chisq.test(m)
fisher.test(m)



# LOAD N=8448 sans NA sur APOE 
# ,  = LOAD
#       1    2    3
# 22    0   30    0
# 23    0  526    1
# 24    0  191    1
# 33    1 4072   10
# 34    1 3203   22
# 44    0  384    6
# NA    0    9    1


#       nc    dup
# 0 e4  4629   11
# 1 e4  3394   23
# 2 e4   384    6

m <- matrix(c(4629,3395,384,11,23,6), ncol = 2)
chisq.test(m)$expected
chisq.test(m)
fisher.test(m)




# EOAD  N=4071 sans NA sur APOE 
# , ,  = EOAD
#       1    2    3
# 22    0    3    0
# 23    0  137    0
# 24    0   78    2
# 33    0 1560    4
# 34    1 1604   14
# 44    0  660    8
# NA    0    6    0

#       nc    dup
# 0 e4  1700    4
# 1 e4  1683   16
# 2 e4   660    8

m <- matrix(c(1700,1683,660,4,16,8), ncol = 2)
chisq.test(m)$expected
chisq.test(m)
fisher.test(m)


