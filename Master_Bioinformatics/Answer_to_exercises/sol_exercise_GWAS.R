#
# LOAD LIBRARY 
#

library(snpStats)

#
# LOAD DATA
#

# genotypes
setwd("Data_for_exercises")
plink <- read.plink("coronary")
geno <- plink$genotypes
geno

# annotation
annotation <- plink$map

#phenotype
feno <- read.delim("coronary.txt", sep="")
rownames(feno) <- feno$id



#
# CHECK ORDER OF INDIVIDUALS
#


identical (rownames(feno), rownames(geno))

sel <- intersect(rownames(feno), rownames(geno))

geno <- geno[sel,]
feno <- feno[sel, ]


identical (rownames(feno), rownames(geno))

geno

#
# QC SNPs
#

info.snps <- col.summary(geno)
head(info.snps)

use <- info.snps$Call.rate > 0.95 &
       info.snps$MAF > 0.05 &
       abs(info.snps$z.HWE < 3.3)    
mask.snps <- use & !is.na(use)

geno.qc.snps <- geno[ , mask.snps]
geno.qc.snps

annotation <- annotation[mask.snps, ]



#
# ASSOCIATION ANALYSIS
#

res <- snp.rhs.tests(bmi ~ 1, data=feno, snp.data=geno.qc.snps,
                     family="Gaussian")

res[1:5,]

bonf.sig <- 1e-7
ps <- p.value(res)
res[ps < bonf.sig & !is.na(ps), ]

#
#  QQ-plot
#

chi2 <- chi.squared(res) 
qq.chisq(chi2)



#
# ASSOCIATION ANALYSIS - ADJUSTED ANALYSIS (population stratification)
#

res.adj <- snp.rhs.tests(bmi ~ ev3 + ev4, family="Gaussian",
                         data=feno, snp.data=geno.qc.snps)


chi2.adj <- chi.squared(res.adj)
qq.chisq(chi2.adj)

