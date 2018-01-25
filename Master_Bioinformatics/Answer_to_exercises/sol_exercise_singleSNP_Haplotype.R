
#
#  EXERCISE 1: single SNP analysis
#

# load library
library(SNPassoc)

# set working directory
setwd("Data_for_exercises")

# load data
data <- read.delim("DM.txt")
head(data)

# prepare SNP data
data.s <- setupSNP(data, 6:14, sep="")

# this does the same ...
ii <- grep("^rs", colnames(data))
ii <- c(ii, grep("lpr", colnames(data)))
data.s <- setupSNP(data, ii, sep="")

summary(data.s)
plotMissing(data.s)

# Association
ans <- WGassociation(RESP, data.s)
plot(ans)
ans

# get SNP name that is significant 
sel <- apply(ans, 1, function(x)  any(x < 0.05 & !is.na(x)))
sel[sel]


# ORs for genetic models
association(RESP ~ rs908867, data.s)

# Adjusted model
association(RESP ~ rs908867 + HDRS, data.s)
association(RESP ~ rs908867 + HDRS + PSICOT + MELANCOL + EPD_PREV, data.s)

# Plot for dominant, recessive and additive models
ans2 <- WGassociation(RESP, data.s, model=c("do", "re", "log"))
plot(ans2)

# Max-statistic
maxstat(data.s, RESP)



#
#  EXERCISE 2: Haplotype analysis
#

# load library

library(LDheatmap)
library(genetics)

# Read and prepare data
data <- read.delim("DM.txt")
data.s <- setupSNP(data, 6:13, sep="")

# Get SNPs as genotype objects
snpGeno <- data.frame(lapply(data.s[, labels(data.s)], genotype))
head(snpGeno)
class(snpGeno[,1])

# Get genetic position  (Internet connection is required)
library(biomaRt)  
mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

snps <- labels(data.s)
snpPos <- getBM(c("chrom_start"),  
                  filters = c("snp_filter"), 
                  values = snps, mart = mart)

# Create Heatmap
MyHeatmap <- LDheatmap(snpGeno, snpPos[,1], LDmeasure = "r",
 title = "Pairwise LD in r^2", add.map = TRUE,
 color = grey.colors(20), name = "myLDgrob", add.key = TRUE, flip=TRUE,
 SNP.name=labels(data.s) )

# Create Heatmap
MyHeatmap <- LDheatmap(snpGeno, 1:8, LDmeasure = "r",
 title = "Pairwise LD in r^2", add.map = TRUE,
 color = grey.colors(20), name = "myLDgrob", add.key = TRUE, flip=TRUE,
 SNP.name=labels(data.s) )



# Sliding window
geno <- make.geno(data.s, labels(data.s))
haplo.score <- list()
y <- as.numeric(data.s$RESP)-1  # This is required in haplo.score function
for (i in 2:4) {
   haplo.score[[i-1]] <- haplo.score.slide(y, geno, 
                          trait.type="binomial",
                          n.slide=i,
                          simulate=TRUE,
                          sim.control=score.sim.control(min.sim=100,
                                       max.sim=200)) 
 }

par(mfrow=c(2,2))
for (i in 2:4) {
   plot(haplo.score[[i-1]])
   title(paste("Sliding Window=", i, sep=""))
 }


# Haplotype estimation
snpsH <- snps[5:6]   # snpsH <- labels(datos.s)[5:6]
snpsH
genoH<-make.geno(data.s, snpsH)

em <- haplo.em(genoH, locus.label = snpsH, miss.val = c(0, NA))
em


# Haplotype association analysis
mod <- haplo.glm(y~genoH,           
             family="binomial", 
             locus.label=snpsH,
             allele.lev=attributes(genoH)$unique.alleles,
             control = haplo.glm.control(haplo.freq.min=0.05))
       

intervals(mod)

