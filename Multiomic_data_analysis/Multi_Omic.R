
## ----load_data-----------------------------------------------------------
load("data/breast_TCGA_subset_multi_omic.RData")
summary(breast_multi)

## ----load_ade4-----------------------------------------------------------
require(ade4)
dim(breast_multi$RNAseq)

## ----pca, cache=TRUE-----------------------------------------------------
breastPCA<-dudi.pca(breast_multi$RNAseq, 
            scannf=FALSE, nf=5)

## ----load_made4----------------------------------------------------------
require(made4)

## ----pca_plot, fig.show='hide'-------------------------------------------
group<-droplevels(breast_multi$clin$ER.Status)
out <- ord(breast_multi$RNAseq, classvec=group)
plot(out, nlab=3, arraylabels=rep("T", 79))

## ----projections_plot, fig.show='hide'-----------------------------------
par(mfrow=c(2,1))
plotarrays(out$ord$co, classvec=group)
plotgenes(out, col="blue")

## ----list_genes----------------------------------------------------------
ax1 <- topgenes(out, axis=1, n=5)
ax2 <- topgenes(out, axis=2, n=5)
cbind(ax1, ax2)

## ----sPCA, cache=TRUE----------------------------------------------------
require(PMA)
dd <- t(breast_multi$RNAseq)
sout <- SPC(dd, sumabsv=3, 
            K=2, orth=TRUE)

## ----sPCA_out, size='scriptsize'-----------------------------------------
rownames(sout$u) <- rownames(dd)
rownames(sout$v) <- colnames(dd)
head(sout$u)
head(sout$v)

## ----sPCA_plot-----------------------------------------------------------
plot(sout$u, type="n", xlab="First sPCA", ylab="Second sPCA")
points(sout$u, col=as.numeric(group))

## ----sPCA_genes, size='scriptsize'---------------------------------------
ss <- sout$v[,1]
ss[ss!=0]
ax1

## ----data_cc-------------------------------------------------------------
require(CCA)
df1 <- t(breast_multi$RNAseq)[,1:1000]
df2 <- t(breast_multi$RPPA)

## ----cc, eval=FALSE------------------------------------------------------
## resCC <- cc(df1, df2)

## ----rcc, cache=TRUE-----------------------------------------------------
resRCC <- rcc(df1, df2, 0.2, 0.1) 

## ----regul_estim, eval=FALSE---------------------------------------------
## regul <- estim.regul(df1, df2)
## resRCC2 <- rcc(df1, df2, regul$lambda1, regul$lambda2)

## ----plot_rcc, fig.show='hide'-------------------------------------------
plt.cc(resRCC)

## ----plot_rcc_out, fig.width=12, echo=FALSE------------------------------
plt.cc(resRCC)

## ----multiCCA------------------------------------------------------------
require(PMA)
ddlist <- list(df1, df2)
perm.out <- MultiCCA.permute(ddlist, 
                             type=c("standard", "standard"),
                             trace=FALSE) 
resMultiCCA <- MultiCCA(ddlist,  
                        penalty=perm.out$bestpenalties, 
                        ws=perm.out$ws.init, 
                        type=c("standard", "standard"), 
                        ncomponents=1, trace=FALSE, standardize=TRUE)

## ----multiCCA_out, size='scriptsize'-------------------------------------
rownames(resMultiCCA$ws[[1]]) <- colnames(df1)
rownames(resMultiCCA$ws[[2]]) <- colnames(df2)
head(resMultiCCA$ws[[1]])
head(resMultiCCA$ws[[2]])

## ----cia, cache=TRUE-----------------------------------------------------
resCIA <- cia(breast_multi$RNAseq, breast_multi$RPPA)

## ----plot_cia, fig.show='hide'-------------------------------------------
plot(resCIA, classvec=group, nlab=3, clab=0, cpoint=3 )

## ----mcia, fig.show='hide', cache=TRUE-----------------------------------
require(omicade4)
resMCIA <- mcia( breast_multi[ c(1,3,4,5,6,7) ] )

## ----plot_mcia, fig.show='hide'------------------------------------------
plot(resMCIA, axes=1:2, sample.lab=FALSE, sample.legend=FALSE, 
     phenovec=group, gene.nlab=2, 
     df.color=c("cyan", "magenta", "red4", "brown","yellow", "orange"),
     df.pch=2:7)

## ----plot_eigen, fig.show='hide'-----------------------------------------
plot(resMCIA$mcoa$cov2,  xlab = "pseudoeig 1", 
     ylab = "pseudoeig 2", pch=19, col="red")
text(resMCIA$mcoa$cov2, labels=rownames(resMCIA$mcoa$cov2), 
     cex=1.4, adj=0)

