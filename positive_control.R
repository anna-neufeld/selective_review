library(dplyr)
library(ggplot2)
library(DuoClustering2018) 
library(scry)
library(glmpca)

set.seed(1)

### Load Zheng data
sce <- sce_full_Zhengmix4eq()
X <- t(counts(sce)) 
dim(X)

### For positive control, keep two true cell types.
cm <- as.data.frame(colData(sce)) 
true_cell_types <- cm$phenoid
table(true_cell_types)
altX <- X[true_cell_types=="cd14.monocytes" | true_cell_types=="b.cells",]
altcm <- cm[true_cell_types=="cd14.monocytes"| true_cell_types=="b.cells",]
true_cell_types <- true_cell_types[true_cell_types=="cd14.monocytes" | true_cell_types=="b.cells"]

### Standard preprocessing
keepGenes <- which(colSums(altX) > 20)
keepCells <- which(rowSums(altX) > 10)
umi <- altX[keepCells, keepGenes]
dim(umi)
altcm <- altcm[keepCells,]
sizefacs <- rowSums(umi)

### Use a common set of 1000 genes for all methods.
### In this case, the top 1000! 
sce <- SingleCellExperiment(assays=list(counts=t(umi)), colData=altcm)
sce <- devianceFeatureSelection(sce, assay="counts", sorted=T, nkeep=1000)
reduced_dim_x <- t(counts(sce))

#### Classical method
sce <- GLMPCA(sce, 30, assay="counts") 
set.seed(1)
clust.naive <- as.factor(kmeans(reducedDims(sce)$GLMPCA, centers=2, nstart=20)$cluster)
alt.p.sizefac <- apply(reduced_dim_x, 2, function(u) summary(glm(u~clust.naive+offset(log(sizefacs)), family="poisson"))$coefficients[2,4])
sum(alt.p.sizefac < 0.01)


#### Full CSI
sce <- scater::logNormCounts(sce, transform="log")
logcounts <- t(logcounts(sce))

### Need an estimate of Sigma. This one takes into account cluster means. 
cluster.preCadet.2 <- kmeans_estimation(logcounts, 2, iter.max=20, seed=121)
sig_est_indep <- diag(apply(logcounts,2,function(u) 
  summary(lm(u~as.factor(cluster.preCadet.2$final_cluster)))$sigma))


library(CADET)
ngenes <- NCOL(logcounts)
genes <- 1:NCOL(logcounts)
cadetres <- rep(NA, ngenes)
for (i in 1:ngenes) {
  print(i)
  res <- kmeans_inference_1f(logcounts,k=2,
    cluster_1=1,cluster_2=2,feat=genes[i],covMat=sig_est_indep, seed=121, iter.max=20)
  cadetres[i] <- res$pval
}
cadet_clust <- res$final_cluster
names(cadetres) <- colnames(logcounts)

### Results for the classical method and full CSI
mclust::adjustedRandIndex(clust.naive, true_cell_types)
mclust::adjustedRandIndex(cadet_clust, true_cell_types)
naive.imp <- names(alt.p.sizefac)[alt.p.sizefac < 0.01]
cadet.imp <- names(cadetres)[cadetres < 0.01]
length(naive.imp)
length(cadet.imp)
length(intersect(naive.imp, cadet.imp))


### Poisson data thinning
set.seed(1)
library(countsplit)
countsplit.poisson <- countsplit(reduced_dim_x)
sce.pois.train <- SingleCellExperiment(assays=list(counts=t(countsplit.poisson[[1]])), colData=altcm)
sce.pois.train <- GLMPCA(sce.pois.train, 30, assay="counts")
clust.pcs <- as.factor(kmeans(reducedDims(sce.pois.train)$GLMPCA, centers=2, nstart=20)$cluster)
selgenes <- rownames(counts(sce.pois.train))
pcs.p.sizefac <-  apply(countsplit.poisson[[2]][,selgenes], 2, function(u) summary(glm(u~clust.pcs+offset(log(sizefacs)), family="poisson"))$coefficients[2,4])


### Poisson data thinnign results
mclust::adjustedRandIndex(clust.pcs, true_cell_types)
naive.imp <- names(alt.p.sizefac)[alt.p.sizefac < 0.01]
pois.imp <- names(pcs.p.sizefac)[pcs.p.sizefac < 0.01]
length(pois.imp)
length(intersect(naive.imp, pois.imp))


### Negative binomial data thinning
### Preprocessing results
set.seed(2)
vst_out <- sctransform::vst(t(reduced_dim_x), n_genes = NULL)
overdisps_est <- rep(Inf, ncol(reduced_dim_x))
names(overdisps_est) <-  colnames(reduced_dim_x)
overdisps_est[rownames(vst_out$model_pars)] <- vst_out$model_pars[,1] 
countsplit.nb <- countsplit(reduced_dim_x, overdisps=overdisps_est)

sce.nb.train <- SingleCellExperiment(assays=list(counts=t(countsplit.nb[[1]])), colData=altcm)
sce.nb.train <- GLMPCA(sce.nb.train, 50, assay="counts") 
clust.nbcs <- as.factor(kmeans(reducedDims(sce.nb.train)$GLMPCA, centers=2, nstart=20)$cluster)
selgenes <- rownames(counts(sce.nb.train))
nbcs.p.sizefacs <-  apply(countsplit.nb[[2]][,selgenes], 2, function(u) try(summary(MASS::glm.nb(u~clust.nbcs+offset(log(sizefacs))))$coefficients[2,4]))

### Negative binomial results
mclust::adjustedRandIndex(clust.nbcs, true_cell_types)
nbcs.imp <- names(nbcs.p.sizefacs)[nbcs.p.sizefacs < 0.01]
length(naive.imp)
length(nbcs.imp)
length(intersect(naive.imp, nbcs.imp))

### Fission
source("fission_inf.R")
fission.p.sizefac <- rep(NA, 1000)
selgenes.poisson <- rownames(counts(sce.pois.train))
for (u in 1:1000) {
  gene <- selgenes.poisson[u]
  X1 <- countsplit.poisson[[1]][,gene]
  X2 <- countsplit.poisson[[2]][,gene]
  pval <- try(fission_test_sizefac(X1,X2,clust.pcs, 0.5, sizefacs))
  fission.p.sizefac[u] <- pval
  print(pval)
}

names(fission.p.sizefac) <- selgenes.poisson

### Fission results
fission <- names(fission.p.sizefac)[fission.p.sizefac < 0.01]
length(naive.imp)
length(fission)
length(intersect(naive.imp, fission))



