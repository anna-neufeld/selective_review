library(dplyr)
library(ggplot2)
library(DuoClustering2018) 
library(scry)
library(glmpca)

set.seed(1)

### Load the full Zheng dataset.
sce <- sce_full_Zhengmix4eq()
X <- t(counts(sce)) 
dim(X)

### Subset to one true cell type.
### Negative control. 
cm <- as.data.frame(colData(sce)) 
true_cell_types <- cm$phenoid
table(true_cell_types)
nullX <- X[true_cell_types=="cd14.monocytes",]
nullcm <- cm[true_cell_types=="cd14.monocytes",]

### Standard preprocessing.
### Don't want cells/genes with not that many counts.
keepGenes <- which(colSums(nullX) > 20)
keepCells <- which(rowSums(nullX) > 10)
umi <- nullX[keepCells, keepGenes]
dim(umi)
nullcm <- nullcm[keepCells,]
sizefacs <- rowSums(umi)

### Let's use the same 1000 genes for every method.
### Let's sort by their estimated biological signal.
### Let's not take the top 1000, because we want negative control. The top 1000 might actually pick up cell cycle or biological signal.
### But let's not take the bottom 1000, because these genes are just all 0s. 
### Let's take numbers 1000 to 1999. 
middle_1000_genes <-order(scry::devianceFeatureSelection(t(umi), fam="poisson"),
                          decreasing=T)[1000:1999]
reduced_dim_x <- umi[,middle_1000_genes]
sce <- SingleCellExperiment(assays=list(counts=t(reduced_dim_x)), colData=nullcm)

#### Classical method. 
sce <- GLMPCA(sce, 30, assay="counts") 
set.seed(1)
clust.naive <- as.factor(kmeans(reducedDims(sce)$GLMPCA, centers=2, nstart=20)$cluster)
null.p.sizefac <- apply(reduced_dim_x, 2, function(u) summary(glm(u~clust.naive+offset(log(sizefacs)), family="poisson"))$coefficients[2,4])

#### Full conditional selective inference.  
#### Tried two 
sce <- scater::logNormCounts(sce, transform="log")
logcounts <- t(logcounts(sce))
sig.diag <- diag(apply(logcounts,2,var))
library(CADET)
ngenes <- NCOL(logcounts)
cadetres_indep <- rep(NA, ngenes)
for (i in 1:ngenes) {
  print(i)
  res_indep <- kmeans_inference_1f(logcounts,k=2,cluster_1=1,cluster_2=2,feat=i,covMat=sig.diag, seed=1, iter.max=20)
  cadetres_indep[i] <- res_indep$pval
}

### Poisson data thinning.
set.seed(10101)
library(countsplit)
countsplit.poisson <- countsplit(reduced_dim_x)
sce.pois.train <- SingleCellExperiment(assays=list(counts=t(countsplit.poisson[[1]])), colData=nullcm)
sce.pois.train <- GLMPCA(sce.pois.train, 30, assay="counts")
clust.pcs <- as.factor(kmeans(reducedDims(sce.pois.train)$GLMPCA, centers=2, nstart=20)$cluster)
pcs.p.sizefac <-  apply(countsplit.poisson[[2]], 2, function(u) summary(glm(u~clust.pcs+offset(log(sizefacs)), family="poisson"))$coefficients[2,4])


### Negative binomial data thinning. 
### Requres first estimating the overdispersion.
set.seed(2)
vst_out <- sctransform::vst(t(reduced_dim_x), n_genes = NULL)
overdisps_est <- rep(Inf, ncol(reduced_dim_x))
names(overdisps_est) <-  colnames(reduced_dim_x)
overdisps_est[rownames(vst_out$model_pars)] <- vst_out$model_pars[,1] 
set.seed(1)
### Now ready to thin. 
countsplit.nb <- countsplit(reduced_dim_x, overdisps=overdisps_est)

### Now ready to do analysis
sce.nb.train <- SingleCellExperiment(assays=list(counts=t(countsplit.nb[[1]])), colData=nullcm)
set.seed(3)
sce.nb.train <- GLMPCA(sce.nb.train, 30, assay="counts", fam="nb2") 
set.seed(521)
clust.nbcs <- as.factor(kmeans(reducedDims(sce.nb.train)$GLMPCA, centers=2, nstart=20)$cluster)
table(clust.nbcs)
nbcs.p.sizefacs <-  apply(countsplit.nb[[2]], 2, function(u) try(summary(MASS::glm.nb(u~clust.nbcs+offset(log(sizefacs))))$coefficients[2,4]))

source("fission_inf.R")
fission.p.sizefac <- rep(NA, 1000)
for (u in 1:1000) {
  if (u%%100==1) {print(u)}
  X1 <- countsplit.poisson[[1]][,u]
  X2 <- countsplit.poisson[[2]][,u]
  pval <- try(fission_test_sizefac(X1,X2,clust.pcs, 0.5, sizefacs))
  fission.p.sizefac[u] <- pval
}


save.image("inf_results_march30_1000.RData")

### Load from here if you just want to remake plot, without slow code!
#load("inf_results_march30.RData")

ggplot(data=NULL)+
  geom_qq(aes(sample=null.p.sizefac, col="aNaive"), distribution=stats::qunif)+
  geom_qq(aes(sample=as.numeric(pcs.p.sizefac), col="bPCS"), distribution=stats::qunif)+
  geom_qq(aes(sample=as.numeric(nbcs.p.sizefacs), col="cNBCS"), distribution=stats::qunif)+
  geom_qq(aes(sample=as.numeric(fission.p.sizefac), col="dFission"), distribution=stats::qunif)+
  geom_qq(aes(sample=as.numeric(cadetres_indep), col="zCadet_seed1"), distribution=stats::qunif)+
  geom_abline()+coord_fixed()+labs(col="Method")+theme_bw()+
  scale_color_discrete(
    name = "Method",
    labels = c("Classical", "Poisson thinning",
               "NB Thinning", "NB Fission", "Full CSI",
               "FullCSI2"))+
  ylab("Observed p-values")+
  xlab("Uniform(0,1) Quantiles")+theme(legend.position="bottom",
                                       legend.title=element_blank(),
                                       legend.text = element_text(size = 14))+
  guides(color=guide_legend(nrow=2, byrow=T))
ggsave("zheng_null_1000.png", width=4, height=4.5)

mean(null.p.sizefac < 0.05)
mean(pcs.p.sizefac < 0.05)
mean(as.numeric(nbcs.p.sizefacs) < 0.05)
mean(fission.p.sizefac < 0.05)
mean(cadetres_indep< 0.05)


