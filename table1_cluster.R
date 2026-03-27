library(dplyr)
library(ggplot2)
library(patchwork)
library(rpart)

p <- 20
n <- 200

mu.null <- matrix(10, nrow=n, ncol=p)
mu.alt <- mu.null
true.cluster <- rep(0,n)
true.cluster[1:80] <- 1
sigma=1
beta <- 2
mu.alt[true.cluster==1,1:5] <- mu.alt[true.cluster==1,1:5] + beta

mu.strong <- mu.null
mu.strong[true.cluster==1,1:5] <- mu.strong[true.cluster==1,1:5] + 2*beta

ntrials <- 2000
alpha <- 0.1

res <- data.frame("t" = NA, "method"=NA, "setting"= NA, "rand" = NA, "cov" = NA, "width" = NA)

for (t in 1:ntrials) {
  if (t %%50 == 1) {print(t)}
  
  set.seed(t)
  
  ### No signal
  X <- apply(mu.null, 2, function(u) rnorm(length(u),u,sigma))
  est.clust.null <- kmeans(X, centers=2)$cluster

  if (min(table(est.clust.null)) > 1) {
    CI.naive  <- t.test(X[,1]~ est.clust.null, conf.level=1-alpha)$conf.int
    res <- res %>% add_row(t=t, method="naive", setting="globalnull", 
                         rand=NA, cov = CI.naive[1] < 0  & CI.naive[2] > 0,
                         width=CI.naive[2]-CI.naive[1])

  }
  
  ### Medium signal
  X <- apply(mu.alt, 2, function(u) rnorm(length(u),u,sigma))
  est.clust.alt <- kmeans(X, centers=2)$cluster
  
  if (min(table(est.clust.alt)) > 1) {
  
    CI.naive  <- t.test(X[,1]~ est.clust.alt,conf.level=1-alpha)$conf.int
    true.mu.alt <- mean(mu.alt[est.clust.alt==1,1])-mean(mu.alt[est.clust.alt==2,1])
    res <- res %>% add_row(t=t, method="naive", setting="alt", 
                         rand=mclust::adjustedRandIndex(est.clust.alt, true.cluster), cov = CI.naive[1] < true.mu.alt  & CI.naive[2] > true.mu.alt,
                         width=CI.naive[2]-CI.naive[1])
  }
  
  ### Strong signal
  X <- apply(mu.strong, 2, function(u) rnorm(length(u),u,sigma))
  est.clust.alt <- kmeans(X, centers=2)$cluster
  
  if (min(table(est.clust.alt)) > 1) {
    
    CI.naive  <- t.test(X[,1]~ est.clust.alt,conf.level=1-alpha)$conf.int
    true.mu.alt <- mean(mu.strong[est.clust.alt==1,1])-mean(mu.strong[est.clust.alt==2,1])
    res <- res %>% add_row(t=t, method="naive", setting="strong", 
                           rand=mclust::adjustedRandIndex(est.clust.alt, true.cluster), cov = CI.naive[1] < true.mu.alt  & CI.naive[2] > true.mu.alt,
                           width=CI.naive[2]-CI.naive[1])
  }
}

res <- res[-1,]
res %>% group_by(setting, method) %>% summarize(mean(cov), mean(rand), mean(rand==1), median(width), mean(width))
