setwd("~/Dropbox/2025 talks or research/Selective Inference Review Paper/New_Section_2")
library(dplyr)
library(ggplot2)
library(patchwork)
library(rpart)

p <- 20
n <- 200
mu.null <- rep(10,n)

ntrials <- 2000
alpha <- 0.1

Sig.cor <- diag(1,p)
Sig.cor[1:5, 1:5] <- 0.9
Sig.cor[6:10, 6:10] <- 0.5
diag(Sig.cor) <- rep(4,p)
beta= 1.5 ## Change back to 1.8 if needed.
betastrong= 4

res <- data.frame("t" = NA, "method"=NA, "setting"= NA, "rand" = NA, "cov" = NA, "width" = NA)

for (t in 1:ntrials) {
  if (t %%50 == 1) {print(t)}
  
  set.seed(t)
  
  ### START WITH NULL SETTING, NAIVE
  y <- rnorm(n, mu.null, sd=1)
  X <- MASS::mvrnorm(n, mu=rep(0,p), Sigma=Sig.cor)
  dat <- data.frame(cbind(y,X))
  tree.naive <- rpart(y~., data=dat, control = rpart.control(maxdepth=1, minbucket=2), cp=-1)
  
  sel.naive = dat[tree.naive$where==3,]
  mu.sel.naive = mu.null[tree.naive$where==3]
  CI.naive <- t.test(sel.naive$y, conf.level=1-alpha)$conf.int
  
  res <- res %>% add_row(t=t, method="naive", setting="null", 
                         rand=NA, cov = CI.naive[1] < mean(mu.sel.naive)  & CI.naive[2] > mean(mu.sel.naive),
                         width=CI.naive[2]-CI.naive[1])
  
  # ### NULL SETTING, SAMPLE SPLIT
  # train <- sample(1:n, size=0.5*n)
  # tree.ss <- rpart(y~., data=dat[train,], control = rpart.control(maxdepth=1, minbucket=2), cp=-1)
  # tree.ss$frame$yval <- 1:NROW(tree.ss$frame)
  # dat.test <-dat[-train,] 
  # preds <- predict(tree.ss, newdata=dat.test)
  # sel.ss <- dat.test[preds==3,]
  # all.preds <- predict(tree.ss, newdata=dat)
  # mu.sel.ss <- mean(mu.null[all.preds==3])
  # if (NROW(sel.ss)>1) {
  # ### THIS IS ACTUALLY FROUGHT, but that is a CAN OF WORMS. 
  #   CI.ss <- t.test(sel.ss$y, conf.level=1-alpha)$conf.int
  # } else {
  #   CI.ss <- c(-Inf, Inf)
  # }
  # 
  # res <- res %>% add_row(t=t, method="sampsplit", setting="null", 
  #                        rand=NA, cov = CI.ss[1] < mean(mu.sel.ss)  & CI.ss[2] > mean(mu.sel.ss),
  #                        width=CI.ss[2]-CI.ss[1])
 
   ### STRONG
  true.group <- X[,2] > 1
  mu.alt <- mu.null
  mu.alt[true.group==1] = mu.alt[true.group==1]+beta
  y <- rnorm(n, mu.alt, sd=1)
  dat.alt <- data.frame(cbind(y,X))
  tree.alt <- rpart(y~., data=dat.alt, control = rpart.control(maxdepth=1, minbucket=2), cp=-1)
  
  sel.naive = dat.alt[tree.alt$where==3,]
  mu.sel.naive = mu.alt[tree.alt$where==3]
  CI.naive <- t.test(sel.naive$y, conf.level=1-alpha)$conf.int
  
  res <- res %>% add_row(t=t, method="naive", setting="strong", 
                         rand=mclust::adjustedRandIndex(true.group, tree.alt$where), cov = CI.naive[1] < mean(mu.sel.naive)  & CI.naive[2] > mean(mu.sel.naive),
                         width=CI.naive[2]-CI.naive[1])
  
  ### ALT SETTING, SAMPLE SPLIT
  # train <- sample(1:n, size=0.5*n)
  # tree.ss <- rpart(y~., data=dat.alt[train,], control = rpart.control(maxdepth=1, minbucket=2), cp=-1)
  # tree.ss$frame$yval <- 1:NROW(tree.ss$frame)
  # dat.test <-dat.alt[-train,] 
  # preds <- predict(tree.ss, newdata=dat.test)
  # sel.ss <- dat.test[preds==3,]
  # all.preds <- predict(tree.ss, newdata=dat.alt)
  # mu.sel.ss <- mean(mu.alt[all.preds==3])
  # if (NROW(sel.ss)>1) {
  #   ### THIS IS ACTUALLY FROUGHT, but that is a CAN OF WORMS. 
  #   CI.ss <- t.test(sel.ss$y, conf.level=1-alpha)$conf.int
  # } else {
  #   CI.ss <- c(-Inf, Inf)
  # }
  # 
  # res <- res %>% add_row(t=t, method="sampsplit", setting="strong", 
  #                        rand=mclust::adjustedRandIndex(true.group, predict(tree.ss, newdata=dat.alt)), cov = CI.ss[1] < mean(mu.sel.ss)  & CI.ss[2] > mean(mu.sel.ss),
  #                        width=CI.ss[2]-CI.ss[1])
  
  ## STRONGER
  true.group <- X[,2] > 1
  mu.alt <- mu.null
  mu.alt[true.group==1] = mu.alt[true.group==1]+betastrong
  y <- rnorm(n, mu.alt, sd=1)
  dat.alt <- data.frame(cbind(y,X))
  tree.alt <- rpart(y~., data=dat.alt, control = rpart.control(maxdepth=1, minbucket=2), cp=-1)
  
  sel.naive = dat.alt[tree.alt$where==3,]
  mu.sel.naive = mu.alt[tree.alt$where==3]
  CI.naive <- t.test(sel.naive$y, conf.level=1-alpha)$conf.int
  
  res <- res %>% add_row(t=t, method="naive", setting="stronger", 
                         rand=mclust::adjustedRandIndex(true.group, tree.alt$where), cov = CI.naive[1] < mean(mu.sel.naive)  & CI.naive[2] > mean(mu.sel.naive),
                         width=CI.naive[2]-CI.naive[1])
  
  ### ALT SETTING, SAMPLE SPLIT
  # train <- sample(1:n, size=0.5*n)
  # tree.ss <- rpart(y~., data=dat.alt[train,], control = rpart.control(maxdepth=1, minbucket=2), cp=-1)
  # tree.ss$frame$yval <- 1:NROW(tree.ss$frame)
  # dat.test <-dat.alt[-train,] 
  # preds <- predict(tree.ss, newdata=dat.test)
  # sel.ss <- dat.test[preds==3,]
  # all.preds <- predict(tree.ss, newdata=dat.alt)
  # mu.sel.ss <- mean(mu.alt[all.preds==3])
  # if (NROW(sel.ss)>1) {
  #   ### THIS IS ACTUALLY FROUGHT, but that is a CAN OF WORMS. 
  #   CI.ss <- t.test(sel.ss$y, conf.level=1-alpha)$conf.int
  # } else {
  #   CI.ss <- c(-Inf, Inf)
  # }
  # 
  # res <- res %>% add_row(t=t, method="sampsplit", setting="stronger", 
  #                        rand=mclust::adjustedRandIndex(true.group, predict(tree.ss, newdata=dat.alt)), cov = CI.ss[1] < mean(mu.sel.ss)  & CI.ss[2] > mean(mu.sel.ss),
  #                        width=CI.ss[2]-CI.ss[1])
  
}

res <- res[-1,]

# res$width_ratio = NA
# res$sel_ratio = NA
# for (t in 1:ntrials) {
#   res[res$t==t & res$method=="sampsplit",]$width_ratio =   res[res$t==t & res$method=="sampsplit",]$width/res[res$t==t & res$method=="naive",]$width
#   res[res$t==t & res$method=="sampsplit",]$sel_ratio =   res[res$t==t & res$method=="sampsplit",]$rand/res[res$t==t & res$method=="naive",]$rand
# }
# save(res, file="treevaluesRes.Rdata")
# 
# ### Reported results
# load("treevaluesRes.Rdata")

res %>% group_by(setting, method) %>% summarize(mean(cov), mean(rand==1), mean(rand), median(width), mean(width))

