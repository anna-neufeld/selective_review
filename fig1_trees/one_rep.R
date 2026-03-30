### No pruning! Two level trees only! Etc . 
oneRepFull <-  function(n,p,sigma_y, seed, alpha = 0.1, beta=0, filename="test.txt", minbucket=1,maxdepth=2, XORlev=1) {
  
  set.seed(seed)
  
  sampsplit_eps <- c(0.5,0.7,0.9,0.95)
  datathin_eps <- sampsplit_eps
  
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  ### This will do for now, but will almost certainly want to edit! 
  mu_y <- beta*I(X[,1] < 0) +
    beta*XORlev*(I(X[,1] < 0 & X[,2] > 0)) +
    + beta*XORlev*(I(X[,1] > 0 & X[,3] > 1))
  
  y <- rnorm(n,mu_y,sigma_y)
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)
  
  ### Build rpart tree.
  base_tree <- rpart::rpart(y~., data=dat,
                            control=rpart.control(minbucket=minbucket,cp=-1, maxcompete=0,maxsurrogate=0, maxdepth=maxdepth), model=TRUE)
  
  ARI <- mclust::adjustedRandIndex(base_tree$where, mu_y)
  
  for (node in unique(base_tree$where)) {
    branch <- getBranch(base_tree,rownames(base_tree$frame)[node])
  
    res <- branchInference(base_tree, branch, "reg", alpha=alpha, sigma_y=sigma_y)
    
    CI <- res$confint
    nu <- (base_tree$where == node)/sum(base_tree$where == node)
    
    sample_signal_1 <- t(nu)%*%dat$y
    true_signal_1 <- t(nu)%*%mu_y
    y1 <- dat$y[base_tree$where == node]
    
    write(paste(c(beta,XORlev, seed, "Tree-Values", maxdepth, CI, length(y1), sample_signal_1,true_signal_1,1, ARI), collapse=" "),
          file=filename,append=TRUE)
    
    #### NAIVE METHOD
    CI_naive <- c(mean(y1) - 1.645*sigma_y/sqrt(length(y1)), mean(y1)  + 1.645*sigma_y/sqrt(length(y1)))
    write(paste(c(beta,XORlev, seed, "naive", maxdepth,  CI_naive, length(y1), sample_signal_1,true_signal_1,1, ARI), collapse=" "),
          file=filename,append=TRUE)
  }
  
  ### Now repeat the process but use sample splitting!!!!!!!!!!!!
  for (eps in sampsplit_eps) {
    n1 <- sample(1:n, size=floor(eps*n))
    n2 <- setdiff(1:n, n1)
    dat1 <- dat[n1,]
    dat2 <- dat[-n1,]
    ### CONSIDER: turning off pruning. 
    split_tree <- rpart::rpart(y~., data=dat1, model=TRUE,
                               control=rpart.control(maxcompete=0,maxsurrogate=0, cp=-1,
                                                     xval=0,maxdepth=maxdepth,minbucket=minbucket))
    
    terminalNodes <- sort(unique(split_tree$where))
    split_tree$frame$yval = 1:NROW(split_tree$frame)
    test_predict = predict(split_tree, newdata=dat2)
    all.preds <- predict(split_tree, newdata=dat)
    
    ARIsplit <- mclust::adjustedRandIndex(all.preds , mu_y)
    
    for (node in unique(split_tree$where)) {
      y1 <- dat2[test_predict==node,]$y
      truemu = mean(mu_y[all.preds==node])
      if (NROW(y1)>0) {
        CI.ss <- c(mean(y1) - 1.645*sigma_y/sqrt(length(y1)), mean(y1)  + 1.645*sigma_y/sqrt(length(y1)))
      } else {
        CI.ss <- c(-Inf, Inf)
      }
      write(paste(c(beta,XORlev, seed, "samplesplit", maxdepth,  CI.ss, length(y1), mean(y1), truemu,eps, ARIsplit), collapse=" "),
            file=filename,append=TRUE)
    }
  }
  
  ### DATATHIN
  for (eps in datathin_eps) {
    res <- datathin::datathin(y, family="gaussian", K=2, epsilon=c(eps, 1-eps), arg=sigma_y^2)
    dat1<- dat2 <- dat
    dat1$y <- res[,,1]
    dat2$y <- res[,,2]
  
    split_tree <- rpart::rpart(y~., data=dat1, model=TRUE,
                               control=rpart.control(maxcompete=0,maxsurrogate=0, cp=-1,
                                                     xval=0,maxdepth=maxdepth,minbucket=minbucket))
    
    terminalNodes <- sort(unique(split_tree$where))
    ARIsplit <- mclust::adjustedRandIndex(split_tree$where,mu_y)
    
    for (node in unique(split_tree$where)) {
      ytest <- dat2[split_tree$where==node,]$y
      truemu = mean(mu_y[split_tree$where==node])
      CI_datathin <- c(mean(ytest) - 1.645*sigma_y*sqrt(1-eps)/sqrt(length(ytest)), mean(ytest)  + 1.645*sigma_y*sqrt(1-eps)/sqrt(length(ytest)))
      CI_datathin <- CI_datathin/(1-eps)
      write(paste(c(beta,XORlev, seed, "datathin", maxdepth, CI_datathin, length(ytest), mean(ytest),truemu,eps, ARIsplit), collapse=" "),
            file=filename,append=TRUE)
    }
  }
  
  #### Finally: RRT METHOD!!!
  #require(reticulate)
  rrt_env <- new.env(parent = emptyenv())
  reticulate::py_require(c("cvxpy", "pandas", "numpy", "scikit-learn"))
  reticulate::source_python("RRT_fixed.py", envir = rrt_env)
  RRT_sigmanoise <- c(sigma_y/5, sigma_y/2, sigma_y, sigma_y*2)
  for (noise_sd in RRT_sigmanoise) {
    out <- rrt_env$run_rrt_node_inference(X, y, sdy=sigma_y, seed=as.integer(seed), noise_sd=noise_sd, alpha=alpha, maxdepth=maxdepth)
    ari_rrt <- mclust::adjustedRandIndex(as.factor(out$muhat), mu_y)
    
    for (i in 1:length(out$node_results)) {
      res <- out$node_results[[i]]
      nu <- res$norm_contrast
      if (res$success) {
      write(paste(c(beta,XORlev, seed, "RRT", maxdepth,  res$lower, res$upper, 
                    sum(nu!=0), t(nu)%*%y, t(nu)%*%mu_y,noise_sd, ari_rrt), collapse=" "),
            file=filename,append=TRUE)
      } else {
        write(paste(c(beta,XORlev, seed, "RRT", maxdepth,  -Inf, Inf,
                     NA, NA, NA,noise_sd, ari_rrt), collapse=" "),
              file=filename,append=TRUE)
      }
    }
  }
  
  ##### DATA FISSION?????????
}