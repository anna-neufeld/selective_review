neg.cond.prof.lik <- function(mymean , X2, X1, eps, sizefac) {
  r.grid <- exp(seq(log(0.001), log(1000), length.out = 100))
  max.lik <- -Inf
  for (r in r.grid) {
    mu_vec = (r+X1)*(sizefac*mymean)*(1-eps)/(r+eps*sizefac*mymean)
    lik <- sum(dnbinom(X2, size=r+X1, mu=mu_vec, log=T))
    if (lik > max.lik) {max.lik <- lik}
  }
  -max.lik
}

neg.cond.prof.lik.alt <- function(mymean , X2, X1, eps, sizefac, clusters) {
  r.grid <- exp(seq(log(0.001), log(1000), length.out = 100))
  max.lik <- -Inf
  for (r in r.grid) {
    mu_vec1 = (r+X1[clusters==1])*(sizefac[clusters==1]*mymean[1])*(1-eps)/(r+eps*sizefac[clusters==1]*mymean[1])
    lik1 <- sum(dnbinom(X2[clusters==1], size=r+X1[clusters==1], mu=mu_vec1, log=T))
    mu_vec2 = (r+X1[clusters==2])*(sizefac[clusters==2]*mymean[2])*(1-eps)/(r+eps*sizefac[clusters==2]*mymean[2])
    lik2 <- sum(dnbinom(X2[clusters==2], size=r+X1[clusters==2], mu=mu_vec2, log=T))
    if (lik1 + lik2 > max.lik) {
      max.lik <- lik1+lik2
    }
  }
  -max.lik
}

fission_test_sizefac <- function(X1, X2, clusters, eps, sizefac) {
  
  mod <- MASS::glm.nb(X1~offset(log(sizefac)))
  r.hat.n <- mod$theta
  mean.init <- as.numeric(exp(mod$coefficients[1])/eps)
  result.null <- optim(
    par =  mean.init,               
    fn = neg.cond.prof.lik,                
    X2 = X2, X1 = X1, eps = eps, 
    sizefac=sizefac,
    method = "L-BFGS-B",               
    lower = 0.01, upper = 100              
  )
  
  result.alt1 <- optim(
    par =  result.null$par,               
    fn = neg.cond.prof.lik,                
    X2 = X2[clusters==1], X1 = X1[clusters==1], eps = eps, 
    sizefac=sizefac[clusters==1],
    method = "L-BFGS-B",               
    lower = 0.01, upper =  100              
  )
  result.alt2 <- optim(
    par =  result.null$par,               
    fn = neg.cond.prof.lik,    
    sizefac=sizefac[clusters==2],
    X2 = X2[clusters==2], X1 = X1[clusters==2], eps = eps, 
    method = "L-BFGS-B",               
    lower = 0.01, upper = 100              
  )
  
  max.log.lik.null <- -result.null$value
  max.log.lik.alt <- -(result.alt1$value+result.alt2$value)
  LRT <- -2*(max.log.lik.null- max.log.lik.alt)
  1-pchisq(LRT, df=2)
}

fission_test_sizefac <- function(X1, X2, clusters, eps, sizefac) {
  
  mod <- MASS::glm.nb(X1~offset(log(sizefac)))
  r.hat.n <- mod$theta
  mean.init <- as.numeric(exp(mod$coefficients[1])/eps)
  result.null <- optim(
    par =  mean.init,               
    fn = neg.cond.prof.lik,                
    X2 = X2, X1 = X1, eps = eps, 
    sizefac=sizefac,
    method = "L-BFGS-B",               
    lower = 0.0001, upper = 100              
  )
  
  result.alt1 <- optim(
    par =  result.null$par,               
    fn = neg.cond.prof.lik,                
    X2 = X2[clusters==1], X1 = X1[clusters==1], eps = eps, 
    sizefac=sizefac[clusters==1],
    method = "L-BFGS-B",               
    lower = 0.0001, upper =  100              
  )
  result.alt2 <- optim(
    par =  result.null$par,               
    fn = neg.cond.prof.lik,    
    sizefac=sizefac[clusters==2],
    X2 = X2[clusters==2], X1 = X1[clusters==2], eps = eps, 
    method = "L-BFGS-B",               
    lower = 0.0001, upper = 100              
  )
  
  max.log.lik.null <- -result.null$value
  max.log.lik.alt <- -(result.alt1$value+result.alt2$value)
  LRT <- -2*(max.log.lik.null- max.log.lik.alt)
  1-pchisq(LRT, df=2)
}


fission_test_sizefac_new <- function(X1, X2, clusters, eps, sizefac) {
  mod <- MASS::glm.nb(X1~offset(log(sizefac)))
  r.hat.n <- mod$theta
  mean.init <- as.numeric(exp(mod$coefficients[1])/eps)
  
  result.null <- optim(
    par =  mean.init,               
    fn = neg.cond.prof.lik,                
    X2 = X2, X1 = X1, eps = eps, 
    sizefac=sizefac,
    method = "L-BFGS-B",               
    lower = 0.0001, upper = 100              
  )
  
  mean.init.2 <- c(result.null$par,result.null$par)

  result.alt <- optim(
    par = mean.init.2,
    fn = neg.cond.prof.lik.alt,
    X2=X2, X1=X1, eps=eps, sizefac=sizefac, clusters=clusters,
    method = "L-BFGS-B",               
    lower = 0.0001, upper = 100    
  )
  
  max.log.lik.null <- -result.null$value
  max.log.lik.alt <-  -result.alt$value
  LRT <- -2*(max.log.lik.null- max.log.lik.alt)
  1-pchisq(LRT, df=1)
}



