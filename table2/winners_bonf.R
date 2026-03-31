library(dplyr)
library(ggplot2)
library(patchwork)

p <- 100
musec3 <- c(4,rep(0,p-1))
ntrials <- 2000


alpha <- 0.1
mult <- -qnorm(alpha/2)
mult_bonf <- -qnorm((alpha/p)/2)

res <- data.frame("t" = NA, "method"=NA, "setting"= NA, "sel" = NA, "cov" = NA, "width" = NA)

for (t in 1:ntrials) {
  if (t %%50 == 1) {print(t)}
  set.seed(t)

  y <- rnorm(p, mean= musec3 , sd=1)
  ihat <- which.max(y)

  
  ### Classical
  CI.naive <- y[ihat] + c(-mult, mult)*1
  res <- res %>% add_row(t=t, method="naive", setting="sec3", sel=ihat, cov = CI.naive[1] <  musec3[ihat] & CI.naive[2] >  musec3[ihat],
                         width=CI.naive[2]-CI.naive[1])
  
  ### Bonferroni
  CI.bonf <- y[ihat] + c(-mult_bonf, mult_bonf)*1
  res <- res %>% add_row(t=t, method="bonf", setting="sec3", sel=ihat, cov = CI.bonf[1] <  musec3[ihat] & CI.bonf[2] >  musec3[ihat],
                         width=CI.bonf[2]-CI.bonf[1])
}

res <- res[-1,]

res %>% group_by(setting, method) %>% summarize(mean(cov), mean(sel==1), median(width))
res %>% group_by(setting, method, sel==1) %>% summarize(mean(cov), mean(sel==1), median(width))
