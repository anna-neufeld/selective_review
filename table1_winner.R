setwd("~/Dropbox/2025 talks or research/Selective Inference Review Paper/New_Section_2")
library(dplyr)
library(ggplot2)
library(patchwork)

p <- 100
mu.null <- rep(10,p)

mu.alt <- mu.null
delta <- 3.5
mu.alt[1] <- mu.alt[1] + delta

mu.stronger <- mu.null
mu.stronger[1] <- mu.stronger[1] + 2*delta


ntrials <- 2000


alpha <- 0.1
mult <- -qnorm(alpha/2)
mult_bonf <- -qnorm((alpha/p)/2)

res <- data.frame("t" = NA, "method"=NA, "setting"= NA, "sel" = NA, "cov" = NA, "width" = NA)

for (t in 1:ntrials) {
  if (t %%50 == 1) {print(t)}
  set.seed(t)
  
  ### Null setting
  y <- rnorm(p, mean=mu.null, sd=1)
  ihat <- which.max(y)
  
  CI.naive <- y[ihat] + c(-mult, mult)*1
  res <- res %>% add_row(t=t, method="naive", setting="null", sel=ihat, cov = CI.naive[1] < mu.null[ihat] & CI.naive[2] > mu.null[ihat],
                         width=CI.naive[2]-CI.naive[1])
  ### Medium signal
  y <- rnorm(p, mean=mu.alt, sd=1)
  ihat <- which.max(y)
  
  CI.naive <- y[ihat] + c(-mult, mult)*1
  res <- res %>% add_row(t=t, method="naive", setting="strong", sel=ihat, cov = CI.naive[1] < mu.alt[ihat] & CI.naive[2] > mu.alt[ihat],
                         width=CI.naive[2]-CI.naive[1])

  
  ### Strong signal
  y <- rnorm(p, mean=mu.stronger, sd=1)
  ihat <- which.max(y)
  
  ### Naive
  CI.naive <- y[ihat] + c(-mult, mult)*1
  res <- res %>% add_row(t=t, method="naive", setting="stronger", sel=ihat, cov = CI.naive[1] < mu.stronger[ihat] & CI.naive[2] > mu.stronger[ihat],
                         width=CI.naive[2]-CI.naive[1])
}

res <- res[-1,]

res %>% group_by(setting, method) %>% summarize(mean(cov), mean(sel==1), median(width))




