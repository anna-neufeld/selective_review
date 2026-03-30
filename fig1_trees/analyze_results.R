library(dplyr)
library(ggplot2)
library(patchwork)
setwd("~/selective_review/")

naivezcol <- "#7CAE00"
sampsplitcol <- "#00BFC4"
treevalcol <- "#F8766D"
datathincol <- "black"

files <- list.files(
  path = ".", 
  pattern = "^Mar26_res_task.*", 
  full.names = TRUE
)

df_list <- lapply(files, function(f) {
  read.table(f, header = FALSE)
})

res <- do.call(rbind, df_list)

names(res) <- c("beta", "XORLev", "seed", "method", "depth", "lower", "upper", "n", "sampsig","truesig", "eps", "rand")

library(dplyr)
library(ggplot2)

res <- res %>% mutate(correctCI = lower < truesig & upper > truesig, length = upper - lower, invalid = n==0)

group_res <- res %>% group_by(method, beta, eps) %>% summarize("cov" = round(mean(correctCI, na.rm=TRUE),3),
                                                               "meanwidth" = mean(length), 
                                                               "medianwidth" = median(length),
                                                               "rand" = mean(rand), "n"=n())

ggplot(data=group_res, aes(x=beta, y=cov, col=method, lty=as.factor(eps)))+geom_line()

ggplot(data=res, aes(x=as.factor(beta),
                     y=rand, col=as.factor(eps)), group=eps)+geom_boxplot()+facet_grid(vars(method))

ggplot(data=res, aes(x=as.factor(beta),
                     y=length, col=as.factor(eps)), group=eps)+geom_boxplot()+facet_grid(vars(method))+
  scale_y_log10()



ggplot(data=res, aes(x=rand, y=length, col=method))+geom_smooth()




