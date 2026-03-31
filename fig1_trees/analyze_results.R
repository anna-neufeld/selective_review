library(dplyr)
library(ggplot2)
library(patchwork)
setwd("~/selective_review/fig1_trees")

files <- list.files(
  path = "mar30", 
  pattern = "^Mar30_res_task.*", 
  full.names = TRUE
)

df_list <- lapply(files, function(f) {
  read.table(f, header = FALSE)
})

res <- do.call(rbind, df_list)

names(res) <- c("beta", "XORLev", "seed", "method", "depth", "lower", "upper", "n", "sampsig","truesig", "eps", "rand")

res %>% pull(seed) %>% unique %>% sort()

res <- res %>% mutate(correctCI = lower < truesig & upper > truesig, length = upper - lower, invalid = n==0)

res <- res %>% mutate(correctCI = lower < truesig & upper > truesig, length = upper - lower, invalid = n==0)
res %>% pull(seed) %>% unique() 


res <- res %>% 
  filter(beta %in% c(1,2,5)) %>% 
  mutate(signal = "Weak signal") %>% 
  mutate(signal = ifelse(beta==2, "Medium signal", signal)) %>%
  mutate(signal = ifelse(beta==5, "Strong signal", signal)) %>%
  mutate(noise = "All") %>% 
  mutate(length2 = ifelse(is.infinite(length),10000,length))

res %>% group_by(method, signal, noise) %>% summarize(n())

res[res$method %in% c("datathin", "samplesplit") & res$eps==0.5,]$noise = "Least"
res[res$method %in% c("datathin", "samplesplit") & res$eps==0.9,]$noise = "Most"
res[res$method %in% c("datathin", "samplesplit") & res$eps==0.95,]$noise = "Very Low"

res[res$method %in% c("datathin", "samplesplit") & res$eps==0.7,]$noise = "Medium"
res[res$method=="RRT" & res$eps==5,]$noise = "Least"
res[res$method=="RRT" & res$eps==2.5,]$noise = "Medium"
res[res$method=="RRT" & res$eps==1.25,]$noise = "Most"
res[res$method=="RRT" & res$eps==0.5,]$noise = "Very Low"

res$method[res$method=="naive"] <- "Classical"
res$method[res$method=="Tree-Values"] <- "Full CSI"
res$method[res$method=="samplesplit"] <- "Splitting"
res$method[res$method=="datathin"] <- "Thinning"
res$method[res$method=="RRT"] <- "RCSI"


res <- res %>% mutate(signal = ordered(signal, levels=c("Weak signal", "Medium signal", "Strong signal"))) %>%
  mutate(noise = ordered(noise, levels=c("All", "Very Low", "Most", "Medium", "Least")))%>%
  mutate(method = ordered(method, levels=c("Classical", 
                                           "Splitting", "Thinning", "Full CSI", "RCSI")))


res$width_var = ifelse(res$noise=="None", 1,0.8)
group_res <- res %>% group_by(method, signal, noise, width_var) %>% summarize("cov" = round(mean(correctCI, na.rm=TRUE),3),
                                                                              "meanwidth" = mean(length), 
                                                                              "medianwidth" = median(length),
                                                                              "rand" = mean(rand), "n"=n())

### Establish coverage. Boring plot. 
p1 <- ggplot(data=group_res %>% 
               filter(noise != "Very Low"), 
             aes(x=method, y=cov,width=width_var,group=interaction(method,noise), 
                 col=as.factor(noise)))+
  geom_boxplot(position = position_dodge2(preserve = "single"))+theme_bw()+
  facet_grid(rows=vars(signal))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  geom_hline(yintercept=0.9, lty=2)+ylab("Coverage")+ggtitle("Coverage")
p1

### Selection ability, please!
p2 <- ggplot(data=res %>% filter(noise != "Very Low"), 
             aes(x=method, y=rand, group=interaction(method,noise), col=as.factor(noise), width=width_var))+
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  facet_grid(rows=vars(signal))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  ylab("Adjusted Rand Index with True Tree")+ggtitle("Selection")

p3 <- ggplot(data=res %>% filter(noise != "Very Low"), 
             aes(x=method, y=length2, group=interaction(method,noise), 
                 width=width_var, col=as.factor(noise)))+
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  facet_grid(rows=vars(signal))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_log10()+
  ylab("Confidence Interval Length")+ggtitle("Length")

p2 + p1 + p3 +plot_layout(guides="collect") &
  labs(col="Information available \n at selection") &
  xlab("")
ggsave(filename="full_fig.png", width=10, height=5)
