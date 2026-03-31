setwd("~/trees_for_cluster")
source("one_rep.R")

#### Note: if you need rpart.utils, will need to install from archive because it is not on CRAN. 
library(treevalues)
library(rpart)
library(intervals)

betas <- rep(c(1,2,5), each=10)
n=200
p=20
sigma_y <- 2.5
maxdepth=2
args <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(args) >= 2) args[2] else "my_results"

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
counter <- (task_id-1)*30+1
task_out <- sprintf("%s_task%04d", out_prefix, task_id)

message("Task ", task_id)

for (j in 1:length(betas)) {
  message("  j = ", j, " beta = ", betas[j])
  print(j)
  filename <- task_out
  try(oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=counter, beta=betas[j],
             filename = filename, minbucket=1, maxdepth=2, XORlev=1, alpha=0.1))
  counter <- counter+1
}

message("Task ", task_id, " done. Output prefix: ", task_out)
