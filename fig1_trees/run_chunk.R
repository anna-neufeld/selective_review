setwd("~/trees_for_cluster")
source("one_rep.R")

#### NOTE: had to install rpart.utils from ARCHIVE because it is off CRAN
library(treevalues)
library(rpart)
library(intervals)
#library(rpart.utils)

betas <- rep(c(1,2,5), each=10)
n=200
p=20
sigma_y <- 2.5
maxdepth=2
args <- commandArgs(trailingOnly = TRUE)
# optional CLI args: chunk_size and out_prefix
##chunk_size <- if (length(args) >= 1) as.integer(args[1]) else 50
out_prefix <- if (length(args) >= 2) args[2] else "my_results"

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
set.seed(task_id)

#start <- (task_id - 1L) * chunk_size + 1L
#end   <- min(task_id * chunk_size, N)
#if (start > N) {
#  message("Task ", task_id, " has no work (start=", start, " > N=", N, "). Exiting.")
#  quit(status = 0)
#}
# Per-task output file to avoid collisions
task_out <- sprintf("%s_task%04d", out_prefix, task_id)

message("Task ", task_id)

for (j in 1:length(betas)) {
  message("  j = ", j, " beta = ", betas[j])
  print(j)
  filename <- task_out
  try(oneRepFull(n=n,p=p,sigma_y=sigma_y, seed=j, beta=betas[j],
             filename = filename, minbucket=1, maxdepth=2, XORlev=1, alpha=0.1))
}

message("Task ", task_id, " done. Output prefix: ", task_out)
