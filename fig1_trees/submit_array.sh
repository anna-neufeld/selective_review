#!/bin/bash
#SBATCH --job-name=trees
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-100
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load R
cd ~/trees_for_cluster
Rscript run_chunk.R 100 Mar30_res
