#!/bin/bash
#SBATCH --job-name=my_r_netcon
#SBATCH --output=netcon_output.txt
#SBATCH --error=analysis_error.txt
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G

#Rscript NetConstruct_Big.r
Rscript NetConstruct_Big_Days_Filter.r
#Rscript test.r