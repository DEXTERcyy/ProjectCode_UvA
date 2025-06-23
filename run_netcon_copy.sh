#!/bin/bash
#SBATCH --job-name=my_r_netcon
#SBATCH --output=netcon_output_copy.txt
#SBATCH --error=analysis_error_copy.txt
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dextercyy@gmail.com

Rscript NetConstruct_Big.r
#Rscript NetConstruct_Big_Days_Filter.r
#Rscript test.r