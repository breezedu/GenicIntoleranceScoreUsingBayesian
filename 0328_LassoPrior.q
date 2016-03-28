#!/bin/sh

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=10
#SBATCH --job-name=0328_Lasso

R CMD BATCH ./0328_LassoPrior.R