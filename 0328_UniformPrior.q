#!/bin/sh

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=2
#SBATCH --job-name=0328_Uniform

R CMD BATCH ./0328_UniformPrior.R