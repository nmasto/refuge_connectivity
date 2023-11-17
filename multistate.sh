#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=10
#SBATCH --account=5-33262

#spack load r-tidyverse
#spack load r-rjags
#spack load r-sf
#spack load r-lubridate

spack env activate nick

Rscript ./scripts/Analyses_Mallard_MultistateModel.R