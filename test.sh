#!/bin/bash

#SBATCH --partition=week             
#SBATCH --time=0-01:00:00             
#SBATCH --nodes=1                     
#SBATCH --ntasks-per-node=8           
#SBATCH --mem=128                     
#SBATCH --job-name="Slurm Test"     
#SBATCH --output=test.out          
#SBATCH --error=test.err          
#SBATCH --mail-type=END     

module load mesa/r24.03.1

export MESA_CACHES_DIR=~/.mesa-cache

./mk && ./rn