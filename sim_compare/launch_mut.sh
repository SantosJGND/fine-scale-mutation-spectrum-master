#!/bin/bash
#SBATCH -n 1
#SBATCH -t 6:00:00

module purge
module load python/3.6.4

python -u mut_count.py