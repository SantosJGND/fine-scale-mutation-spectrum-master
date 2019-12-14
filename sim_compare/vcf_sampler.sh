#!/bin/bash
#SBATCH -n 1
#SBATCH -t 6:00:00

module purge
module load vcftools/0.1.12a
module load tabix/0.2.6
module load python/3.6.4

vcf_file="/home/jgarc235/Rhesus/Picard/rheMac10_data/rheMac8_lift/rheMac8_full_lift_rheMac10.vcf.gz"
N=2
L=1000000

assembly="rheMac10"
batch="test"
out_dir="/home/jgarc235/Sim_compare/data/"

python -u ordered_extractions.py -v $vcf_file -n $N -l $L -a $assembly \
-b $batch -o $out_dir
