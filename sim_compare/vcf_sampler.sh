#!/bin/bash
#SBATCH -n 1
#SBATCH -t 6:00:00

module purge
module load vcftools/0.1.12a
module load tabix/0.2.6
module load python/3.6.4

vcf_file="/home/jgarc235/Human/data_vcf/phase1_chr1_filtered.vcf.gz.recode.vcf.gz"
N=2
L=1000000

assembly="hg19"
batch="Harris"
out_dir="/home/jgarc235/Human/sim_compare/data/"
diffs="/home/jgarc235/Human/mutation_study/data/hg19_panTro4_diffs/hg19_panTro4_diffs_chr"
ids="/home/jgarc235/Human/ind_assignment/ind_assignments_onePop.txt"

chrom="1"

python -u ordered_extractions.py -v $vcf_file -c $chrom -n $N -l $L -a $assembly \
-b $batch -d $diffs -i $ids -o $out_dir

