#!/bin/bash -l
#PBS -l walltime=06:00:00,pmem=2580mb,nodes=1:ppn=24
#PBS -m abe
#PBS -M lynch197@umn.edu
#PBS -q mesabi
#PBS -e /home/blekhman/jlynch/lab/lime/lime/profile/mprof_lime_mpi_mesabi.stderr
#PBS -o /home/blekhman/jlynch/lab/lime/lime/profile/mprof_lime_mpi_mesabi.stdout

# do this before module load intel impi
source anaconda/bin/activate mesabimsi

# built anaconda mesabimpi env with module load intel impi
module load intel impi
module list

# change directory so mprof output will be written to the profile directory
cd ~/lab/lime/lime/profile
# -1 on the command line is written by mprof as -- -1 so instead use 20000 which guarantees all SNPs will be processed

mpirun -n 24 \
    mprof run --interval=1.0 --nopython python ~/lab/lime/lime/lime.py \
        ~/lab/glowing-happiness/lime/project/hmp/16S_laf_sadj/data/taxa/finalotu.v35.anterior_nares.txt_laf_sadj_combined_collapsed.txt \
        ~/lab/glowing-happiness/lime/project/hmp/snp/exonic_snps_lasso_mpi.vcf \
        ~/lab/lime/lime/profile/mprof_lime_mpi_output.vcf \
        arcsinsqrt \
        20000 \
        100 \
        no_permutation
