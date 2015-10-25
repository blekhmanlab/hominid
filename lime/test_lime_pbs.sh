#!/bin/bash -l
#PBS -l walltime=00:20:00,pmem=2580mb,nodes=1:ppn=3
#PBS -m abe
#PBS -M lynch197@umn.edu
#PBS -q mesabi
#PBS -e /home/blekhman/jlynch/lab/lime/test_lime_pbs.stderr
#PBS -o /home/blekhman/jlynch/lab/lime/test_lime_pbs.stdout

source miniconda/bin/activate venv2_impi
module load impi

mpirun -n 3 \
    python lime.py \
        ~/lab/lime/example/lime_example_taxon_table_input.txt \
        ~/lab/lime/example/lime_example_snp_input.vcf \
        ~/lab/lime/example/lime_example_output.vcf \
        arcsinsqrt \
        -1 \
        100 \
        no_permutation
