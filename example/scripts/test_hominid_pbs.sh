#!/bin/bash -l
#PBS -l walltime=00:20:00,pmem=2580mb,nodes=1:ppn=3
#PBS -m abe
#PBS -M lynch197@umn.edu
#PBS -q mesabi
#PBS -e /home/blekhman/jlynch/lab/hominid/test_hominid_pbs.stderr
#PBS -o /home/blekhman/jlynch/lab/hominid/test_hominid_pbs.stdout

source miniconda/bin/activate venv2_impi
module load impi

echo $PYTHONPATH
export PYTHONPATH=~/lab/hominid/hominid:$PYTHONPATH

mpirun -n 3 \
    python ~/lab/hominid/hominid/hominid.py \
        ~/lab/hominid/example/hominid_example_taxon_table_input.txt \
        ~/lab/hominid/example/hominid_example_snp_input.vcf \
        ~/lab/hominid/example/hominid_example_output.vcf \
        arcsinsqrt \
        -1 \
        100 \
        no_permutation \
        --maf-lower-cutoff 0.2
