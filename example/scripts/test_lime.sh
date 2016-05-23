#!/usr/bin/env bash

export PYTHONPATH=../../lime:$PYTHONPATH

mpirun -n 3 \
    python ../../lime/lime.py \
        ../data/lime_example_taxon_table_input.txt \
        ../data/lime_example_snp_input.vcf \
        ../lime_example_output.vcf \
        arcsinsqrt \
        -1 \
        100 \
        no_permutation
