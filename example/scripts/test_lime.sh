#!/usr/bin/env bash

export PYTHONPATH=../../lime:$PYTHONPATH

mpirun -n 3 \
    python ../../lime/lime.py \
        ../data/lime_example_taxon_table_input.txt \
        ../data/lime_example_snp_input.rvcf \
        ../lime_example_output.rvcf \
        arcsinsqrt \
        -1 \
        100 \
        no_permutation \
        --maf-lower-cutoff 0.2
