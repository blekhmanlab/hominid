#!/usr/bin/env bash

# run from lime package directory

mpirun -n 3 \
    python lime.py \
        ../example/lime_example_taxon_table_input.txt \
        ../example/lime_example_snp_input.vcf \
        ../example/lime_example_output.vcf \
        arcsinsqrt \
        -1 \
        100 \
        no_permutation
