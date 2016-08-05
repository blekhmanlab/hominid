#!/usr/bin/env bash

export PYTHONPATH=../../lime:$PYTHONPATH

python ../../lime/stability_selection.py \
    ../lime_example_output.rvcf \
    ../data/lime_example_taxon_table_input.txt \
    ../stability_selection_example_output.rvcf \
    0.3 \
    arcsinsqrt \
    -1 \
    --maf-lower-cutoff 0.2
