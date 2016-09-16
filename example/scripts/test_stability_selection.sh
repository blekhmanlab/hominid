#!/usr/bin/env bash

export PYTHONPATH=../../hominid:$PYTHONPATH

python ../../hominid/stability_selection.py \
    ../hominid_example_output.rvcf \
    ../data/hominid_example_taxon_table_input.txt \
    ../stability_selection_example_output.rvcf \
    0.3 \
    arcsinsqrt \
    -1 \
    --maf-lower-cutoff 0.2
