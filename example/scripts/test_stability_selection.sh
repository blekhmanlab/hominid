#!/usr/bin/env bash

hominid_stability_selection \
    ../hominid_example_output.rvcf \
    ../data/hominid_example_taxon_table_input.txt \
    ../stability_selection_example_output.rvcf \
    0.3 \
    arcsinsqrt \
    -1 \
    --maf-lower-cutoff 0.2
