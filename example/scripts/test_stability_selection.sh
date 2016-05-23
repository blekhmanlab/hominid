#!/usr/bin/env bash

python stability_selection.py \
    ../example/lime_example_output.vcf \
    ../example/lime_example_taxon_table_input.txt \
    ../example/stability_selection_example_output.vcf \
    0.3 \
    arcsinsqrt \
    -1
