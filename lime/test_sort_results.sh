#!/usr/bin/env bash

python sort_results.py \
  ../example/stability_selection_example_output.vcf \
  ../example/lime_example_taxon_table_input.txt \
  arcsinsqrt \
  0.05 \
  0.5 \
  10 \
  --extra-columns=GENE,ID