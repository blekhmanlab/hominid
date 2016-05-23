#!/usr/bin/env bash

export PYTHONPATH=../../lime:$PYTHONPATH

python ../../lime/sort_results.py \
  ../stability_selection_example_output.vcf \
  ../data/lime_example_taxon_table_input.txt \
  arcsinsqrt \
  0.05 \
  0.5 \
  10 \
  --extra-columns=GENE,ID