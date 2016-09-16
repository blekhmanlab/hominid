#!/usr/bin/env bash

export PYTHONPATH=../../hominid:$PYTHONPATH

python ../../hominid/sort_results.py \
  ../stability_selection_example_output.rvcf \
  ../data/hominid_example_taxon_table_input.txt \
  arcsinsqrt \
  0.05 \
  0.5 \
  10 \
  --extra-columns=GENE,ID
