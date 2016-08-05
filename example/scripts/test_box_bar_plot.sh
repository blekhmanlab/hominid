#!/usr/bin/env bash

export PYTHONPATH=../../lime:$PYTHONPATH

python ../../lime/box_bar_plot.py \
  ../stability_selection_example_output.rvcf \
  ../data/lime_example_taxon_table_input.txt \
  arcsinsqrt \
  ../box_bar_plots \
  0.5 \
  10
