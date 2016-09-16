#!/usr/bin/env bash

export PYTHONPATH=../../hominid:$PYTHONPATH

python ../../hominid/box_bar_plot.py \
  ../stability_selection_example_output.rvcf \
  ../data/hominid_example_taxon_table_input.txt \
  arcsinsqrt \
  ../box_bar_plots \
  0.5 \
  10
