#!/usr/bin/env bash

python box_bar_plot.py \
  ../example/stability_selection_example_output.txt \
  ../example/lime_example_taxon_table_input.txt arcsinsqrt \
  ../example/box_bar_plots \
  0.5 \
  10