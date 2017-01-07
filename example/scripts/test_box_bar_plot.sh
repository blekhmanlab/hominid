#!/usr/bin/env bash

export OMP_NUM_THREADS=1  # each process launches only one thread

hominid_box_bar_plot \
  ../stability_selection_example_output.rvcf \
  ../data/hominid_example_taxon_table_input.txt \
  arcsinsqrt \
  ../box_bar_plots \
  0.5 \
  10
