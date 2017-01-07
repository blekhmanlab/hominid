#!/usr/bin/env bash

export OMP_NUM_THREADS=1  # each process launches only one thread

(  hominid_sort_results \
  ../stability_selection_example_output.rvcf \
  ../data/hominid_example_taxon_table_input.txt \
  arcsinsqrt \
  0.05 \
  0.5 \
  10 \
  --extra-columns=GENE,ID \
) > ../sort_results_example_output
