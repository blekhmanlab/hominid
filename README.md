# lime
Python MPI program using LASSO regression to find associations between host genetics and microbiome.

# Requirements
lime is a Python 2.7 MPI program intended to run on a cluster but it will run anywhere with a working MPI implementation and mpi4py installed.

 - mpi4py
 - numpy
 - pandas
 - SciPy
 - scikit-learn
 - scikits.bootstrap

# Use
See the [wiki pages](https://github.com/jklynch/lime/wiki) for 
instructions on how to run the software, and the 
[example directory](https://github.com/jklynch/lime/tree/master/example) for 
sample bash scripts and input data files.

The software for the Lasso analysis pipeline is in 
the [lime directory](https://github.com/jklynch/lime/tree/master/lime).
