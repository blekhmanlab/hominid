# HOMINID
Python MPI program to identify associations between host genetic variation
and microbiome taxonomic composition.

**Summary**: Recent studies have uncovered a strong effect of host genetic variation on the  composition of host-associated microbiota. Here, we present HOMINID, a computational approach based on Lasso linear regression, that given host genetic variation and microbiome composition data, identifies host SNPs that are correlated with microbial taxa abundances. By using HOMINID on data from the Human Microbiome Project, we identified 13 human SNPs in which genetic variation is correlated with microbiome taxonomic composition in 15 body sites.

We also present a tool for visualization of host-microbiome association network identified in HOMINID, currently including toy data representing all SNP-microbe associations with a nominal p-value <= 0.1.  Online visualization tool at http://z.umn.edu/genemicrobe

Contact: blekhmanlab@gmail.com

# Running the HOMINID software

This README describes the installation process and how to test the
`HOMINID` software on included test data. Once `HOMINID` is installed
and known to work read the following documents for instructions on using
it with your own data:

1. [HOMINID analysis pipeline](https://github.com/blekhmanlab/hominid/wiki/HOMINID-pipeline)
2. [Running `hominid` on your data](https://github.com/blekhmanlab/hominid/wiki/Running-hominid-on-your-data)
3. [Running `hominid_stability_selection` on your data](https://github.com/blekhmanlab/hominid/wiki/Running-stability-selection-on-your-data)
4. [Running `hominid_sort_results` on your data](https://github.com/blekhmanlab/hominid/wiki/Running-sort-results-on-your-data)

## Requirements
`HOMINID` is a Python 2.7/3.5+ MPI program. It is intended to run on a cluster,
but it will run anywhere with a working MPI implementation and
[mpi4py](http://mpi4py.readthedocs.io/en/stable/). `HOMINID` has been
tested only on Linux operating systems.

The required Python packages will be automatically installed. They are:

 - mpi4py (version 2.0 or greater)
 - numpy
 - pandas
 - SciPy
 - scikit-learn (version 0.19.1)
 - scikits.bootstrap

The optional plotting script requires R 3.2+ and rpy2. These can be
installed using the Anaconda Python distribution with the `r-essentials`
packages as shown below.

## Computer requirements
`HOMINID` is a multiprocess program that benefits from multiple cores and multiple
processors.  It can run on any hardware that supports mpi4py from laptops to clusters.

## Install
It is recommended that hominid be installed in a Python virtual environment.
These instructions are specifically for the [Miniconda3](https://conda.io/miniconda.html)
distribution, which has been tested with `HOMINID`.

### Install a MPI Implementation
A MPI implementation is available on most clusters so this step is generally
necessary only on laptop and desktop computers. `HOMINID` is known to run with OpenMPI on Ubuntu 14.04 and Ubuntu 16.04. These commands will install OpenMPI on Ubuntu and similar Debian-based Linux distributions:

```
$ sudo apt update
$ sudo apt install mpi-default-dev openmpi-bin
```

### Create a Python Virtual Environment
1. Download and install Miniconda3 from a terminal with these commands:
```
$ wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
$ chmod u+x miniconda3.sh
$ ./miniconda3.sh
```
  + You will be asked to accept the license as part of the installation process.
  + Accept the default installation directory unless you have good reason not to. The default directory is ~/miniconda3.
  + The installer will ask to modify your PATH variable to include `conda`. Choose `yes`. If you choose `no` you will need to use full pathnames in subsequent steps.

2. Close the terminal and open a new one so the update to your PATH variable takes effect. In the new terminal update `conda`:
```
$ conda update conda
```

3. Create a new virtual environment, which will be installed in the ~/miniconda3/envs directory. Here the virtual environment is named `hom` but another name will work.
```
$ conda create -n hom pip
```
The `r-essentials` packages support the box plot script and are not
strictly required.

4. Install `HOMINID` in the virtual environment with the following commands:
```
$ git clone https://github.com/blekhmanlab/hominid.git
$ cd hominid
$ source activate hom
(hom) $ pip install -r requirements.txt
(hom) $ conda install r-essentials
```

Once `HOMINID` has been installed with `pip` the scripts can be executed from any directory by name as follows:

```
(hom)$ hominid
(hom)$ hominid_stability_selection
(hom)$ hominid_sort_results
(hom)$ hominid_box_bar_plot
```
The Python programs are in directory `hominid/hominid`.
Test scripts are in directory `hominid/example/scripts`.
Test input data files are in directory `hominid/example/data`.

## Test the Installation
The installation can be tested using the included test scripts and
test data with the following steps.

1. Run `HOMINID` on the sample data. Change directory to the `example/scripts` directory and run `test_hominid.sh`.
```
(venv_hominid)$ cd example/scripts
(venv_hominid)$ ./test_hominid.sh
```
In `test_hominid.sh`, the option `-n 3` to `mpirun` specifies that 3 processes will be used. Change this if you want to use a different number of processes. Performance will be reduced if more processes are specified than available cores. A minimum of two processes must be specified.
The test output is written to `hominid/example/hominid_example_output.rvcf`.
Many lines will be printed to `stderr` so you can watch hominid's progress.

2. Run hominid on the sample data, with permuted sample IDs:
```
(venv_hominid)$ ./test_hominid_permute.sh
```
Output is written to `hominid/example/hominid_example_output_permute.rvcf`.

3. Run stability selection to find associated OTUs/taxa/covariates.
```
(venv_hominid)$ ./test_stability_selection.sh
```
Output is written to `hominid/example/stability_selection_example_output.rvcf`.

4. Combine the Lasso regression results with the microbiome abundances:
```
(venv_hominid)$ ./test_sort_results.sh
```
Output is written to `hominid/example/sort_results_example_output`.
