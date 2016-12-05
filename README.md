lenscov
==

#### The package
This package contains scripts for computing CMB and lensing covariances
as described in Peloton et al. 2016 (<a href="https://arxiv.org/abs/1611.01446" title=1611.01446">1611.01446</a>)
It has routines for:
* Spectra
    * computing CMB lensed power spectra using series-expansion
    * computing CMB lensed power spectra using correlation functions
    * computing (analytically) the N0 bias
    * computing (analytically) the N1 bias
* Covariances
    * computing (analytically) the covariance matrix of
      CMB lensed spectra up to second order in lensing potential power-spectrum
    * computing (analytically) the covariance matrix of
      lensing potential spectrum
    * computing (analytically) the cross-covariance matrix of
      CMB lensed spectra and lensing potential spectrum up to
      second order in lensing potential power-spectrum
* Additional
    * computing derivatives of spectra using correlation functions
    * computing wigner (small) d-matrix efficiently

The code performs computation in the full-sky formalism, and
the code is fully-parallelized using mpi4py (but some parts work also in serial).
In addition, some Fortran parts make use of openmp (computation of N1).

### Before starting
This code has the following dependencies:
* numpy, pylab, etc (required)
* scipy (required, for C interfacing)
* f2py (required, for Fortran interfacing)
* mpi4py (optional, to use parallel computing)

Make sure you update your PYTHONPATH to use the code.
Just add in your bashrc:
```bash
COVPATH=/path/to/the/package
export PYTHONPATH=$PYTHONPATH:$COVPATH:$COVPATH/src
```

### Compilation of C and Fortran
The code is primarily written in Python, although the low-level routines
(wigner3-j, nested loops, etc) are implemented in Fortran and C,
and interfaced with Python. The compilation of C is done on-the-fly
through scipy.weave. The compilation of Fortran is done either using the
setup file, by running:
```bash
python setup.py build_ext --inplace --fcompiler=gfortran
```
Depending on system, you may need to specify a different fortran compiler.
The compilation has been tested successfully with gfortran (gcc version 5.1.0) and ifort.
For the purist, we also provide a Makefile for a direct compilation.

### Usage
You can find a launcher (run_covariances.sh) to have a quick usage.
We also provide a batch file for usage on cluster.
You have pre-defined experimental configurations in misc.py.

### License
GNU License (see the LICENSE file for details) covers all files
in the lenscov repository unless stated otherwise.
