#!/usr/bin/env python
# Copyright (C) 2016 Peloton
"""Distutils based setup script for lenscov.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. To compile Fortran code for internal
computation inplace, just type the command (see README also):

    python setup.py build_ext --inplace --fcompiler=gfortran

For most purposes, you do not need to do anything else.
For the easiest installation just type the command (you'll probably need
root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py bench -> will run the complete benchmark suite
    python setup.py audit -> will run pyflakes checker on source code

To get a full list of avaiable commands, read the output of:

    python setup.py --help-commands
"""
import os
import numpy as np
import distutils

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('src',parent_package,top_path)

    if distutils.version.StrictVersion(np.version.version) > distutils.version.StrictVersion('1.6.1'):
        config.add_extension('loop_lensing', sources = ['src/loop_lensing.f90'],
                             libraries=[], f2py_options=[],
                             extra_f90_compile_args=['-ffixed-line-length-1000', '-O3'],
                             extra_compile_args=[''], extra_link_args=['src/*.o'],)
        config.add_extension('correlation_functions', sources = ['src/correlation_functions.f90'],
                             libraries=[], f2py_options=[],
                             extra_f90_compile_args=['-ffixed-line-length-1000', '-O3'],
                             extra_compile_args=[''], extra_link_args=['src/*.o'],)
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup

    ## Wigner 3-j symbols
    os.system("cd src/; gfortran -fPIC -c *.f -lgfortran -lifcore; cd ..")

    setup(name='lenscov',
        configuration=configuration,
        version='1.0.0',
        author='Julien Peloton',
        author_email='j.peloton@sussex.ac.uk',
        packages=['lenscov'],
        description='CMB and lensing covariances',)
