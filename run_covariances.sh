#!/bin/bash
## ju@sussex 2015/2016
## Launcher script
## To use the library files, add <path_to_this_file>/src to your PYTHONPATH
## If you just need to recompile fortran codes, use: ./run_covariances.sh recompile
## If you want to use at NERSC (or any cluster), please use the launch_nersc script instead.

main_script=covariances.py

## Input CAMB files
unlensed_file=additional_files/planck_lensing_wp_highL_bestFit_20130627_scalCls.dat
lensed_file=additional_files/planck_lensing_wp_highL_bestFit_20130627_lensedCls.dat

## Runmode: What do you want to compute?
runmode=covariances_phixCMB

## Choice between Planck, Core++, CMB-S4, and Ideal
exp=_test_

if [ "$1" == "recompile" ]
	then
	cd src
	make clean
	make
	exit
elif [ "$1" == "help" ]
	then
	python ${main_script} --help
	exit
elif [ "$1" == "mpi" ]
	then
	time mpirun -np $2 python ${main_script} -input_unlensed_spectra ${unlensed_file} -input_lensed_spectra ${lensed_file} -exp $exp -runmode N1 --mpi
	wait
	time mpirun -np $2 python ${main_script} -input_unlensed_spectra ${unlensed_file} -input_lensed_spectra ${lensed_file} -exp $exp -runmode N0 --mpi
	wait
	time mpirun -np $2 python ${main_script} -input_unlensed_spectra ${unlensed_file} -input_lensed_spectra ${lensed_file} -exp $exp -runmode covariances_CMBxCMB --mpi
	wait
	time mpirun -np $2 python ${main_script} -input_unlensed_spectra ${unlensed_file} -input_lensed_spectra ${lensed_file} -exp $exp -runmode trispB --mpi
	wait
	time mpirun -np $2 python ${main_script} -input_unlensed_spectra ${unlensed_file} -input_lensed_spectra ${lensed_file} -exp $exp -runmode covariances_phixCMB --mpi
	wait
	time mpirun -np $2 python ${main_script} -input_unlensed_spectra ${unlensed_file} -input_lensed_spectra ${lensed_file} -exp $exp -runmode covariances_phixphi --mpi
else
	echo 'In order to run the script, just execute ./run_covariances.sh <option>'
	echo '<option> can be: recompile, help, or mpi #proc'
	exit
fi
