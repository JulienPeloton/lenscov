## use ifort and f2py at NERSC
ifeq (${NERSC_HOST}, edison)
	FF = ifort
	FPY = f2py
	OPT = --opt=-O3 -lgfortran -lifcore
else ifeq (${NERSC_HOST}, cori)
	FF = ifort
	FF = gfortran
	FPY = f2py
	OPT = --opt=-O3
else
	FF = gfortran
	FPY = f2py-2.7
	OPT = --opt=-ffixed-line-length-none --opt=-O3
endif

all: loop_lensing

loop_lensing:
		${FF} -fPIC -c *.f -lgfortran -lifcore
		${FPY} -c loop_lensing.f90 *.o -m loop_lensing ${OPT}
		${FPY} -c correlation_functions.f90 *.o -m correlation_functions ${OPT}

clean:
		-rm *.o *.so
