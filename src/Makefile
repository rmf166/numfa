# GNU Fortran Compiler
F90    = gfortran
FFLAGS = -O0 -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fcheck=all \
         -Wall -Wextra -fimplicit-none -Wuninitialized -pedantic -Warray-temporaries

#FFLAGS = -O3 -march=native -funroll-loops -fopenmp

## Intel Fortran Compiler
#F90    = ifort
#FFLAGS = -g -C -traceback -mkl

solver.exe : main.f90
	$(F90) $(FFLAGS) -o solver.exe main.f90

clean:
	rm -f solver.exe *.mod *.o

