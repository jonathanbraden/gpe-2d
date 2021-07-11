# Default to GNU compiler without OpenMP, no warnings and no profiling
FC = gfortran
OMP = T
PROF = F
WARN = T
BENCH = F
BENCH_FINE = F
BLAS = F

# Default FFTW Lib location
#FFTW_DIR = /share/apps/fftw-3.3.5/
#FFTW_DIR = /opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/
FFTW_DIR = /usr/local/
FFTW_INC = -I$(FFTW_DIR)include/
FFTW_LIB = -L$(FFTW_DIR)lib/

# Intel Compiler
ifeq ($(FC),ifort)
	FFLAGS= -r8 -fpp
ifeq ($(OMP),T)
	FFLAGS+= -openmp
endif	
	FOPT= -O3 -ipo -xHost #-fast -no-prec-div -parallel
ifeq ($(WARN),T)
	FFLAGS +=
endif
ifeq ($(PROF),T)
	FFLAGS +=
endif

# GNU Compiler
else ifeq ($(FC),gfortran)
	FFLAGS=-fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none 
ifeq ($(OMP),T)
	FFLAGS += -fopenmp
endif
	FOPT=-funroll-loops -O3 #-fblas  (check this blas one to replace matmul calls)
ifeq ($(WARN),T)
	FFLAGS+= -Wall -fbounds-check -pedantic
endif
ifeq ($(PROF),T)
	FFLAGS+= -g -pg
endif
endif

ifeq ($(OMP),T)
	FFTW = -lfftw3_omp -lfftw3
else
	FFTW = -lfftw3
endif

FLIBS = $(FFTW) 
ifeq ($(BLAS),T)
	FLIBS+=-lblas
endif
FLIBS+= -lm

# Define compiler macros
MACROS=
ifeq ($(OMP),T)
	MACROS+= -DUSEOMP=1
endif
ifeq ($(BLAS),T)
	MACROS+= -DUSEBLAS=1
endif

OBJS = constants.o utils.o fftw_mod_cmplx.o lattice.o hamiltonian.o  yoshida.o # output.o

gpe-2d: %: $(OBJS) driver.o
	$(FC) $(FFLAGS) $(FOPT) $(FFTW_INC) $(MACROS) -o gpe-2d driver.o $(OBJS) $(FFTW_LIB) $(FLIBS) $(THREAD_LIB)

%.o: %.f90
	$(FC) $(FFLAGS) $(FOPT) $(FFTW_INC) $(MACROS) -c $< -o $@ $(FFTW_LIB) $(FLIBS) $(THREAD_LIB)

.PHONY : clean

clean:
	rm -f *~
	rm -f *.o
	rm -f *.mod
	rm -f gpe-2d
