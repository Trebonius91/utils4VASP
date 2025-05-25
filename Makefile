# Makefile for compiling all utils4VASP fortran programs

# Compiler and flags
#FC = ifort  # Intel Fortran
FC = gfortran  # GNU Fortran
FFLAGS = -O2 -m64

# MKL link line (for dynamic linking, LP64 interface)
MKLROOT ?= $(shell echo $$MKLROOT)
MKL_LINK = -L$(MKLROOT)/lib/intel64 \
	-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.so \
	$(MKLROOT)/lib/intel64/libmkl_sequential.so \
	$(MKLROOT)/lib/intel64/libmkl_core.so -Wl,--end-group -lpthread -lm -ldl

# Source and executable lists
SRCS = $(wildcard fortran/*.f90)
EXES = $(SRCS:.f90=)

# Default target: build all executables
all: $(EXES)

# Pattern rule for building each executable
fortran/%: fortran/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(MKL_LINK)

# Clean target
clean:
	find fortran/ -maxdepth 1 -type f ! -name "*.f90" -exec rm -f {} +

.PHONY: all clean
