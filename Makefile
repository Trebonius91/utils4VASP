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
SRCS = $(wildcard setup/*.f90 evaluation/*.f90 ML-FF/*.f90)
EXES = $(SRCS:.f90=)

# Default target: build all executables
all: $(EXES)

# Pattern rule for building each executable
setup/%: setup/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(MKL_LINK)

evaluation/%: evaluation/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(MKL_LINK)

ML-FF/%: ML-FF/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(MKL_LINK)



# Clean target
clean:
	find setup/ -maxdepth 1 -type f ! -name "*.f90" -exec rm -f {} +
	find evaluation/ -maxdepth 1 -type f ! -name "*.f90" -exec rm -f {} +
	find ML-FF/ -maxdepth 1 -type f ! -name "*.f90" -exec rm -f {} +

install:
	mkdir -p $(HOME)/bin
	cp setup/gen_incar.py $(HOME)/bin
	cp setup/gen_poscar.py $(HOME)/bin
	cp setup/analyze_poscar.py $(HOME)/bin
	cp setup/modify_poscar.py $(HOME)/bin
	cp setup/build_adsorb.py $(HOME)/bin
	cp setup/split_freq $(HOME)/bin
	cp evaluation/modify_xdatcar $(HOME)/bin
	cp evaluation/analyze_md $(HOME)/bin
	cp evaluation/analyze_dft $(HOME)/bin
	cp evaluation/check_geoopt.py $(HOME)/bin
	cp evaluation/manage_neb.py $(HOME)/bin
	cp ML-FF/mlff_select $(HOME)/bin
	cp ML-FF/vasp2trainset $(HOME)/bin
	cp ML-FF/mlp_quality $(HOME)/bin

.PHONY: all clean install
