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

# Choose either MKL or, if not available, Lapack and BLAS
LFLAGS = $(MKL_LINK)
#LFLAGS = -llapack -lblas

# Location of the directory for system-wide installation (default: home/bin)
BINDIR = $(HOME)/bin

# Source and executable lists
SRCS = $(wildcard setup/*.f90 evaluation/*.f90 ML-FF/*.f90)
EXES = $(SRCS:.f90=)

# Default target: build all executables
all: $(EXES)

# Pattern rule for building each executable
setup/%: setup/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(LFLAGS)

evaluation/%: evaluation/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(LFLAGS)

ML-FF/%: ML-FF/%.f90
	$(FC) $(FFLAGS) -o $@ $< $(LFLAGS)



# Clean target
clean:
	find setup/ -maxdepth 1 -type f ! -name "*.f90" ! -name "*.py" ! -name "*.sh"  -exec rm -f {} +
	find evaluation/ -maxdepth 1 -type f ! -name "*.f90" ! -name "*.py" ! -name "*.sh" -exec rm -f {} +
	find ML-FF/ -maxdepth 1 -type f ! -name "*.f90" ! -name "*.py" ! -name "*.sh" -exec rm -f {} +

install:
	mkdir -p $(BINDIR)
	cp setup/gen_incar.py $(BINDIR)
	cp setup/gen_poscar.py $(BINDIR)
	cp setup/analyze_poscar.py $(BINDIR)
	cp setup/modify_poscar.py $(BINDIR)
	cp setup/build_adsorb.py $(BINDIR)
	cp setup/split_freq $(BINDIR)
	cp evaluation/modify_xdatcar $(BINDIR)
	cp evaluation/analyze_md $(BINDIR)
	cp evaluation/analyze_dft $(BINDIR)
	cp evaluation/check_geoopt.py $(BINDIR)
	cp evaluation/manage_neb.py $(BINDIR)
	cp ML-FF/mlff_select $(BINDIR)
	cp ML-FF/vasp2trainset $(BINDIR)
	cp ML-FF/mlip_quality $(BINDIR)
	cp management/prepare_slurm.sh $(BINDIR)

.PHONY: all clean install
