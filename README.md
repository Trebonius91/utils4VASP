# utils4VASP
**Utility scripts and programs for VASP calculations of bulk and interface systems**

## License

utils4VASP is distributed unter the terms of the [MIT license](https://opensource.org/licenses/mit-license):

```
   utils4VASP - Utility scripts and programs for VASP calculations
                of bulk and interface systems  

   Copyright (c) 2025 by Julien Steffen (julien.steffen@fau.de)
                         Andreas Mölkner (andreas.moelkner@fau.de)
                         Maximilian A. Bechtel (maxi.bechtel@fau.de)

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
```

## General overview

This repository contains a list of Python3 scripts and Fortran programs that can be used to set up, 
manage and evaluate VASP calcuations, both periodic DFT and ML-FF. 

A description of all scripts/programs as well as an overview of important VASP calculations and how to do them is given in the [utils4VASP-Wiki](https://github.com/Trebonius91/utils4VASP/wiki)!

If you use utils4VASP for your research, we would be happy if you could cite the [utils4VASP paper](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c01491).

The scripts and programs are grouped by the area of application (setup, evaluation, ML-FF, management) and stored into subdirectories of those names in the utils4VASP directory.
Fortran programs have no file ending, ".py" are Python scripts, ".sh" are Bash shell scripts

The Fortran programs are compiled by the Makefile given in the main directory. 
By default, the compilation is done with gfortran, utilizing the intel MKL libraries for Lapack and BLAS routines.
Change the compiler (FC = ..) and linker flags (LFLAGS = ...) if needed for your system.

For the Python routines, only the numpy and scipy packages are needed besides the standard os, sys, re, random etc packagesthat are part of the Python3 installation itself.

Compilation is done as usual with the following command in the utils4VASP main folder:

    make

Further, all scripts and programs can be copied to a suitable directory to make them accessible systemwide. By default, this is the $HOME/bin directory, which can be changed by altering the (BINDIR = ...) flag in the Makefile. Installation then is done by:

    make install

Compiled fortran programs can be removed with:

    make clean 


Currently included are:

## Setup:
 - **gen_incar.py** : Generate INCAR file templates for several job types
 - **gen_poscar.py** : Generate POSCAR of alloy and surface structures
 - **analyze_poscar.py** : Analyze POSCAR, generate KPOINTS and POTCAR
 - **modify_poscar.py** : Modify POSCAR: multiply, transform, shift, insert etc.
 - **build_adsorb.py** : Place adsorbates on surfaces, set translation, rotation
 - **split_freq** : Divide frequency calculations for large molecules

## Evaluation:
 - **modify_xdatcar** : Modify trajectory files: multiply, shift, pick etc.
 - **analyze_md** : Analyze MD trajectories for RDFs, diffusion, density etc.
 - **analyze_dft** : Analyze DFT caculations (Bader cahrges, STM, CLS, pDOS)
 - **check_geoopt.py** : Monitor geometry optimizations with selective dynamics
 - **manage_neb.py** : Setup, monitor and restart NEB calculations

## ML-FF
 - **mlff_select** : Heuristic selection of local reference configurations
 - **eval_vasp_ml.py** : Visualize results of VASP ML-FF on the fly learnings
 - **vasp2trainset** : Generate ML-FF training sets from VASP calculations
 - **mlp_quality** : Determine quality of MLPs for VASP validation set
 - **check_conv.py** : Check SCF convergence of vasp2trainset single points

## Management
 - **md_long.sh** : Automated restart of MD calculations with slurm
 - **opt_long.sh** : Automated restart of geometry optimizations with slurm
 - **ml_long.sh** : Automated restart of VASP ML-FF on the fly learnings
 - **mace_long.sh** : Automated restart of MACE MD trajectories with ASE
 - **manage_MDs.sh** : Management-tool for MD simulations plus example MD.in.

