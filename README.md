# utils4VASP
**Utility scripts and programs for VASP calculations of bulk and interface systems**

by

Julien Steffen, julien.steffen@fau.de

Andreas Mölkner, andreas.moelkner@fau.de

Maximilian Bechtel, maxi.bechtel@fau.de

## General overview

This repository contains a list of scripts and programs that can be used to set up, manage and evaluate VASP
calcuations. It was initially designed to support calculations of liquid surface catalysts, 
but the scripts and programs should be general enough to use them for arbitrary VASP calculations
of surface or bulk systems.

A description of all scripts/programs as well as an overview of important VASP calculations and how to do them is given in the [utils4VASP-Wiki](https://github.com/Trebonius91/utils4VASP/wiki)!

The scripts and programs are grouped by the area of application (setup, evaluation, ML-FF, management).
Fortran programs have no file ending, ".py" are Python scripts, ".sh" are Bash shell scripts

The Fortran programs are compiled by the Makefile given in the main directory. 
Modify it if needed and then type 'make' to build all Fortran programs.

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

