# utils4VASP
**Utility scripts and programs for VASP calculations of bulk and interface systems**

by

Julien Steffen, julien.steffen@fau.de

Andreas Mölkner, andreas.moelkner@fau.de

## General overview

This repository contains a list of scripts and programs that can be used to set up, manage and evaluate VASP
calcuations. It was initially designed to support calculations of liquid surface catalysts, 
but the scripts and programs should be general enough to use them for arbitrary VASP calculations
of surface or bulk systems.

A description of all scripts/programs as well as an overview of important VASP calculations and how to do them is given in the [utils4VASP-Wiki](https://github.com/Trebonius91/utils4VASP/wiki)!

The scripts and programs are grouped by the programming language of their implementation and alphabetically within their sections.

Currently included are:

## Python scripts:

 - **built_adsorbates.py** : Place adsorbate atoms or molecules on or in surface slabs. Positions and rotations can be controlled by input files.
 - **build_alloy.py** :  Build unit cells of element alloys on a simple cubic grid (bulk, slab and cluster possible)
 - **check_vac.py (x)** : Checks if single atoms have departed from a surface during ML (utility for **ml_long.sh**!)
 - **eval_neb.py** : Evaluates a nudged elastic band calculation, no matter if already started or finished
 - **integrate_dens.py (x)** : Calculates elemental concentrations in SCALMS from elemental density distributions
 - **manage_mlff_md.py** : Starts and supervises a ML-FF trajectory of a surface slab, restarts if errors occur
 - **modify_poscar.py** : Perform several simple operations on a POSCAR file, such as shifting its atoms, multiply the cel, transform from cartesian to internal coordinates and the other way round

## Bash scripts

 - **manage_MDs.sh** : Management-tool for MD simulations plus example MD.in. 
 - **ml_long.sh** : Do VASP machine learning force field on-the-fly learning trajectories for arbitrary long times, even if the calculation cluster has a walltime limit.
 - **opt_long.sh** : Do a VASP geometry optimization run for arbitrary long times, even if the calculation cluster has a wallime limit

## Fortran programs:

 - **analyze_slab** : Analyze trajectories of surface slab simulations
 - **cut_unitcell** : Cut an arbitrary shaped surface slab unit cell from a given larger surface slab
 - **eval_bader** : Evaluates and visualizes Bader charge calculations
 - **eval_stm** : Generate STM images from VASP partial charge calculations
 - **manage_cls** : Prepares and evaluates core level shift calculations for many atoms in a system
- **modify_xdatcar** : Modify XDATCAR trajectory files: shifts, multiplications, writing of xyz files
 - **mlff_select** : Selects basis functions for a VASP machine learning force field from given ML_AB files
 - **partial_dos** : Extract the partial density of states for certain atoms, elements or orbitals from a DOSCAR file
 - **split_freq** : Splits a large frequency calculation into arbitrary many parts and combine the results after finishing
 - **rdf_pca** : Calculates the weight of chosen components in a time-dependent radial distribution function given from td_rdf
 - **td_rdf** : Calculates time-dependent radial-distribution functions from a XDATCAR trajectory
 - **vasp2trainset**: Generate machine-learned potential training set files from VASP calculations
 - **mlp_quality**: Evaluate the quality of optimized machine-learned potentials against VASP training/verification set
  
