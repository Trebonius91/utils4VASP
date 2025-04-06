#!/usr/bin/env python3
#           
#    gen_incar: Generate the INCAR file for different types of 
#      of typical VASP calculations, that are also described in 
#      the utils4VASP wiki
#    Part of VASP4CLINT
#     Julien Steffen, 2025 (julien.steffen@fau.de)
#        
      
import sys
import os
import re 
import numpy as np
from numpy import linalg as LA


print('''
 This script generates INCAR files for a large variety of different 
 typical VASP calculations that are described in the utils4VASP Wiki.
 The INCAR files are given together with explanations for each keyword,
 allowing the user to change settings from the defaut values if needed.
 The desired calculation is given by a command, of which the following
 are available:
  -single_point : A simple energy+gradient calculation for a given structure
  -geoopt : A geometry optimization to find the next local minimum
  -aimd : An ab-initio molecular dynamics (AIMD) calculation
  -freq : A normal mode frequency calculation (with IR intensities)
  -solvation : A single point calculation with implicit solvation
  -dos : A density of states calculation, resolved to atoms and orbitals
  -band : An electronic band structure calculation for a given path
  -bader : A Bader partial charge calculation
  -cls_is : A core level shift (CLS) calculation (initial state)
  -cls_fs : A core level shift (CLS) calculation (final state)
  -stm : Calculation of the needed data for a simulated STM image
  -neb : A nudged elastic band calculation with the VTST addon
  -ts_opt : A transition state (TS) optimization with the VTST addon
  -steered_aimd : A steered molecular dynamics simulation (AIMD)
  -steered_mlff : A steered molecular dynamics simulation (ML-FF) 
  -mlff_train : A on-the-fly ML-FF training calculation (new or continue)
  -mlff_select : A reselection of ML-FF local reference calculations
  -mlff_refit : Refit a given ML-FF to obtain the fast version of it
  -mlff_md : A ML-FF molecular dynamics calculation, similar to AIMD
''')


header = '''
# Enter the name of your system here: 
SYSTEM = New system

# The PBE exchange correlation functional, mostly used
GGA = PE

# Number of CPUs per orbital, should be divider of available CPUs per node
NCORE = 4

# If the calculation shall be spin-polarized 
ISPIN = 1 # no spin polarization
# ISPIN = 2 # spin polarization

# Maximum number of electronic SCF steps
NELM = 500

# Minimum number of electronic SCF steps
NELMIN = 6

# Electronic cutoff in eV, depends on ENMAX in POTCAR! Set larger for geoopts
ENCUT = 400

# The electronic convergence criterion, energy difference in eV between two cycles
EDIFF = 1E-07

# The used smearing method for electronic states
ISMEAR = 0 # Gaussian smearing, for isolators or unknown systems
# ISMEAR = 1 # Methfessel-Paxton, for metals

# The width of the electronic smearin in eV
SIGMA = 0.04 # standard value for Gaussian smearing
# SIGMA = 0.15 # standard value for Methfessel-Paxton smearing

# The Grimme empirical dispersion correction scheme with Becke-Johnson damping
IVDW = 12

# Write no wave function to file to save memory
LWAVE = .FALSE.

# Evaluation of projection operators, depending on system size
LREAL = Auto  # For "large" systems (more than approx. 30 atoms)
# LREAL = .FALSE. # For "small" systems
'''

print(header)
