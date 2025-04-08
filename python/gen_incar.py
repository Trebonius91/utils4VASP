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

# If the wave function shall be read from file, deactivated with 0
ISTART = 0

# Take a superposition of atomic charge densities to construct initial density
ICHARG = 2

# The PBE exchange correlation functional, mostly used
GGA = PE

# The precision of the FFT grid
PREC = Normal  # good for most calculations
# PREC = High  # higher precision 

# Number of CPUs per orbital, should be divider of available CPUs per node
NCORE = 4

# If the calculation shall be spin-polarized 
ISPIN = 1 # no spin polarization
# ISPIN = 2 # spin polarization

# The SCF algorithm, either Davidson or RMM-DIIS or a mixture
ALGO = Normal  # Davidson (most stable)
# Algo = Fast  # Davidson + RMM-DIIS
# Algo = VeryFast # only RMM-DIIS 

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
# ISMEAR = -5 # Tetrahedron smearing, for DOS or energy calculations

# The width of the electronic smearin in eV
SIGMA = 0.04 # standard value for Gaussian smearing
# SIGMA = 0.15 # standard value for Methfessel-Paxton smearing

# The Grimme empirical dispersion correction scheme with Becke-Johnson damping
IVDW = 12

# Write no wave function to file to save memory
LWAVE = .FALSE.

# Write no charge density to file to save memory
LCHARG = .FALSE.

# Evaluation of projection operators, depending on system size
LREAL = Auto  # For "large" systems (more than approx. 30 atoms)
# LREAL = .FALSE. # For "small" systems
'''

incar = open("INCAR","w")
for arg in sys.argv[1:]:
   if arg == "-single_point":
      print(" An INCAR for a single point calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a single point energy+gradient calculation")
         print(header)
         print('''
# Set a static calculation, no geoopt or MD
IBRION = -1               

# Set the number of MD steps to zero
NSW = 0
         ''')
      sys.stdout=original_stdout
   elif arg == "-geoopt": 
      print(" An INCAR for a geometry optimization is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a geometry optimization")
         print(header)
         print('''
# Set the algorithm for the optimization, quasi-Newton or conjugate-gradient
# IBRION = 1 # Quasi-Newton, faster, but can explode for bad start structures      
IBRION = 2 # conjugate gradient, slower but more robust                

# How the degrees of freedom in the system shall be treated               
ISIF = 0 # Only atomic positions optimized, cell shape constant
# ISIF = 3  # Cell shape and volume optimized
# ISIF = 8  # Relax only volume, cell shape constant               

# Set the maximum number of optimization steps to a large value
NSW = 1000

# The criterion for geometric convergence
EDIFFG = -0.02  # the smallest gradient component, frequently used value
# EDIFFG = 1E-5 # Energy change between two steps, only use for pathlogical systems

# The size of the geometry optimization step
POTIM = 0.5 

# Number of remembered ionic steps for quasi Newton algorithm
POTIM = 15

# Assumed symmetry of the system
ISYM = -1 # No symmetry, all complex systems with adsorbates etc.
# ISYM = 2 # highly-symmetric systems, such as clean metal surfaces               
         ''')
      sys.stdout=original_stdout
   elif arg == "-aimd":
      print(" An INCAR for ab-initio molecular dynamics (AIMD) is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a ab-initio molecular dynamics (AIMD) calculation")
         print(header)
         print('''
# Set a molecular dynamics calculation
IBRION = 0

# Set the thermostat for NVT (or activate barostat for NpT)               
# MDALGO = 1 # Stochastic Andersen thermostat
MDALGO = 2 # deterministic Nose-Hoover thermostat
# MDALGO = 3 # Langevin thermostat/barostat, NpT is activated!

# Number of MD steps
NSW = 1000
              
# The MD time step in fs. If H is involved no larger than 1 fs!
POTIM = 1.0

# The initial temperature in K
TEBEG = 300

# The final temperature in K (if different to TEBEG, a linear ramp is applied)
TEEND = 300

# The mass of the Nose-Hoover extra degree of freedom
SMASS = 0               

# How new orbitals are predicted from old ones (previous MD step)
IWAVPR = 12

# PAW charge densities up to d-electrons are passed through Broyden mixer  
LMAXMIX = 4

# Number of previous MD steps stored in the Broyden mixer               
MAXMIX = 100

# Only for NpT: Friction coefficients of Langevin thermostat, one value per element!
# LANGEVIN_GAMMA = 5.0 5.0 (...)

# Only for NpT: Friction coefficient for lattics degree of freedom
# LANGEVIN_GAMMA_L = 5.0

# Only for NpT: Ficticious mass of lattice degree of freedom
# PMASS = 5               
         ''')
      sys.stdout=original_stdout
   elif arg == "-freq":
      print(" An INCAR for Hessian matrix/frequency calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a Hessian matrix/frequency calculation")
         print(header)
         print('''
# Activates a numerical frequency calculation
IBRION = 5

# The elongation for numerical calculation in each degree of freedom
POTIM = 0.015

# Activates the dipole moment corrections (for IR intensities)
LDIPOL = .TRUE.

# The direction of the dipole moment corrections
IDIPOL = 3 # For surface slabs (z-axis = surface normal)
# IDIPOL = 4 # For isolated molecules in gas phase/vacuum
               
# The location of the dipole corretion, change if convergence problems occur!
DIPOL = 0.5 0.5 0.5

# No symmetry applied, important for more than one node!
ISYM = -1               
         ''')               
   elif arg == "-solvation":
      print(" An INCAR for an implicit solvation calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is an implicit solvation calculation")
         print(header)
         print('''
# Activates a numerical f
         ''')
      sys.stdout=original_stdout
   elif arg == "-dos":
      print(" An INCAR for a density of states calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a density of states calculation")
         header=header.replace('ISMEAR = 0 ', '# ISMEAR = 0 ')
         header=header.replace('# ISMEAR = -5 ', 'ISMEAR = -5 ')
         header=header.replace('PREC = Normal ', '# PREC = Normal ')
         header=header.replace('# PREC = High ', 'PREC = High ')
         print(header)
         print('''
# Activates possible resolution of DOS into single atoms/orbitals
LORBIT = 11

# Resolution, on how many points the DOS will be calculated
NEDOS = 1001

# The lowest energy (KS-eigenvalue), starting point of DOS plot
# EMIN = -10.0               

# The highest energy (KS-eigenvalue), end point of DOS plot
# EMAX = 10.0
         ''')
   elif arg == "-band":
      print(" An INCAR for a band structure calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a band structure calculation")
         header=header.replace('# Take a superposition of atomic charge densities to construct initial density',
                    '# Set the charge density fixed ') 
         header=header.replace("ICHARG = 2","ICHARG = 11")
         header=header.replace("SIGMA = 0.04 # standard value for Gaussian smearing",
                     "SIGMA = 0.01 # small value for band structure calculations")         
         print(header)
         print('''
         ''')
      sys.stdout=original_stdout
   elif arg == "-bader":
      print(" An INCAR for a Bader charge calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a Bader charge calculation")
         header=header.replace('# Write no charge density to file to save memory',
                     '# Write the charge density to file')
         header=header.replace('LCHARG = .FALSE.','LCHARG = .TRUE.')
         print(header)
         print('''
# Write the all-electron charge density to file
LAECHG = .TRUE.
         ''')
      sys.stdout=original_stdout



   else: 
      print(" Please give one of the options for INCAR!")
      sys.exit(1)
sys.stdout=original_stdout


print(" ")
print(" gen_incar.py exited normally.")
print(" ")
