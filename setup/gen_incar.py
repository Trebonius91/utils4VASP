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


#
#    only print the overview of utils4VASP scripts/programs and stop
#
for arg in sys.argv:
   if arg == "-overview":
      print('''
utils4VASP: Setup and Evaluation of DFT and MLP simulations with VASP
The following scripts and programs are contained:
Setup:
 - gen_incar.py: Generate INCAR file templates for several job types
 - gen_poscar.py: Generate POSCAR of alloy and surface structures
 - analyze_poscar.py: Analyze POSCAR, generate KPOINTS and POTCAR
 - modify_poscar.py: Modify POSCAR: multiply, transform, shift, insert etc.
 - build_adsorb.py: Place adsorbates on surfaces, set translation, rotation
 - split_freq: Divide frequency calculations for large molecules
Evaluation:
 - modify_xdatcar: Modify trajectory files: multiply, shift, pick etc.
 - analyze_md: Analyze MD trajectories for RDFs, diffusion, density etc.
 - analyze_dft: Analyze DFT calculations (Bader charges, STM, CLS, pDOS)
 - check_geoopt.py: Monitor geometry optimizations with selective dynamics
 - manage_neb.py: Setup, monitor and restart NEB calculations
ML-FF:
 - mlff_select: Heuristic selection of local reference configurations
 - vasp2trainset: Generate ML-FF training sets from VASP calculations
 - mlp_quality: Determine quality of MLPs for VASP validation set
Management:
 - md_long.sh: Automated restart of MD calculations with slurm
 - opt_long.sh: Automated restart of geometry optimizations with slurm
 - ml_long.sh: Automated restart of VASP ML-FF on the fly learnings
 - mace_long.sh: Automated restart of MACE MD trajectories with ASE
''')
      sys.exit(0)

#
#    Print general information and all possible keywords of the program    
#
print('''
  *** utils4VASP -- utility scripts and programs for VASP ***
 gen_incar.py: Generates INCAR files for a large variety of different 
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
  -steered : A steered molecular dynamics simulation (AIMD)
  -mlff_train : A on-the-fly ML-FF training calculation (new or continue)
  -mlff_select : A reselection of ML-FF local reference calculations
  -mlff_refit : Refit a given ML-FF to obtain the fast version of it
  -mlff_md : A ML-FF molecular dynamics calculation, similar to AIMD
 Give the -overview command for an overview of utils4VASP.
''')

#
#    Print INCAR header for periodic DFT calculations    
#
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
# ISMEAR = -1 # Fermi, for AIMD
# ISMEAR = -5 # Tetrahedron smearing, for DOS or energy calculations

# The width of the electronic smearin in eV
SIGMA = 0.04 # standard value for Gaussian smearing
# SIGMA = [8.61733E-5 * Temp]  # scales with temperature for Fermi smearing
# SIGMA = 0.15 # standard value for Methfessel-Paxton smearing

# The Grimme empirical dispersion correction scheme with Becke-Johnson damping
IVDW = 12

# Decides which degrees of freedom are flexible 
ISIF = 2 # only atoms are movable, stress tensor is calculated
# ISIF = 3 # Atoms and cell shape fully movable 
# ISIF = 7 # Only volume of cell flesible, constant shape
# ISIF = 8 # Volume and atoms movable, constant shape

# Decides if symmetry is used to speedup the calculations
ISYM = 2 # use symmetry (single points and symmetric geoopts)
# ISYM = 0 # Turn off most symmetry (for MD)
# ISYM = -1 # Turn off symmetry completely (for freq and complex geoopts)

# Write no wave function to file to save memory
LWAVE = .FALSE.

# Write no charge density to file to save memory
LCHARG = .FALSE.

# Evaluation of projection operators, depending on system size
LREAL = Auto  # For "large" systems (more than approx. 30 atoms)
# LREAL = .FALSE. # For "small" systems
'''

incar = open("INCAR","w")
if (len(sys.argv) < 2):
   print(" Please give one of the possible calculations!")
   print(" ")
   sys.exit(1)

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
         header=header.replace('ISYM = 2 ', '# ISYM = 2 ')
         header=header.replace('# ISYM = -1 ', 'ISYM = -1 ')
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
         ''')
      sys.stdout=original_stdout
   elif arg == "-aimd":
      print(" An INCAR for ab-initio molecular dynamics (AIMD) is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is an ab-initio molecular dynamics (AIMD) calculation")
         header=header.replace('ISYM = 2 ', '# ISYM = 2 ')
         header=header.replace('# ISYM = 0 ', 'ISYM = 0 ')
         header=header.replace('ISMEAR = 0 ', '# ISMEAR = 0 ')
         header=header.replace('# ISMEAR = -1 ', 'ISMEAR = -1 ')
         header=header.replace('SIGMA = 0.04 ', '# SIGMA = 0.04 ')
         header=header.replace('# SIGMA = [8.61733E-5 ', 'SIGMA = [8.61733E-5 ')
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
         header=header.replace('ISYM = 2 ', '# ISYM = 2 ')
         header=header.replace('# ISYM = -1 ', 'ISYM = -1 ')
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
      sys.stdout=original_stdout          
   elif arg == "-solvation":
      print(" An INCAR for an implicit solvation calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is an implicit solvation calculation")
         header=header.replace('# Take a superposition of atomic charge densities to construct initial density', 
                           '# Read the initial density/wavefunction from WAVECAR file ')
         header=header.replace('ISTART = 0', 'ISTART = 1')
         print(header)
         print('''
# Activates the solvation calculation
LSOL = .TRUE.

# The relative permittivity of the solvent (e.g., for water: 78.4)
EB_k = [value]               
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
      sys.stdout=original_stdout
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
         header=header.replace('PREC = Normal ', '# PREC = Normal ')
         header=header.replace('# PREC = High ', 'PREC = High ')

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
   elif arg == "-cls_is":
      print(" An INCAR for a core level shift calculation (initial state) is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a core level shift (initial state) calculation")
         header=header.replace('PREC = Normal ', '# PREC = Normal ')
         header=header.replace('# PREC = High ', 'PREC = High ')
         print(header)
         print('''
# Calculate the energies for the core level state electrons
ICORELEVEL = 1
         ''')
      sys.stdout=original_stdout
   elif arg == "-cls_fs":
      print(" An INCAR for a core level shift calculation (final state) is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a core level shift (final state) calculation")
         header=header.replace('PREC = Normal ', '# PREC = Normal ')
         header=header.replace('# PREC = High ', 'PREC = High ')
         print(header)
         print('''
# Activates the final state mode for the core level calculation
ICORELEVEL = 2

# Number of species for which the core level shall be calculated (number in POTCAR)
CLNT = [number]

# Main quantum number of the excited electron (1: first shell, 2: second shell, ...)
CLN = [number]

# Angular quantum number L of the excited electron (0:s, 1:p, 2:d, ...)
CLL = [number]

# The number of electrons to be excited from the orbital, 0.5 is good
CLZ = 0.5               
         ''')
      sys.stdout=original_stdout
   elif arg == "-stm":
      print(" An INCAR for a STM image calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a STM image calculation")
         print(header)
         print('''
# Evaluate partial charge density
LPARD = .TRUE.

# Read in the existent WAVECAR file instead of SCF calcuation
ISTART = 1

# The energy range for partial charge density 
# (first value experimental STM current (*-1), second zero)
EINT = -1 0

# Calculate partial charge density in an energy range given by EINT
NBMOD = -3
         ''')
      sys.stdout=original_stdout
   elif arg == "-neb":
      print(" An INCAR for a nudged elastic band (NEB) calculation is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a nudged elastic band (NEB) calculation")
         print(header)
         print('''
# Number of NEB images (usually 8, 16 or 32)
IMAGES = [number]

# The spring constant between two images along the path, -5 is default
SPRING = -5

# Activates the NEB method
ICHAIN = 0

# Activates the climbing image method, TS will be optimized as well
LCLIMB = .TRUE.

# Turns on the SS-NEB routine
LNEBCELL = .FALSE.

# Use a special optimization routine from the Henkelman group (Quick min)
IOPT = 3

# Disable VASP default optimizers for IOPT = 3
IBRION = 3

# Sets the VASP geometry optimization step to zero (all will be done by VTST)
POTIM = 0

# Dynamical NEB time step for optimization
TIMESTEP = 0.1               
         ''')
      sys.stdout=original_stdout
   elif arg == "-ts_opt":
      print(" An INCAR for a transition state (TS) optimization is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a transition state (TS) optimization")
         print(header)
         print('''
# Activates the dimer method
ICHAIN = 2

# MD with zero time step
IBRION = 3

# Zero time step for MD optimization
POTIM = 0

# The dimer separation, twice distance between the images
DdR = 0.005

# Number of rotation steps per dimer translation
DRotMax = 4

# Lower threshold for dimer rotation (rotational force)
DFNMin = 0.01

# Upper threshold for dimer rotation
DFNMax = 1.0               
         ''')
      sys.stdout=original_stdout
   elif arg == "-steered":
      print(" An INCAR for steered molecular dynamics (ab inito) is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a ab-initio steered molecular dynamics calculation")
         header=header.replace('ISYM = 2 ', '# ISYM = 2 ')
         header=header.replace('# ISYM = 0 ', 'ISYM = 0 ')
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

# Strengths of the umbrella force constant(s) in eV/Ang, one for each constraint    
SPRING_K = [value1] [value2] ...

# Initial positions of umbrella potentials, for each chosen coordinate (Ang or rad)
SPRING_R0 = [value1] [value2] ...               

# The movement speed(s) of umbrella positions, in Ang/fs or rad/fs      
SPRING_V0 = [value1] [value2] ...
         ''')
      sys.stdout=original_stdout   
   elif arg == "-mlff_train":
      print(" An INCAR for on the fly machine-learning force field training.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is an on-the-fly machine-learning force field training")
         header=header.replace('ISYM = 2 ', '# ISYM = 2 ')
         header=header.replace('# ISYM = 0 ', 'ISYM = 0 ')
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

# Activates the ML-FF mode (master keyword)
ML_LMLFF = .TRUE.

# Activates the learning, new (no ML_AB) or continuation (ML_AB there)
ML_MODE = train

# Number of basis functions (local reference configurations) per element
ML_MB = 4000  # increase if complex system, decrease if memory problems
 
# Number of collected reference configurations (full system screenshots)
ML_MCONF = 4000  # increase if complex system, decrease if memory problems

# Calculation will not be aborted if ML-MCONF is reached, but configs. are discarded
ML_LBASIS_DISCARD = .TRUE.

# Cutoff for the distance (two-body) descriptor in Angstroms
ML_RCUT1 = 7.0 # might be changed, default is 5.0 

# Cutoff for the angle (three-body) descriptor in Angstroms
ML_RCUT2 = 7.0 # might be changed, default is 5.0
         ''')
      sys.stdout=original_stdout
   elif arg == "-mlff_select":
      print(" An INCAR for a selection of ML-FF training data is written.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a reselection of ML-FF training data")
         print('''
# Activates the ML-FF mode (master keyword)
ML_LMLFF = .TRUE.

# Activates the selection of training data (from ML_AB file) mode
ML_MODE = select

# Number of basis functions (local reference configurations) per element
ML_MB = 4000  # increase if complex system, decrease if memory problems
 
# Number of collected reference configurations (full system screenshots)
ML_MCONF = 4000  # increase if complex system, decrease if memory problems

# Cutoff for the distance (two-body) descriptor in Angstroms
ML_RCUT1 = 7.0 # might be changed, default is 5.0 

# Cutoff for the angle (three-body) descriptor in Angstroms
ML_RCUT2 = 7.0 # might be changed, default is 5.0 

         ''')         
      sys.stdout=original_stdout

   elif arg == "-mlff_refit":
      print(" An INCAR for a refit of a ML_FF to the fast mode.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a refit of a ML-FF to the fast mode")
         print('''
# Activates the ML-FF mode (master keyword)
ML_LMLFF = .TRUE.

# Activates the refit mode for the current training data (from ML_AB file)
ML_MODE = refit

# Activates the fast mode, the refitted ML-FF will be fast
ML_LFAST = .TRUE. 

# If the ML_AB file originates from the mlff_select program, take all configs.
# ML_EPS_LOW = 1E-20               
         ''')
      sys.stdout=original_stdout
   elif arg == "-mlff_md":
      print(" An INCAR for a MD simulation with an ML-FF.")
      original_stdout=sys.stdout
      with incar as outfile:
         sys.stdout = outfile
         print("# This is a MD simulation with a ML-FF")
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

# How new orbitals are predicted from old ones (simple extrapolation for ML)
IWAVPR = 11

# Only for NpT: Friction coefficients of Langevin thermostat, one value per element!
# LANGEVIN_GAMMA = 5.0 5.0 (...)

# Only for NpT: Friction coefficient for lattics degree of freedom
# LANGEVIN_GAMMA_L = 5.0

# Only for NpT: Ficticious mass of lattice degree of freedom
# PMASS = 5           

# Activates the ML-FF mode (master keyword)
ML_LMLFF = .TRUE.

# Activates the run mode, evaluation of energies and gradients
ML_MODE = run

# Activates the fast mode, the refitted ML-FF will be fast
ML_LFAST = .TRUE.

# Always needed for fast ML-FFs, else, an error is thrown!
ML_IALGO_LINREG = 4

# Reduced output in each time step (e.g., no radial distribution function)
ML_OUTPUT_MODE = 0

# Only each Nth frame is printed out to OUTCAR and XDATCAR
# ML_OUTBLOCK = 20 # comment in if you don't need each step               
         ''')
      sys.stdout=original_stdout



   else: 
      print(" Please give one of the options for INCAR!")
      sys.exit(1)


print(" ")
print(" gen_incar.py exited normally.")
print(" ")
