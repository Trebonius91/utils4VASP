#!/usr/bin/env python
#           
#    check_geoopt: Check the progress of VASP geometry optimizations
#      for selective and nonselective dynamics, print out energies,
#      gradients and largest gradient components as well as 
#      nonconverged atoms
#    Part of VASP4CLINT
#     Julien Steffen, 2025 (julien.steffen@fau.de)
# 

#    BASED ON:
# *****************************************************************************
#   grad2[.py]
#   Description: a tool producing summarized output of VASP calculations
#
#   Copyright 2008-2012 Peter Larsson
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# *****************************************************************************

#import commands
import os
import sys
import math
import re
import subprocess



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
check_geoopt.py: evaluates the process of a VASP geometry optimization,
by returning the total energy, largest force component, volume, spin
etc. of each ionic update step and writes them to file check_geoopt.log.
For selective dynamics, only force components of activated degrees of 
freedom are considered.

Usage: check_geoopt.py [OUTCAR-file]

Tip: this script can also be used to analyze AIMD calculations by applying
this script to their OUTCAR file.

Give -overview  to print an overview of utils4VASP scripts/programs

''')

from optparse import OptionParser

def get_number_of_atoms(where):
   os.system("grep \"NIONS\" " + where + " > grep.log")
   grep_in = open("grep.log","r")

   with grep_in as infile:
      line = infile.readline()
      line = line.rstrip("\n")
      return int(line.split()[11])

#   return int(str(subprocess.getstatusoutput("grep \"NIONS\" " + where)).split()[13])

def get_ediff(where):
   os.system("grep \"  EDIFF\" " + where + " > grep.log")
   grep_in = open("grep.log","r")

   with grep_in as infile:
      line = infile.readline()
      line = line.rstrip("\n")
      return float(line.split()[2])

#   return float(str(subprocess.getstatusoutput("grep \"  EDIFF\" " + where)).split()[2])

# Some ANSI colors
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
nelmax=300
# *****************************************
#             Main program
# *****************************************
	
parser = OptionParser()
parser.add_option("-v","--verbose",dest="verbose",help="Print extra info",action="store_true")
(options,args) = parser.parse_args()

print("This script prints out energies and gradients for (selective) dynamics!")
print(" Besides the OUTCAR, a POSCAR file must be present!")


try:
   poscar = open("POSCAR","r")
except IOError:
   sys.stderr.write(FAIL)
   sys.stderr.write("There was a problem opening the POSCAR file. Does it exist at all?")
   sys.stderr.write(ENDC+"\n")
   sys.exit(1)

#   Initialize arrays with selective dynamics flags of atoms/coordinates
select_x = []
select_y = []
select_z = []
with poscar as infile:
   line= infile.readline()
   line= infile.readline()
   line= infile.readline()
   line= infile.readline()
   line= infile.readline()
   line= infile.readline()
   line= infile.readline()
   line = line.rstrip("\n")
   atnums = line.split()
   attypes = len(atnums)
   natoms = 0
   for i in atnums:
      natoms = natoms + int(i)
   line= infile.readline()
   line= infile.readline()
   for i in range(natoms):
      line= infile.readline()
      atline = line.split()
#
#     Determine the flags for selective dynamics. If they do not exist,
#     set all to T
#
      try: 
         select_x.append(atline[3])
         select_y.append(atline[4])
         select_z.append(atline[5])
      except Exception:
         select_x.append("T")
         select_y.append("T")
         select_z.append("T")

no_conv = [True]*natoms
try:
   outcar = open(args[0],"r")
except IOError:
   sys.stderr.write(FAIL)
   sys.stderr.write("There was a problem opening the OUTCAR file. Does it exist at all?")
   sys.stderr.write(ENDC+"\n")
   sys.exit(1)

if outcar != None:
   outcarfile = args[0]
   outcarlines = outcar.readlines()

	#Find max iterations
   #nelmax = int(subprocess.run("grep NELM " + outcarfile).split()[2][0:-1])
   os.system("grep NELM " + outcarfile + " > grep1.log")
   grep_in = open("grep1.log","r")
   gradcrit = 0.0

   with grep_in as infile:
      line = infile.readline()
      line = infile.readline()
      line = line.rstrip("\n")
      print(line)
#      if line.split()[0] == "NELM": 
#          nelmax = int((line.split()[2])[:-1])
#      else:
#         print(line.split())
      nelmax = 200
   os.system("grep EDIFFG " + outcarfile + " > grep.log")
   grep_in = open("grep.log","r")

   with grep_in as infile:
  #    line = infile.readline()
      line = infile.readline()
      line = line.rstrip("\n")
      if line.split()[0] == "EDIFFG":
         gradcrit = abs(float(line.split()[2]))




   natoms = get_number_of_atoms(outcarfile)
   ediff = math.log10(float(get_ediff(outcarfile)))

   re_energy = re.compile("energy  without entropy=")
   re_iteration = re.compile("Iteration")
   re_timing = re.compile("LOOP:")
   re_force = re.compile("TOTAL-FORCE")
   re_mag = re.compile("number of electron")
   re_volume = re.compile("volume of cell")
   re_forcecrit = re.compile("EDIFFG")

   lastenergy = 0.0
   energy = 0.0
   steps = 1
   iterations = 0
   cputime = 0.0
   totaltime = 0.0
   dE = 0.0
   magmom = 0.0
   spinpolarized = False
   volume = 0.0
   #average = 0.0
   #maxforce = 0.0

   i = 0
   grad_out = open("check_geoopt.log","w")
   print(" Step-No.   energy    average force     max. force     volume")
   print(" Step-No.     energy   average force   max. force      volume",file=grad_out)
   for line in outcarlines:
      if re_iteration.search(line):
         iterations = iterations + 1
      if re_force.search(line):
         # Calculate forces here...
         forces = []
         magnitudes = []
         for j in range(0,natoms):
            parts = outcarlines[i+j+2].split()
            if select_x[j] == "T":
               x = float(parts[3])
            else:
               x = 0.0
            if select_y[j] == "T":
               y = float(parts[4])
            else:
               y = 0.0 
            if select_z[j] == "T":
               z = float(parts[5])
            else:
               z = 0.0
            forces.append([x,y,z]) 
            mag_act = math.sqrt(x*x + y*y + z*z)
            if mag_act <= gradcrit:
               no_conv[j] = False
            else:
               no_conv[j] = True 

            magnitudes.append(mag_act)
         average = sum(magnitudes)/natoms
         maxforce = max(magnitudes)
         maxpos = magnitudes.index(maxforce)
      if re_mag.search(line):
         parts = line.split()
         if len(parts) > 5 and parts[0].strip() != "NELECT":
            spinpolarized = True
            magmom = float(parts[5])

      if re_timing.search(line):
         cputime = cputime + float(line.split()[6])/60.0
             				
      if re_volume.search(line):
         parts = line.split() 
         if len(parts) > 4:
            volume = float(parts[4])
      

      if re_energy.search(line):
         lastenergy = energy 
         if line.split()[4] == "=":
            energy = float(line.split()[7])
         elif line.split()[0] == "ML":
            energy = float(line.split()[8])
         else:
            energy = float(line.split()[6])

         dE = math.log10(abs(energy-lastenergy+1.0E-12))

         # Construct output string
         try:
            stepstr = str(steps).rjust(4)
            energystr = "Energy: " + ("%3.6f" % (energy)).rjust(12)

            logdestr = "Log|dE|: " + ("%1.3f" % (dE)).rjust(6)					
            iterstr = "SCF: " + ("%3i" % (iterations))
            avgfstr="Avg|F|: " + ("%2.3f" % (average)).rjust(6)
            maxfstr="Max|F|: " + ("%2.3f" % (maxforce)).rjust(6)
            timestr="Time: " + ("%3.2fm" % (cputime)).rjust(6)
            volstr="Vol.: " + ("%3.1f" % (volume)).rjust(5)
         except NameError:
            print("Cannot understand this OUTCAR file...try to read ahead")
            continue

         if iterations == nelmax:
            sys.stdout.write(FAIL)
            #print "         ^--- SCF cycle reached NELMAX. Check convergence!"

         if (dE < ediff):
            sys.stdout.write(OKGREEN)
         
         if spinpolarized:
            # sys.stdout.write(("Step %3i  Energy: %+3.6f  Log|dE|: %+1.3f  Avg|F|: %.6f  Max|F|: %.6f  SCF: %3i  Mag: %2.2f  Time: %03.2fm") % (steps,energy,dE,average,maxforce,iterations,magmom,cputime))
            magstr="Mag: " + ("%2.2f" % (magmom)).rjust(6)
            print("%s  %s  %s  %s  %s  %s  %s  %s  %s" % (stepstr,energystr,logdestr,iterstr,avgfstr,maxfstr,volstr,magstr,timestr))
            print("%4i %3.7f  %2.5f %2.5f  %3.7f" % (steps,energy,average,maxforce,volume), file=grad_out)
         else:
            print("%s  %s  %s  %s  %s  %s  %s  %s" % (stepstr,energystr,logdestr,iterstr,avgfstr,maxfstr,volstr,timestr))
            print("%4i       %3.7f   %2.5f    %2.5f   %3.7f " % (steps,energy,average,maxforce,volume), file=grad_out)
            # sys.stdout.write(("Step %3i  Energy: %+3.6f  Log|dE|: %+1.3f  Avg|F|: %.6f  Max|F|: %.6f  SCF: %3i  Time: %03.2fm") % (steps,energy,dE,average,maxforce,iterations,cputime))

         sys.stdout.write(ENDC)

         steps = steps + 1
         iterations = 0
         totaltime = totaltime + cputime
         cputime = 0.0
	
      i = i + 1
for i in range(natoms):
   if no_conv[i]:
      print("Atom " + str(i) + " is not converged yet.")
print ("The maximum force acts on atom " + str(maxpos) + ".")
print ("File 'check_geoopt.log' for plots was written!")
