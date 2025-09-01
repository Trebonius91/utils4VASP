#!/usr/bin/env python3
#
#    gen_poscar: Build VASP POSCAR files for bulk and surface slab
#      systems with different unit cell geometries.
#    Part of VASP4CLINT
#     Julien Steffen, 2025 (julien.steffen@fau.de)
#

import sys
import os
import numpy as np
import re
import math
import random

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
gen_poscar.py: Manages the buildup of POSCAR files for bulk and surface
systems with different unit cell geometries.
An arbitrary number of elements can be added, and the atoms can
be placed either regularly or by chance. Different crystal lattices
and surface facets can be chosen.

The script must be called with a number of command line parameters.
 -overview : Print overview of utils4VASP scripts/programs
 -lattice=[option] : Which kind of crystal lattice, available are:
    fcc (face-centered cubic)
    hcp (hexagonal close package) 
    sc  (simple cubic).
 -facet=[option] : If a surface slab shall be build: which facet
   shall be used, available are: 
    111, 100, 110  (for fcc), 
    0001, 10-10    (for hcp)
    100, 110       (for sc)
   If not given, the first item in the respective list will be chosen.
 -el_symbols=[list] : List of element symbols for the contained elements
 -el_freq=[list] : Relative frequency of the elements in the list 
   (arbitrary numbers are possible, internally normalized)
 -lat_const=[value] : Lattice constant of the system (default: 3 Ang.)  
 -unit_x=[number] : Number of atoms along the x axis (or a axis)
 -unit_y=[number] : Number of atoms along the y axis (or b acis)
 -unit_z=[number] : Number of atoms along the z axis (or c axis)
The following arguments are optional and have default values 
 -z_vac=[value] : Size of the vacuum along z-axis (default: 15 Ang.)
 -x_vac=[value] : Size of the vacuum along x-axis (default: 0 Ang.)
 -y_vac=[value] : Size of the vacuum along y-axis (default: 0 Ang.)
 -z_shift=[value] : Shift of the system along x (default: 0 Ang.)
 -x_shift=[value] : Shift of the system along x (default: 0 Ang.)
 -y_shift=[value] : Shift of the system along x (default: 0 Ang.)

Example: gen_poscar-py gen_poscar.py -lattice=sc -facet=110 -unit_x=3 
      -unit_y=3 -unit_z=6 -el_symbols=Au,Ag,In -el_freq=0.4,0.3,0.3
                 ''')

#
#    Default values for central command line parameters
#
geom_name="none"
facet_name="default"
unit_x = 1
unit_y = 1
unit_z = 1
#
#    Default values for optional command line parameters
#

unit_len = 3.00  # The lattice constant of a unit cell (Angstrom)
z_vacuum = 15.0  # Length of vacuum in z-direction (Angstrom)
x_vacuum = 0.0   # Length of vacuum in x-direction (Angstrom)
y_vacuum = 0.0   # Length of vacuum in y-direction (Angstrom)
x_shift = 0.0    # Shift of the system along x-direction (Angstrom)
y_shift = 0.0    # Shift of the system along y-direction (Angstrom)
z_shift = 0.0    # Shift of the system along z-direction (Angstrom)
#    Offset of atom positions in the grid: in the center of their "cells"
offset = unit_len/2.0

elnum=0
unit_num=0
el_sym_list=""
el_num_list=""
el_freq_list=[]

for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-lattice":
         geom_name=actval
      if param == "-facet":
         facet_name=actval
      if param == "-el_symbols":
         el_sym_list=actval.split(",")
      if param == "-el_freq":
         el_freq_list=float_arr = list(map(float, actval.split(",")))
      if param == "-unit_x":
         unit_x=int(actval)
      if param == "-unit_y":
         unit_y=int(actval)
      if param == "-unit_z":
         unit_z=int(actval)
      if param == "-z_vac":
         z_vacuum=float(actval)
      if param == "-x_vac":
         x_vacuum=float(actval)
      if param == "-y_vac":
         y_vacuum=float(actval)
      if param == "-x_shift":
         x_shift=float(actval)
      if param == "-y_shift":
         y_shift=float(actval)
      if param == "-z_shift":
         z_shift=float(actval)         
      if param == "-lat_const":
         unit_len=float(actval) 
         

#
#    Determine the number of atoms per element from the frequency
#
elnum=len(el_sym_list)
if elnum == 0:
   print("Please give at least one element with -el_symbols!\n")
   sys.exit(1)
natoms=unit_x*unit_y*unit_z
if elnum == 1:
   el_freq_list=[]
   el_freq_list.append("1.0")
else:
   norm=(sum(el_freq_list))
   if norm == 0:
      print("Plese give the relative element frequencies with -el_freq!\n")
      sys.exit(1)
   el_freq_list = np.array(el_freq_list) / norm
 
el_num_list=[]
for i in range(elnum):
   try:
      el_num_list.append(int(float(el_freq_list[i])*natoms))
   except Exception:
      print("Please give a full list of element frequencies!\n")
      sys.exit(1)
   if el_num_list[i] == 0:
      print("Species",i+1," has no atoms!\n")
      sys.exit(1)

while sum(el_num_list) < natoms:
   el_num_list[el_num_list.index(max(el_num_list))]=el_num_list[el_num_list.index(max(el_num_list))]+1
while sum(el_num_list) > natoms:
   el_num_list[el_num_list.index(max(el_num_list))]=el_num_list[el_num_list.index(max(el_num_list))]-1

#
#    Define the geometry and facet from the input and default values
#
if geom_name == "fcc":
   if facet_name == "111":
      fac_name = "111"
   elif facet_name == "100":
      fac_name = "100"
   elif facet_name == "110":
      fac_name = "110"
   elif facet_name == "default":
      fac_name = "111"
   else: 
      print("Please give a valid surface facet! (-facet)\n")
      sys.exit(1)
elif geom_name == "hcp":
   if facet_name == "0001":
      fac_name = "0001"
   elif facet_name == "10-10":
      fac_name = "10-10"
#   elif facet_name == "10-11":
#      fac_name = "10-11"
   elif facet_name == "default":
      fac_name = "0001"
   else:
      print("Please give a valid surface facet! (-facet)")
      sys.exit(1)
elif geom_name == "sc":
   if facet_name == "100":
      fac_name = "100"
   elif facet_name == "110":
      fac_name = "110"
#   elif facet_name == "111":
#      fac_name = "111"
   elif facet_name == "default":
      fac_name = "100"
   else:
      print("Please give a valid surface facet! (-facet)")
      sys.exit(1)
else:
   print("Please give a valid crystal lattice! (-lattice)")
   sys.exit(1)      


print("Given settings:")
print(" - Lattice type of the unit cell: ",geom_name)
print(" - Chosen surface facet (z-direction): ",fac_name)
print(" - Number of elements (-elnum):              ",str(elnum))
print(" - Number of atoms along x/a axis (-unit_x)  ",str(unit_x))
print(" - Number of atoms along y/b axis (-unit_y)  ",str(unit_y))
print(" - Number of atoms along z/c axis (-unit_z)  ",str(unit_z))
print(" - Vacuum along x-axis (-z_vac):             ",str(x_vacuum))
print(" - Vacuum along y-axis (-y_vac):             ",str(y_vacuum))
print(" - Vacuum along z-axis (-z_vac):             ",str(z_vacuum))
print(" - Length of a unit cell/Ang. (-unit_len):   ",str(unit_len))


print(" - Total number of atoms:                    ",str(natoms))
for i in range(elnum):
   print(" - Element ",i+1,":                              ",str(el_sym_list[i]),"  (",str(el_num_list[i]),"atoms)")


#
#    Set the unit cell size and shape of the POSCAR depending 
#    on the used geometry and facet
#
a_vec=np.zeros((3))
b_vec=np.zeros((3))
c_vec=np.zeros((3))

if geom_name == "fcc":
   if fac_name == "111":
      a_vec[0]=unit_x*unit_len
      b_vec[0]=(unit_y*unit_len+x_vacuum)/2.0
      b_vec[1]=(unit_y*unit_len*np.sqrt(3)+y_vacuum)/2.0   
      c_vec[2]=unit_z*np.sqrt(2.0)/np.sqrt(3.0)*unit_len+z_vacuum
   elif fac_name == "100":
      a_vec[0]=unit_x*unit_len+x_vacuum
      b_vec[1]=unit_y*unit_len+y_vacuum
      c_vec[2]=unit_z/np.sqrt(2.0)*unit_len+z_vacuum
   elif fac_name == "110":
      a_vec[0]=unit_x*unit_len*np.sqrt(2.0)+x_vacuum
      b_vec[1]=unit_y*unit_len+y_vacuum
      c_vec[2]=unit_z/2.0*unit_len+z_vacuum
elif geom_name == "hcp":
   if fac_name == "0001":
      a_vec[0]=unit_x*unit_len+x_vacuum
      b_vec[0]=(unit_y*unit_len+x_vacuum)/2.0
      b_vec[1]=(unit_y*unit_len*np.sqrt(3)+y_vacuum)/2.0
      c_vec[2]=unit_z*np.sqrt(2.0)/np.sqrt(3.0)*unit_len+z_vacuum
   elif fac_name == "10-10":
      a_vec[0]=unit_x*unit_len+x_vacuum
      b_vec[1]=unit_y*unit_len*1.584+y_vacuum
      c_vec[2]=unit_z*0.5*0.866025403*unit_len+z_vacuum
      if unit_z % 2 != 0:
         print(" Warning: A complete hcp 10-10 unit cell comprises of")
         print(" two primitive z-units in the scope of this script!")
elif geom_name == "sc":
   if fac_name == "100":
      a_vec[0]=unit_x*unit_len+x_vacuum
      b_vec[1]=unit_y*unit_len+y_vacuum
      c_vec[2]=unit_z*unit_len+z_vacuum
   elif fac_name == "110":
      a_vec[0]=unit_x*unit_len+x_vacuum
      b_vec[1]=unit_y*unit_len+y_vacuum
      c_vec[2]=unit_z*unit_len+z_vacuum


x_len=np.linalg.norm(a_vec)
y_len=np.linalg.norm(b_vec)
z_len=np.linalg.norm(c_vec)

offset = offset/z_len

print(" - Length of the x-axis (Ang.):              ",str(x_len))
print(" - Length of the y-axis (Ang.):              ",str(y_len))
print(" - Length of the z-axis (Ang.):              ",str(z_len))
print(" - Shift along x-axis (-x_shift):            ",str(x_shift))
print(" - Shift along y-axis (-y_shift):            ",str(y_shift))
print(" - Shift along z-axis (-z_shift):            ",str(z_shift))

#
#    Assign an element species to each atom in the system, which are 
#    thought to be in a 1D list, which will later be combined with 
#    the positions in the unit cell grid
#

el_list = list(range(natoms)) 
#
#    Create assigned list with element species
#
el_labels = np.repeat(el_sym_list,el_num_list)

#
#    Shuffle the list randomly
#
random.shuffle(el_labels)
#
#    Pair objects with species (might be needed sometimes)
#
assignment = list(zip(el_list, el_labels))


#
#    Now build the positions of the atoms, depending on the 
#    used unit cell shape, in internal coordinates first
#

pos_all=np.zeros((natoms,3))
pos_act=np.zeros((3))
inc=0
if geom_name == "fcc":
   if fac_name == "111":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 1:
                   pos_act[0]=(0.33333333333+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.33333333333+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 2:
                   pos_act[0]=(0.66666666666+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.66666666666+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 2:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1
   elif fac_name == "100":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 1:
                   pos_act[0]=(0.5+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.5+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 1:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1
   elif fac_name == "110":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(0.5+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.5+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 1:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 1:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1
if geom_name == "hcp":
   if fac_name == "0001":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(-0.3333333333+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.66666666666+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 1:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 1:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1
   elif fac_name == "10-10":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.1666666+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 1:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.5+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(-0.1666666+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 2:
                   pos_act[0]=(0.5+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.1666666+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 3:
                   pos_act[0]=(0.5+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.5+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(-0.1666666+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 3:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1
if geom_name == "sc":
   if fac_name == "100":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 0:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1
   elif fac_name == "110":
      for i in range(unit_x):
         for j in range(unit_y):
            z_count=0
            for k in range (unit_z):
                if z_count == 0:
                   pos_act[0]=(0.0+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                elif z_count == 1:
                   pos_act[0]=(0.5+i)/unit_x*(x_len-x_vacuum)/x_len
                   pos_act[1]=(0.0+j)/unit_y*(y_len-y_vacuum)/y_len
                   pos_act[2]=(0.0+k)/unit_z*(z_len-z_vacuum)/z_len+offset
                z_count=z_count+1
                if z_count > 1:
                   z_count=0
                pos_all[inc,:]=pos_act
                inc=inc+1



#
#    Now write the actual POSCAR!
#
original_stdout=sys.stdout
with open("POSCAR","w") as f:
   sys.stdout = f
   print("New POSCAR written by gen_poscar (",geom_name,",",fac_name,"facet(z-direction))")
   print("1.0")
   print(a_vec[0],a_vec[1],a_vec[2])
   print(b_vec[0],b_vec[1],b_vec[2])
   print(c_vec[0],c_vec[1],c_vec[2])   
   print(" ".join(el_sym_list))
   print(" ".join(map(str,el_num_list)))
   print("Direct")
   if elnum == 1:
      for i in range(natoms):
         print(pos_all[i][0], pos_all[i][1], pos_all[i][2])
   else:
      for i in range(elnum):
         for j in range(natoms):
            if assignment[j][1] == el_sym_list[i]:
               print(pos_all[j][0], pos_all[j][1], pos_all[j][2])

sys.stdout = original_stdout

print("""
gen_poscar.py has finished! A POSCAR file with your system was written!
""")

