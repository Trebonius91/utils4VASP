#!/usr/bin/env python3
#
#    build_alloy: Build VASP POSCAR files with atoms positioned on 
#      a regular cubic grid, designed for simulation of metal alloy
#      systems (bulk or surface slab)
#    Part of VASP4CLINT
#     Julien Steffen, 2023 (julien.steffen@fau.de)
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
 - analyze_poscar.py: Analyze POSCAR, generate KPOINTS and POTCAR
 - build_alloy.py: Build regular alloy structures of various shapes
 - modify_poscar.py: Modify POSCAR: multiply, transform, shift, insert etc.
 - cut_unitcell: Generate surface slab unit cells of arbitrary shape
 - build_adsorb.py: Place adsorbates on surfaces, set translation, rotation
 - split_freq: Divide frequency calculations for large molecules
Evaluation:
 - modify_xdatcar: Modify trajectory files: multiply, shift, pick etc.
 - analyze_bulk: Analyze bulk MD trajectories for RDFs, diffusion etc.
 - analyze_slab: Analyze slab MD trajectories for RDFs, density etc.
 - check_geoopt.py: Monitor geometry optimizations with selective dynamics
 - manage_cls: Prepare, evaluate core level energy calculations
 - eval_bader: Evaluate and visualize Bader charge calculations
 - eval_stm: Generate STM images with different settings
 - partial_dos: Plot atom and orbital resolved density of states
 - manage_neb.py: Setup, monitor and restart NEB calculations
ML-FF:
 - mlff_select: Heuristic selection of local reference configurations
 - eval_vasp_ml.py: Visualize results of VASP ML-FF on the fly learnings
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
This script manages the buildup of POSCAR files for bulk and surface
systems with different unit cell geometries.
An arbitrary number of elements can be added, and the atoms can
be placed either regularly or by chance. Simple cubic grids as well
as well-known unit cells (hpc, fcc) and surface configurations 
(001, 011, 111, 0001) can be generated.

The script must be called with a number of command line parameters.
 -overview : Print overview of utils4VASP scripts/programs
 -geometry=[option] : Which kind of close-packed structure shall be 
   taken. Currently available are: fcc (face-centered cubic), hcp 
   (hexagonal close pakackage) and sc (simple cubic).
 -facet=[option] : If a surface slab shall be build: which facet
   shall be used, available are: 111, 100, 110 (for fcc), 
   0001, 10-10, 10-11 (hcp), 100, 110, 111 (sc). If this option is 
   not given, always the first of the options for the package 
   is chosen. 
 -el_symbols=[list] : List of element symbols for the contained elements
 -elnum=[number] : The number of different elements in the system (max: 4)
 -unit_x=[number] : Number of atoms along the x axis (or a axis)
 -unit_y=[number] : Number of atoms along the y axis (or b acis)
 -unit_z=[number] : Number of atoms along the z axis (or c axis)
 -natoms=[list] : List of abundancies of elements, one number for each
The following arguments are optional and have default values 
 -z_vac=[value] : Size of the vacuum along z-axis (default: 20 Ang.)
 -x_vac=[value] : Size of the vacuum along x-axis (default: 0 Ang.)
 -y_vac=[value] : Size of the vacuum along y-axis (default: 0 Ang.)
 -z_shift=[value] : Shift of the system along x (default: 0 Ang.)
 -x_shift=[value] : Shift of the system along x (default: 0 Ang.)
 -y_shift=[value] : Shift of the system along x (default: 0 Ang.)
 -unit_len=[value]: Distance between two atoms (defaut: 2.54 Ang.)

Example: build_alloy.py -elnum=2 -unit_num=5 -el_symbols=Ga,Pt 
                 -natoms=171,9 -zvac=25.0
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

unit_len = 2.54  # The dimensions of a unit cell in x and y 
z_vacuum = 20.0  # Length of vacuum in z-direction (Angstrom)
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


for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-geometry":
         geom_name=actval
      if param == "-facet":
         facet_name=actval
      if param == "-el_symbols":
         el_sym_list=actval.split(",")
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
      if param == "-unit_len":
         unit_len=float(actval) 
         
#if elnum <= 0:
#   print("At least one element must be given! (elnum > 0)")
#   sys.exit(1)

#if elnum >= 5:
#   print("No more than four elements can be given! (elnum < 5)")
#   sys.exit(1)


#if len(el_sym_list) != elnum:
#   print("Please give a valid list of element symbols! (el_symbols)")
#   sys.exit(1)
 
#if len(el_num_list) != elnum:           
#   print("Please give a valid list of element atom numbers! (natoms)")
#   sys.exit(1)

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
      print("Please give a valid surface facet! (-facet)")
      sys.exit(1)
elif geom_name == "hcp":
   if facet_name == "0001":
      fac_name = "0001"
   elif facet_name == "10-10":
      fac_name = "10-10"
   elif facet_name == "10-11":
      fac_name = "10-11"
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
   elif facet_name == "111":
      fac_name = "111"
   elif facet_name == "default":
      fac_name = "100"
   else:
      print("Please give a valid surface facet! (-facet)")
      sys.exit(1)
else:
   print("Please give a valid unit cell geometry! (-geometry)")
   sys.exit(1)      


print("Given settings:")
print(" - Geometry of the unit cell: ",geom_name)
print(" - Chosen surface facet (z-direction): ",fac_name)
print(" - Number of elements (-elnum):              ",str(elnum))
print(" - Number of atoms along x/a axis (-unit_x)  ",str(unit_x))
print(" - Number of atoms along y/b axis (-unit_y)  ",str(unit_y))
print(" - Number of atoms along z/c axis (-unit_z)  ",str(unit_z))
print(" - Vacuum along x-axis (-z_vac):             ",str(x_vacuum))
print(" - Vacuum along y-axis (-y_vac):             ",str(y_vacuum))
print(" - Vacuum along z-axis (-z_vac):             ",str(z_vacuum))
print(" - Length of a unit cell/Ang. (-unit_len):   ",str(unit_len))

#
#    Set element symbols and numbers, calculate total number of atoms
#
el1=el_sym_list[0]
natoms1=int(el_num_list[0])
natoms=natoms1
if (elnum > 1):
   el2=el_sym_list[1]
   natoms2=int(el_num_list[1])
   natoms=natoms1+natoms2
if (elnum > 2):
   el3=el_sym_list[2]
   natoms3=int(el_num_list[2])
   natoms=natoms1+natoms2+natoms3
if (elnum > 3):
   el4=el_sym_list[3]
   natoms4=int(el_num_list[3])
   natoms=natoms1+natoms2+natoms3+natoms4


print(" - Total number of atoms:                    ",str(natoms))
print(" - Element 1:                                ",str(el1), "  (",str(natoms1)," atoms)")
if elnum > 1:
   print(" - Element 2:                                ",str(el2), "  (",str(natoms2)," atoms)")
if elnum > 2:
   print(" - Element 3:                                ",str(el3), "  (",str(natoms3)," atoms)")
if elnum > 3:
   print(" - Element 4:                                ",str(el4), "  (",str(natoms4)," atoms)")


#
#    Number of z-repetitions: ceil of total natoms divided by number of atoms 
#    per z-layer
#
z_layers = int(natoms/(unit_num**2))
if z_layers*unit_num**2 < natoms:
   z_layers= z_layers + 1

print (" - Number of atoms in x dimension:           ",str(unit_num))
print (" - Number of atoms in y dimension:           ",str(unit_num))
print (" - Number of atoms in z dimension:           ",str(z_layers))


xlen = unit_len * unit_num + x_vacuum
ylen = unit_len * unit_num + y_vacuum
zlen = unit_len * z_layers + z_vacuum


print(" - Length of the x-axis (Ang.):              ",str(xlen))
print(" - Length of the y-axis (Ang.):              ",str(ylen))
print(" - Length of the z-axis (Ang.):              ",str(zlen))
print(" - Shift along x-axis (-x_shift):            ",str(x_shift))
print(" - Shift along y-axis (-y_shift):            ",str(y_shift))
print(" - Shift along z-axis (-z_shift):            ",str(z_shift))
#
#    Build up the coordinates of the atoms 
#
at1_xyz = np.zeros((3,natoms1))
if elnum > 1:
   at2_xyz = np.zeros((3,natoms2))
   at1_frac = float(natoms1)/float(natoms)
if elnum > 2:
   at2_frac = float(natoms1+natoms2)/float(natoms)
   at3_xyz = np.zeros((3,natoms3))  
if elnum > 3:
   at3_frac = float(natoms1+natoms2+natoms3)/float(natoms)
   at4_xyz = np.zeros((3,natoms4))
   

done = False
at1_act = 0
at2_act = 0
at3_act = 0
at4_act = 0
for i in range (z_layers):  # z axis
   for j in range (unit_num):  # x-axis
      for k in range (unit_num): # y-axis
         pos_act = [offset+unit_len*j+x_shift,offset+unit_len*k+y_shift,
                  offset+unit_len*i+z_shift] 
         randnum=random.random()
         if at1_act+at2_act+at3_act >= natoms:
             done = True
             break
#
#    If ony one element is there, the buildup is really easy
#
         if (elnum == 1):
            at1_act = at1_act + 1  
            at1_xyz[0][at1_act-1] = pos_act[0]
            at1_xyz[1][at1_act-1] = pos_act[1]
            at1_xyz[2][at1_act-1] = pos_act[2]
 
#
#    If the random number is in the lower part of the interval, place 
#    a Ga, atom, unless already all have been built
#
         if (elnum == 2):
            if randnum <= at1_frac:
               at1_act = at1_act + 1
               if at1_act <= natoms1:
                  at1_xyz[0][at1_act-1] = pos_act[0]
                  at1_xyz[1][at1_act-1] = pos_act[1]
                  at1_xyz[2][at1_act-1] = pos_act[2]
               else:
                  at2_act = at2_act + 1
                  at1_act = at1_act - 1
                  at2_xyz[0][at2_act-1] = pos_act[0]
                  at2_xyz[1][at2_act-1] = pos_act[1]
                  at2_xyz[2][at2_act-1] = pos_act[2]
            else:    
               at2_act = at2_act + 1
               if at2_act <= natoms2:
                  at2_xyz[0][at2_act-1] = pos_act[0]
                  at2_xyz[1][at2_act-1] = pos_act[1]
                  at2_xyz[2][at2_act-1] = pos_act[2]
               else:
                  at1_act = at1_act + 1
                  at2_act = at2_act - 1
                  at1_xyz[0][at1_act-1] = pos_act[0]
                  at1_xyz[1][at1_act-1] = pos_act[1]
                  at1_xyz[2][at1_act-1] = pos_act[2]
         elif (elnum == 3):    
            if randnum <= at1_frac:
               at1_act = at1_act + 1
               if at1_act <= natoms1:
                  at1_xyz[0][at1_act-1] = pos_act[0]
                  at1_xyz[1][at1_act-1] = pos_act[1]
                  at1_xyz[2][at1_act-1] = pos_act[2]
               else:
                  if at3_act < natoms3:
                     at3_act = at3_act + 1
                     at1_act = at1_act - 1
                     at3_xyz[0][at3_act-1] = pos_act[0]
                     at3_xyz[1][at3_act-1] = pos_act[1]
                     at3_xyz[2][at3_act-1] = pos_act[2]
                  else:  
                     at2_act = at2_act + 1
                     at1_act = at1_act - 1
                     at2_xyz[0][at2_act-1] = pos_act[0]
                     at2_xyz[1][at2_act-1] = pos_act[1]
                     at2_xyz[2][at2_act-1] = pos_act[2]

            elif randnum <= at2_frac:
               at2_act = at2_act + 1
               if at2_act <= natoms2:
                  at2_xyz[0][at2_act-1] = pos_act[0]
                  at2_xyz[1][at2_act-1] = pos_act[1]
                  at2_xyz[2][at2_act-1] = pos_act[2]
               else:
                  if at1_act < natoms1: 
                     at1_act = at1_act + 1
                     at2_act = at2_act - 1
                     at1_xyz[0][at1_act-1] = pos_act[0]
                     at1_xyz[1][at1_act-1] = pos_act[1]
                     at1_xyz[2][at1_act-1] = pos_act[2]
                  else:
                     at3_act = at3_act + 1
                     at2_act = at2_act - 1
                     at3_xyz[0][at3_act-1] = pos_act[0]
                     at3_xyz[1][at3_act-1] = pos_act[1]
                     at3_xyz[2][at3_act-1] = pos_act[2]
            else:
               at3_act = at3_act + 1
               if at3_act <= natoms3:
                  at3_xyz[0][at3_act-1] = pos_act[0]
                  at3_xyz[1][at3_act-1] = pos_act[1]
                  at3_xyz[2][at3_act-1] = pos_act[2]
               else:
                  if at1_act < natoms1:
                     at1_act = at1_act + 1
                     at3_act = at3_act - 1
                     at1_xyz[0][at1_act-1] = pos_act[0]
                     at1_xyz[1][at1_act-1] = pos_act[1]
                     at1_xyz[2][at1_act-1] = pos_act[2]
                  else:
                     at2_act = at2_act + 1
                     at3_act = at3_act - 1
                     at2_xyz[0][at2_act-1] = pos_act[0]
                     at2_xyz[1][at2_act-1] = pos_act[1]
                     at2_xyz[2][at2_act-1] = pos_act[2]
         elif (elnum == 4):
            if randnum <= at1_frac:
               at1_act = at1_act + 1
               if at1_act <= natoms1:
                  at1_xyz[0][at1_act-1] = pos_act[0]
                  at1_xyz[1][at1_act-1] = pos_act[1]
                  at1_xyz[2][at1_act-1] = pos_act[2]
               else:
                  if at4_act < natoms4:
                     at4_act = at4_act + 1
                     at1_act = at1_act - 1
                     at4_xyz[0][at4_act-1] = pos_act[0]
                     at4_xyz[1][at4_act-1] = pos_act[1]
                     at4_xyz[2][at4_act-1] = pos_act[2]
                  elif at3_act < natoms3:
                     at3_act = at3_act + 1
                     at1_act = at1_act - 1
                     at3_xyz[0][at3_act-1] = pos_act[0]
                     at3_xyz[1][at3_act-1] = pos_act[1]
                     at3_xyz[2][at3_act-1] = pos_act[2]
                  else:
                     at2_act = at2_act + 1
                     at1_act = at1_act - 1
                     at2_xyz[0][at2_act-1] = pos_act[0]
                     at2_xyz[1][at2_act-1] = pos_act[1]
                     at2_xyz[2][at2_act-1] = pos_act[2]
            elif randnum <= at2_frac:
               at2_act = at2_act + 1
               if at2_act <= natoms2:
                  at2_xyz[0][at2_act-1] = pos_act[0]
                  at2_xyz[1][at2_act-1] = pos_act[1]
                  at2_xyz[2][at2_act-1] = pos_act[2]
               else:
                  if at1_act < natoms1:
                     at1_act = at1_act + 1
                     at2_act = at2_act - 1
                     at1_xyz[0][at1_act-1] = pos_act[0]
                     at1_xyz[1][at1_act-1] = pos_act[1]
                     at1_xyz[2][at1_act-1] = pos_act[2]
                  elif at3_act < natoms3:
                     at3_act = at3_act + 1
                     at1_act = at1_act - 1
                     at3_xyz[0][at3_act-1] = pos_act[0]
                     at3_xyz[1][at3_act-1] = pos_act[1]
                     at3_xyz[2][at3_act-1] = pos_act[2]
                  else:
                     at4_act = at4_act + 1
                     at2_act = at2_act - 1
                     at4_xyz[0][at4_act-1] = pos_act[0]
                     at4_xyz[1][at4_act-1] = pos_act[1]
                     at4_xyz[2][at4_act-1] = pos_act[2]
            elif randnum <= at3_frac:
               at3_act = at3_act + 1
               if at3_act <= natoms3:
                  at3_xyz[0][at3_act-1] = pos_act[0]
                  at3_xyz[1][at3_act-1] = pos_act[1]
                  at3_xyz[2][at3_act-1] = pos_act[2]
               else:
                  if at1_act < natoms1:
                     at1_act = at1_act + 1
                     at3_act = at3_act - 1
                     at1_xyz[0][at1_act-1] = pos_act[0]
                     at1_xyz[1][at1_act-1] = pos_act[1]
                     at1_xyz[2][at1_act-1] = pos_act[2]
                  elif at2_act < natoms2:
                     at2_act = at2_act + 1
                     at3_act = at3_act - 1
                     at2_xyz[0][at2_act-1] = pos_act[0]
                     at2_xyz[1][at2_act-1] = pos_act[1]
                     at2_xyz[2][at2_act-1] = pos_act[2]
                  else:
                     at4_act = at4_act + 1
                     at3_act = at3_act - 1
                     at4_xyz[0][at4_act-1] = pos_act[0]
                     at4_xyz[1][at4_act-1] = pos_act[1]
                     at4_xyz[2][at4_act-1] = pos_act[2]
            else:
               at4_act = at4_act + 1
               if at4_act <= natoms4:
                  at4_xyz[0][at4_act-1] = pos_act[0]
                  at4_xyz[1][at4_act-1] = pos_act[1]
                  at4_xyz[2][at4_act-1] = pos_act[2]
               else:
                  if at1_act < natoms1:
                     at1_act = at1_act + 1
                     at4_act = at4_act - 1
                     at1_xyz[0][at1_act-1] = pos_act[0]
                     at1_xyz[1][at1_act-1] = pos_act[1]
                     at1_xyz[2][at1_act-1] = pos_act[2]
                  elif at2_act < natoms2:
                     at2_act = at2_act + 1
                     at4_act = at4_act - 1
                     at2_xyz[0][at2_act-1] = pos_act[0]
                     at2_xyz[1][at2_act-1] = pos_act[1]
                     at2_xyz[2][at2_act-1] = pos_act[2]
                  else:
                     at3_act = at3_act + 1
                     at4_act = at4_act - 1
                     at3_xyz[0][at3_act-1] = pos_act[0]
                     at3_xyz[1][at3_act-1] = pos_act[1]
                     at3_xyz[2][at3_act-1] = pos_act[2]
                     



      if done:
         break
   if done:
      break

#
#    Now write the new POSCAR file 
#

original_stdout=sys.stdout
with open("POSCAR","w") as f:
   sys.stdout = f

   if elnum == 1:
      print("alloy system: " + el1 + str(natoms1))
   if elnum == 2:
      print("alloy system: " + el1 + str(natoms1) + el2 + str(natoms2)) 
   if elnum == 3:
      print("alloy system: " + el1 + str(natoms1) + el2 + str(natoms2) + el3 +str(natoms3)) 
   if elnum == 4:
      print("alloy system: " + el1 + str(natoms1) + el2 + str(natoms2) + el3 +str(natoms3) + el4 +str(natoms4))
      
   print(" 1 ")
   print("  " + str(xlen) + "  0.0   0.0")
   print(" 0.0 " + str(ylen) + " 0.0")
   print(" 0.0    0.0  " + str(zlen))
   if elnum == 1:
      print("  " + el1)
      print("  " + str(natoms1))
   if elnum == 2:
      print("  " + el1 + "  " + el2)
      print("  " + str(natoms1) + "  " + str(natoms2))
   elif elnum == 3:
      print("  " + el1 + "  " + el2 + "  " + el3) 
      print("  " + str(natoms1) + "  " + str(natoms2) +  "  " + str(natoms3)) 
   elif elnum == 4:
      print("  " + el1 + "  " + el2 + "  " + el3 + "  " + el4)
      print("  " + str(natoms1) + "  " + str(natoms2) +  "  " + str(natoms3) +  "  " + str(natoms4))


   print("Cartesian")
   for i in range(natoms1):
      print("  " + str(at1_xyz[0][i]) + " " + str(at1_xyz[1][i]) + " " + str(at1_xyz[2][i]))
   if elnum > 1:   
      for i in range(natoms2):
         print("  " + str(at2_xyz[0][i]) + " " + str(at2_xyz[1][i]) + " " + str(at2_xyz[2][i])) 
   if elnum > 2:
      for i in range(natoms3):
         print("  " + str(at3_xyz[0][i]) + " " + str(at3_xyz[1][i]) + " " + str(at3_xyz[2][i]))
   if elnum > 3:
      for i in range(natoms4):
         print("  " + str(at4_xyz[0][i]) + " " + str(at4_xyz[1][i]) + " " + str(at4_xyz[2][i]))


sys.stdout = original_stdout

print("""
build_alloy.py has finished! A POSCAR file with your system was written!
""")

