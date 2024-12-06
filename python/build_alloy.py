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

print('''
This script manages the buildup of unary/binary/ternary metal 
alloy systems in a randomized manner, i.e., the number 
of element atoms representing the overall composition
is placed randomically in a regular cubic grid.

The script must be called with a number of command line parameters.
 -elnum=[number] : The number of different elements in the system (max: 4)
 -unit_num=[number] : Number of atoms in the grid along x/y axes
 -el_symbols=[list] : List of element symbols, must agree with elnum
 -natoms=[list] : List of number of atoms per element given abive
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
      if param == "-elnum":
         elnum=int(actval)
      if param == "-unit_num":
         unit_num=int(actval)
      if param == "-el_symbols":
         el_sym_list=actval.split(",")
      if param == "-natoms":
         el_num_list=actval.split(",")
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
         
if elnum <= 0:
   print("At least one element must be given! (elnum > 0)")
   sys.exit(1)

if elnum >= 5:
   print("No more than four elements can be given! (elnum < 5)")
   sys.exit(1)

if unit_num <= 0:
   print("At least one unit must be given! (unit_num > 0)")
   sys.exit(1)

if len(el_sym_list) != elnum:
   print("Please give a valid list of element symbols! (el_symbols)")
   sys.exit(1)
 
if len(el_num_list) != elnum:           
   print("Please give a valid list of element atom numbers! (natoms)")
   sys.exit(1)


print("Given settings:")
print(" - Number of elements (-elnum):              ",str(elnum))
print(" - Number of units along x/y (-unit_num):    ",str(unit_num))
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
      print("SCALMS system: " + el1 + str(natoms1))
   if elnum == 2:
      print("SCALMS system: " + el1 + str(natoms1) + el2 + str(natoms2)) 
   if elnum == 3:
      print("SCALMS system: " + el1 + str(natoms1) + el2 + str(natoms2) + el3 +str(natoms3)) 
   if elnum == 4:
      print("SCALMS system: " + el1 + str(natoms1) + el2 + str(natoms2) + el3 +str(natoms3) + el4 +str(natoms4))
      
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

