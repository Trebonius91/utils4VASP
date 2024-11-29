#!/usr/bin/env python3
#   
#    modify_poscar: Modify a given POSCAR file in several different
#      ways such as: shift its content, multiply it in different 
#      coordinate directions, convert it from frac 2 cart and the 
#      other way around, write it to a xyz file
#    Part of VASP4CLINT
#     Julien Steffen, 2024 (julien.steffen@fau.de)
#

import sys
import os
import re
import numpy as np
from numpy import linalg as LA
from scipy.spatial import distance_matrix

print('''
 This script takes a POSCAR file in direct or cartesian coordinates 
 and performs several different operations on it, depending 
 on the keyword and its specifiers given in the command line. 
 The list of possible options:
  -shift=a,b,c : Shift the unitcell contents along the given components 
     of a vector in direct coordinates along the coordinate axes.
     Example: shift=0.1,0.0,0.5 
  -multiply=a,b,c: Multiply the unit cell along a, b and c axes, where 
     the arguments are integers. Example: -multiply=2,2,1 (doubling the 
     cell and its content along x and y, keeping it along z)
  -frac2cart : If the POSCAR is in direct/fractional coordinates, convert
     it to cartesian coordinates
  -cart2frac : The POSCAR is in cartesian coordinates, convert it to
     fractional coordinates 
  -map2unit : Remove all image flags (i.e., direct coordinates larger 1
     or smaller 0) and shift all atoms to the central unit cell
  -freeze=elems : Set all atoms belonging to a list of elements to F F F,
     for example for frequency calculations, where a surface is kept fix.
     Example: freeze=Pt,O (Pt and O atoms will be kept fix) 
  -select_el=elems : Select all elements whose atoms shall be modified,
     for example with a shift command. Same syntax as for freeze.
  -select_ind=numbers : Select all atoms by index (start with 1) which
     shall be modified, for example with a shift command. Both single 
     numbers and ranges can be given, for example: select_ind=1,7-13,56
  -insert_struc=a,b,c : Insert the structure given by POSCAR_insert into
     the current POSCAR. The center of mass of the inserted structure is 
     moved to the position given by a,b,c in direct coordinates.
  -bottom=value : Set all atoms below the z-coordinate 'value' to F F F,
     for example for optimizations on surfaces, where the lower two layers
     are kept fix. Depending on the input format of the POSCAR, the value
     must be given in direct or cartesian coordinates.
     Example: bottom=3.5 (atoms below 3.5 A are frozen, cartesian)
  -remove_dist=value : The distance (in Angstroms) below which all atoms 
     in the POSCAR are removed, which are nearer to one of the atoms in 
     the POSCAR_insert structure (only for insert_struc): default: 2.5 Ang
  -phasetrans=axis : a unitcell for a phase transition sampling will be 
     build. Must be done in combination with -multiply, for 'axis' either
     a, b or c must be given (as letter). Along the given axis, the 
     multiplication must be done by an even number. The atoms in the lower
     part of the cell (along the chosen axis) will be kept fix.
     Example -phasetrans=c (with -multiply=2,2,4)
  -writexyz : Write the coordinates of the POSCAR file to a xyz file  
 More than one job can be done at once, the ordering of operation is the 
 same as the ordering of keywords above 
''')

multiply_job=False
shift_job=False
phasetrans_job=False
frac2cart=False
freeze=False
cart2frac=False
writexyz=False
bottom=False
insert=False
select_el=False
select_ind=False
map2unit=False
dist_remove=2.5

# Read in the command line arguments
for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-multiply":
         multiply_list=actval.split(",")
         multiply_job=True
      if param == "-shift":
         shift_list=actval.split(",")
         shift_job=True
      if param == "-phasetrans":
         phasetrans_axis=actval.split(",") 
         phasetrans_job=True
      if param == "-freeze":
         freeze_list=actval.split(",")
         freeze=True
      if param == "-insert_struc":
         insert_list=actval.split(",")
         insert=True
      if param == "-select_el":
         select_list=actval.split(",")
         select_el=True
      if param == "-select_ind":
         select_ind_list=actval.split(",")
         select_ind=True
      if param == "-remove_dist":
         dist_remove=float(actval)
      if param == "-bottom":
         freeze_bottom=float(actval)
         bottom=True

   else:
      param=arg
      if param == "-frac2cart":
         frac2cart=True
      if param == "-cart2frac":
         cart2frac=True
      if param == "-writexyz":
         writexyz=True
      if param == "-map2unit":
         map2unit=True 

if ((not multiply_job) and (not shift_job) and (not frac2cart) and (not cart2frac) 
    and (not writexyz) and (not freeze) and (not bottom) and(not insert)):
   print("Please give a least one valid keyword!")
   sys.exit(1)

# Read in the POSCAR file
poscar_name="POSCAR"

poscar_in = open(poscar_name,"r")

# array for selective dynamics specifiers
coord_select= []  

cartesian=False
selective=False
with poscar_in as infile:
   line = infile.readline()
   line = line.rstrip("\n")
   poscar_comm = line # the comment line 
   line = infile.readline().rstrip("\n")
   sys_scale = float(line)  # the global lattice scaling factor
   # read in the lattice vectors a,b and c

   a_vec=np.zeros(3)
   b_vec=np.zeros(3)
   c_vec=np.zeros(3)

   line = infile.readline().rstrip("\n")
   a_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   b_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   c_read = line.rstrip().split()[0:3]
   for i in range(3):
      a_vec[i]=float(a_read[i])
      b_vec[i]=float(b_read[i])
      c_vec[i]=float(c_read[i])
   # read in the element ordering
   coord_vec=np.zeros(3)
   coord_vec[0]=a_vec[0]
   coord_vec[1]=b_vec[1]
   coord_vec[2]=c_vec[2]
   line = infile.readline().rstrip("\n")
   elements = line.rstrip().split()
   nelem = len(elements)
   # read in the number of elements 
   line = infile.readline().rstrip("\n")
   line_split = line.rstrip().split()

   elem_num=[]
   natoms=0
   names=[]
   for i in range(nelem):
      elem_num.append(int(line_split[i]))
      # total number of atoms
      natoms=natoms+elem_num[i]
      for j in range(elem_num[i]):
          names.append(elements[i])


   natoms=int(natoms)
   # read in the list of atoms 
   xyz = np.zeros((natoms,3))
   # check if selective dynamics was used
   line = infile.readline()
   line_split = line.rstrip().split()
   if (line_split[0] == "Selective" or line_split[0] == "selective"
        or line_split[0] == "Select" or line_split[0] == "select"):
      selective=True
      print(" The POSCAR file has selective dynamics.")
      print(" ")
   if selective:
      line = infile.readline().rstrip("\n")
   line_parts = line.split()
   if (line_parts[0] != "Direct" and line_parts[0] != "direct"  
                and line_parts[0] != " Direct" and line_parts[0] != " direct"):
      cartesian=True
      print(" The POSCAR has cartesian coordinates.")
   else:
      print(" The POSCAR has direct coordinates.")

   print(" ")
   for i in range(natoms):
      line = infile.readline().rstrip("\n")
      xyz_read = line.rstrip().split()[0:3]
      if selective:
         select_read = line.rstrip().split()[3:6]
         coord_select.append(select_read[0] + " " + select_read[1] + " " + select_read[2]) 
      for j in range(3):
         xyz[i][j]=float(xyz_read[j]) 
# If the freeze mode is activated, set the selective dynamics for all chosen elements to F
#  and for all others to T, no matter, which selectivity existed beforehand

if freeze:
   coord_select = [] 
   if not selective: 
      selective = True     
   for name in names:
      freeze_act=False  
      for frozen in freeze_list: 
         if name == frozen:
            freeze_act=True 
      if freeze_act:
         coord_select.append(" F F F ")
      else:   
         coord_select.append(" T T T ")
if ((not multiply_job) and (not shift_job) and (not frac2cart) and (not cart2frac)):
   xyz_new=xyz
   select_new=coord_select

# For elements selected for shift or other operations: generate boolean mask with 
#  indices of atoms that are one of the chosen elements

select_mask=[]
for name in names:
   select_mask.append(True)


#
#    Select by element: Set all atoms of one or several elements true
#
if select_el:
   select_mask=[] 
   for name in names:
      select_act=False
      for selected in select_list:
         if name == selected:
            select_act=True
      select_mask.append(select_act)
#
#    Select by index: Set either single atoms or whole ranges of atoms true
#
if select_ind:
   select_mask=[]
   for name in names:
      select_mask.append(False)
   for line in select_ind_list:
      multiply_list=line.split("-")   
      if len(multiply_list) == 1:
         index1=int(multiply_list[0])
         select_mask[index1-1]=True
      if len(multiply_list) == 2:
         index1=int(multiply_list[0])
         index2=int(multiply_list[1])+1
         for i in range (index1,index2):
            select_mask[i-1]=True

# Define coordinate transformation as functions in order to use them intermediately!
# F1: FRAC2CART
def trans_frac2cart(xyz,natoms,a_vec,b_vec,c_vec):
   xyz_cart = np.zeros((natoms,3))
   for i in range(natoms):
      xyz_cart[i][0]=(xyz[i][0]*a_vec[0]+xyz[i][1]*b_vec[0]+\
                     xyz[i][2]*c_vec[0])*sys_scale
      xyz_cart[i][1]=(xyz[i][0]*a_vec[1]+xyz[i][1]*b_vec[1]+\
                     xyz[i][2]*c_vec[1])*sys_scale
      xyz_cart[i][2]=(xyz[i][0]*a_vec[2]+xyz[i][1]*b_vec[2]+\
                     xyz[i][2]*c_vec[2])*sys_scale

   return xyz_cart

# F2: CART2FRAC
def trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec):
   xyz_frac = np.zeros((natoms,3))
   act_vec = np.zeros((3))
   vec_mat=np.zeros((3,3))    
#
#   Build direct from cartesian coordinates by inverting the 
#   unit cell vector matrix and multiplying it to the coordinate vector
#

   vec_mat[0][0]=a_vec[0]
   vec_mat[1][0]=a_vec[1]
   vec_mat[2][0]=a_vec[2]
   vec_mat[0][1]=b_vec[0]
   vec_mat[1][1]=b_vec[1]
   vec_mat[2][1]=b_vec[2]
   vec_mat[0][2]=c_vec[0]
   vec_mat[1][2]=c_vec[1]
   vec_mat[2][2]=c_vec[2]

   mat_inv=np.linalg.inv(vec_mat)
   
   for i in range(natoms):      
      xyz_frac[i][0]=(xyz[i][0]*mat_inv[0][0]+xyz[i][1]*mat_inv[0][1]+\
                     xyz[i][2]*mat_inv[0][2])*sys_scale
      xyz_frac[i][1]=(xyz[i][0]*mat_inv[1][0]+xyz[i][1]*mat_inv[1][1]+\
                     xyz[i][2]*mat_inv[1][2])*sys_scale
      xyz_frac[i][2]=(xyz[i][0]*mat_inv[2][0]+xyz[i][1]*mat_inv[2][1]+\
                     xyz[i][2]*mat_inv[2][2])*sys_scale

      
   return xyz_frac

# If insert is activated, map the atoms to the central unit cell
if insert:
   map2unit=True

# 0: MAP ALL ATOMS TO CENTRAL UNIT CELL ################
if map2unit:
   print(" Map all atoms to the central unit cell, i.e., remove image flags ...") 

   if selective:
      select_new=[]

   # Perform the actual shifting of the coordinates
   xyz_new=np.zeros((natoms,3))
   # If cartesian coordinates: translate the coordinates first to direct, then 
   #  perform the shift, and finally translate it back to cartesian!
   if cartesian:
      xyz = trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec)
      for i in range(natoms):
         for j in range(3):
            xyz_new[i][j] = xyz[i][j]
            while xyz_new[i][j] < 0.0:
               xyz_new[i][j] = xyz_new[i][j] + 1.0
            while xyz_new[i][j] > 1.0:
               xyz_new[i][j] = xyz_new[i][j] - 1.0
         if selective:
            select_new.append(coord_select[i])
      xyz_new = trans_frac2cart(xyz_new,natoms,a_vec,b_vec,c_vec)
   else:
      for i in range(natoms):
         for j in range(3):
            xyz_new[i][j] = xyz[i][j] 
            while xyz_new[i][j] < 0.0:
               xyz_new[i][j] = xyz_new[i][j] + 1.0
            while xyz_new[i][j] > 1.0:
               xyz_new[i][j] = xyz_new[i][j] - 1.0

         if selective:
            select_new.append(coord_select[i])
   print(" done!\n")

# A: SHIFT OF UNIT CELL ###############################
if shift_job:
   print(" Shift the unit cell along the vector: a=",shift_list[0],", b=",shift_list[1],
        " c=",shift_list[2],"  ...")
# iterating till the range
   shift_vec = []
   for i in range(3):
      ele = float(shift_list[i])
      shift_vec.append(ele) # adding the element



   if len(shift_vec) < 3:
      print("Please give three values for the coordinates!")
      exit(1)
   for i in range(3):
      if shift_vec[i] < -1.0:
         print("ERROR! All three shift coordinates must be greater than -1!")
         exit(1)
      if shift_vec[i] > 1.0:
         print("ERROR! All three shift coordinates must be smaller than 1!")
         exit(1)
   

   if selective:
      select_new=[]

   # If the unit cell was mapped, used the updated cell
   if map2unit:
      xyz=xyz_new


   # Perform the actual shifting of the coordinates
   xyz_new=np.zeros((natoms,3))
   # If cartesian coordinates: translate the coordinates first to direct, then 
   #  perform the shift, and finally translate it back to cartesian!
   if cartesian:
      xyz = trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec)
      for i in range(natoms):
         for j in range(3):
            if select_mask[i]: 
               xyz_new[i][j] = xyz[i][j] + shift_vec[j]
            else:  
               xyz_new[i][j] = xyz[i][j] 
            if xyz_new[i][j] < 0.0:
               xyz_new[i][j] = xyz_new[i][j] + 1.0
            if xyz_new[i][j] > 1.0:
               xyz_new[i][j] = xyz_new[i][j] - 1.0
         if selective:
            select_new.append(coord_select[i])
      xyz_new = trans_frac2cart(xyz_new,natoms,a_vec,b_vec,c_vec)

   else:
      for i in range(natoms):
         for j in range(3):
            if select_mask[i]: 
               xyz_new[i][j] = xyz[i][j] + shift_vec[j]
            else:   
               xyz_new[i][j] = xyz[i][j] 
            if xyz_new[i][j] < 0.0:
               xyz_new[i][j] = xyz_new[i][j] + 1.0
            if xyz_new[i][j] > 1.0:
               xyz_new[i][j] = xyz_new[i][j] - 1.0 

         if selective:
            select_new.append(coord_select[i])  
   print(" done!\n")

# B: MULTIPLY THE UNITCELL #################################
if multiply_job:
   print(" Multiply the unit cell: ",multiply_list[0],"x along a, ",multiply_list[1],
                 "x along b, ",multiply_list[2],"x along c ...") 
   mult_vec = np.zeros(3)
   for i in range(3):
      ele = int(multiply_list[i])
      mult_vec[i] = ele # adding the element
#   If the previous job (shift of unit cell) was already done, overwrite the xyz array   
   if shift_job or map2unit:
      xyz=xyz_new

#   If a phase transition shall be studied, hold the lower half of the multiplied cell
#   along the chosen axis frozen. Throw error, if the multiplication along this axis 
#   has no even value
   if phasetrans_job:
      if phasetrans_axis[0] == "a":
         if (mult_vec[0] % 2) == 0: 
            print(" A phase transition will be built, the lower ",int(mult_vec[0]/2),
                  "cells along the a-axis")
            print("  will be kept frozen.")
         else: 
            print(" Please give an even number of unit cells along the axis (a) that shall")
            print("  be used to sample the phase transition!")
            sys.exit(1)
      if phasetrans_axis[0] == "b":
         if (mult_vec[1] % 2) == 0:
            print(" A phase transition will be built, the lower ",int(mult_vec[1]/2),
                  "cells along the b-axis") 
            print("  will be kept frozen.")
         else: 
            print(" Please give an even number of unit cells along the axis (b) that shall")
            print("  be used to sample the phase transition!")
            sys.exit(1)
      if phasetrans_axis[0] == "c":
         if (mult_vec[2] % 2) == 0:
            print(" A phase transition will be built, the lower ",int(mult_vec[2]/2),
                  "cells along the c-axis") 
            print("  will be kept frozen.")
         else: 
            print(" Please give an even number of unit cells along the axis (c) that shall")
            print("  be used to sample the phase transition!")
            sys.exit(1)
            

#   First, determine the size and the number of atoms in the new unit cell      
  
   factor=int(mult_vec[0])*int(mult_vec[1])*int(mult_vec[2])
   natoms_old=natoms
   natoms=natoms*factor

   for i in range(nelem):
      elem_num[i]=elem_num[i]*factor


#  Scale the atom coordinates
   xyz_new=np.zeros((natoms,3))
#  The individual atoms: their elements   
   names_new=[]
   if selective:
      select_new=[]
   if phasetrans_job:
      select_new=[]
   if cartesian:
      pos_new=0
      pos_old=0
      for i in range(natoms_old):
         for j in range(int(mult_vec[0])):
            for k in range(int(mult_vec[1])):
                for l in range(int(mult_vec[2])):
                   xyz_new[pos_new,0]=(xyz[pos_old,0]+float(j)*a_vec[0]+float(k)*b_vec[0]+
                                       float(l)*c_vec[0])
                   xyz_new[pos_new,1]=(xyz[pos_old,1]+float(j)*a_vec[1]+float(k)*b_vec[1]+
                                       float(l)*c_vec[2])
                   xyz_new[pos_new,2]=(xyz[pos_old,2]+float(j)*a_vec[2]+float(k)*b_vec[2]+
                                       float(l)*c_vec[2])
                   pos_new=pos_new+1
#   For the case of phase transitions, freeze the lower part of the multiplied cell
                   if phasetrans_job:
                      if phasetrans_axis[0] == "a":
                         if j < int(mult_vec[0]/2):
                            select_new.append(" F F F ")
                         else:
                            select_new.append(" T T T ")
                      if phasetrans_axis[0] == "b":
                         if k < int(mult_vec[1]/2):
                            select_new.append(" F F F ")
                         else:
                            select_new.append(" T T T ")
                      if phasetrans_axis[0] == "c":
                         if l < int(mult_vec[2]/2):
                            select_new.append(" F F F ")
                         else:
                            select_new.append(" T T T ")
                   else:
                      if selective:
                         select_new.append(coord_select[pos_old])
                   names_new.append(names[pos_old])

         pos_old=pos_old+1


   else:
      pos_new=0
      pos_old=0
      for i in range(natoms_old):
         for j in range(int(mult_vec[0])):
            for k in range(int(mult_vec[1])):
                for l in range(int(mult_vec[2])):
                   xyz_new[pos_new,0]=(xyz[pos_old,0]+1.0*float(j))/float(mult_vec[0])
                   xyz_new[pos_new,1]=(xyz[pos_old,1]+1.0*float(k))/float(mult_vec[1])
                   xyz_new[pos_new,2]=(xyz[pos_old,2]+1.0*float(l))/float(mult_vec[2])
                   pos_new=pos_new+1
#   For the case of phase transitions, freeze the lower part of the multiplied cell
                   if phasetrans_job:
                      if phasetrans_axis[0] == "a":
                         if j < int(mult_vec[0]/2):
                            select_new.append(" F F F ")
                         else:
                            select_new.append(" T T T ")
                      if phasetrans_axis[0] == "b":
                         if k < int(mult_vec[1]/2):
                            select_new.append(" F F F ")            
                         else:
                            select_new.append(" T T T ")
                      if phasetrans_axis[0] == "c":
                         if l < int(mult_vec[2]/2):
                            select_new.append(" F F F ")            
                         else:
                            select_new.append(" T T T ")     
                   else:    
                      if selective:
                         select_new.append(coord_select[pos_old])
                   names_new.append(names[pos_old])
         pos_old=pos_old+1

#   Scale the unit cell vectors
   a_vec=np.multiply(a_vec,mult_vec[0])
   b_vec=np.multiply(b_vec,mult_vec[1])
   c_vec=np.multiply(c_vec,mult_vec[2])

   names=names_new
#   For a phase transition job, always print selective coordinates   
   if phasetrans_job:
      selective = True

   if selective:
      coord_select=select_new
   print(" done!\n")

# C: Insert a second structure into the coordinates of POSCAR
if insert:
   print(" Insert a structure from file POSCAR_insert into POSCAR...")
   
# Read in the POSCAR_insert file
   poscar_name="POSCAR_insert"

   poscar_insert = open(poscar_name,"r")

# array for selective dynamics specifiers
   coord_select= []

   cartesian_insert=False
   selective_insert=False
   with poscar_insert as infile:
      line = infile.readline()
      line = line.rstrip("\n")
      poscar_comm = line # the comment line
      line = infile.readline().rstrip("\n")
      sys_scale_insert = float(line)  # the global lattice scaling factor
   # read in the lattice vectors a,b and c

      a_vec_insert=np.zeros(3)
      b_vec_insert=np.zeros(3)
      c_vec_insert=np.zeros(3)

      line = infile.readline().rstrip("\n")
      a_read = line.rstrip().split()[0:3]
      line = infile.readline().rstrip("\n")
      b_read = line.rstrip().split()[0:3]
      line = infile.readline().rstrip("\n")
      c_read = line.rstrip().split()[0:3]
      for i in range(3):
         a_vec_insert[i]=float(a_read[i])
         b_vec_insert[i]=float(b_read[i])
         c_vec_insert[i]=float(c_read[i])
   # read in the element ordering
      coord_vec=np.zeros(3)
      coord_vec[0]=a_vec_insert[0]
      coord_vec[1]=b_vec_insert[1]
      coord_vec[2]=c_vec_insert[2]
      line = infile.readline().rstrip("\n")
      elements_insert = line.rstrip().split()
      nelem_insert = len(elements_insert)
   # read in the number of elements
      line = infile.readline().rstrip("\n")
      line_split = line.rstrip().split()

      elem_num_insert=[]
      natoms_insert=0
      names_insert=[]
      for i in range(nelem_insert):
         elem_num_insert.append(int(line_split[i]))
      # total number of atoms
         natoms_insert=natoms_insert+elem_num_insert[i]
         for j in range(elem_num_insert[i]):
             names_insert.append(elements_insert[i])


      natoms_insert=int(natoms_insert)
   # read in the list of atoms
      xyz_insert = np.zeros((natoms_insert,3))
   # check if selective dynamics was used
      line = infile.readline()
      line_split = line.rstrip().split()
      if (line_split[0] == "Selective" or line_split[0] == "selective"
           or line_split[0] == "Select" or line_split[0] == "select"):
         selective_insert=True
      if selective_insert:
         line = infile.readline().rstrip("\n")
      line_parts = line.split()
      if (line_parts[0] != "Direct" and line_parts[0] != "direct"
                and line_parts[0] != " Direct" and line_parts[0] != " direct"):
         cartesian_insert=True
         print(" The POSCAR_insert has cartesian coordinates.")
      else:
         print(" The POSCAR_insert has direct coordinates.")

      for i in range(natoms_insert):
         line = infile.readline().rstrip("\n")
         xyz_read = line.rstrip().split()[0:3]
         if selective:
            select_read = line.rstrip().split()[3:6]
            coord_select_insert.append(select_read[0] + " " + select_read[1] + " " + select_read[2])
         for j in range(3):
            xyz_insert[i][j]=float(xyz_read[j])

   # Shift all atoms of the inserted structure into the central unit cell
   xyz_insert_new=np.zeros((natoms_insert,3))
   # If cartesian coordinates: translate the coordinates first to direct, then 
   #  perform the shift, and finally translate it back to cartesian!
   if cartesian_insert:
      xyz_insert = trans_cart2frac(xyz_insert,natoms_insert,a_vec_insert,b_vec_insert,c_vec_insert)
      for i in range(natoms_insert):
         for j in range(3):
            xyz_insert_new[i][j] = xyz_insert[i][j]
            while xyz_insert_new[i][j] < 0.0:
               xyz_insert_new[i][j] = xyz_insert_new[i][j] + 1.0
            while xyz_insert_new[i][j] > 1.0:
               xyz_insert_new[i][j] = xyz_insert_new[i][j] - 1.0
         if selective:
            select_insert_new.append(coord_select_insert[i])
      xyz_insert_new = trans_frac2cart(xyz_insert_new,natoms_insert,a_vec_insert,b_vec_insert,c_vec_insert)
   else:
      for i in range(natoms_insert):
         for j in range(3):
            xyz_insert_new[i][j] = xyz_insert[i][j]
            while xyz_insert_new[i][j] < 0.0:
               xyz_insert_new[i][j] = xyz_insert_new[i][j] + 1.0
            while xyz_insert_new[i][j] > 1.0:
               xyz_insert_new[i][j] = xyz_insert_new[i][j] - 1.0

         if selective:
            select_insert_new.append(coord_select_insert[i])
   print(" done!\n")


#  Shift the content of POSCAR_insert by the shift vector given in the command
#  Since the shift vector is defined relative to the shape of the POSCAR in
#  which POSCAR_insert shall be inserted, it must be translated into the 
#  direct coordinates of POSCAR_insert!
   shift_outer = np.zeros((1,3))
   for i in range(3):
      shift_outer[0][i] = float(insert_list[i])

#  Translate to cartesian coordinates
   shift_cart=trans_frac2cart(shift_outer,1,a_vec,b_vec,c_vec)
#  Translate to direct coordinates of POSCAR_insert
   shift_insert=trans_cart2frac(shift_cart,1,a_vec_insert,b_vec_insert,c_vec_insert)
#  If needed, translate the coordinates of the inserted structure to direct
   if cartesian_insert:
      xyz_insert_new = trans_cart2frac(xyz_insert_new,natoms_insert,a_vec_insert,b_vec_insert,c_vec_insert)


#  We want to move the center of POSCAR insert, so correct to the center
   for i in range(3):
      shift_insert[0][i]=shift_insert[0][i]-0.5

   
#  Now shift the atoms of POSCAR_insert by the shift vector of it

   for i in range(natoms_insert):
      for j in range(3):
         xyz_insert_new[i][j] = xyz_insert_new[i][j] + shift_insert[0][j]

#  Now translate the inserted structure back to cartesian coordinates

   xyz_insert_new = trans_frac2cart(xyz_insert_new,natoms_insert,a_vec_insert,b_vec_insert,c_vec_insert)

#  If needed, translate the large POSCAR structure to cartesian coordinates
#  If one of the previous jobs were done, overwrite the xyz array
   if shift_job or multiply_job or map2unit:
      xyz=xyz_new

   if not cartesian:
      xyz=trans_frac2cart(xyz,natoms,a_vec,b_vec,c_vec)

#  Now, plug the inserted structure into the POSCAR structure 
#   Generate a distance matrix of all atoms in POSCAR_insert and POSCAR, remove all
#   atoms in POSCAR that are below a certain threshold

   dist_mat = distance_matrix(xyz,xyz_insert_new,p=2)

   dist_mins = dist_mat.min(axis=1)
 
   min_mask=dist_mins < dist_remove

#  Finally, overwrite arrays and values for POSCAR with new numbers for inserted structure  
   index=0
   index2=0
   xyz_tmp=np.zeros((natoms,3))
   elem_num_tmp=[0]*nelem
   natoms_tmp=0
   names_tmp=[""]*natoms

   for i in range (nelem):
      for j in range(elem_num[i]):
         if not min_mask[index]:
            xyz_tmp[index2][0]=xyz[index][0]
            xyz_tmp[index2][1]=xyz[index][1]
            xyz_tmp[index2][2]=xyz[index][2]
            names_tmp[index2]=names[index]
            index2=index2+1 
            elem_num_tmp[i]=elem_num_tmp[i]+1
            natoms_tmp=natoms_tmp+1
         index=index+1
#  Now, add the elements and atom numbers of POSCAR_insert to the structure
  
   for i in range (nelem_insert):
      found=False
      for j in range (nelem):
         if elements[j] == elements_insert[i]:
            elem_num_tmp[j] = elem_num_tmp[j]+elem_num_insert[i]
            found=True
            break
      if not found:
         nelem=nelem+1
         elements.append(elements_insert[i])
         elem_num_tmp.append(elem_num_insert[i])   

   natoms=natoms_tmp+natoms_insert
   xyz_new=np.zeros((natoms,3))
   names=[""]*natoms

   index=0
   for i in range(nelem):
      for j in range(natoms_tmp):
         if names_tmp[j] == elements[i]:
            names[index]=names_tmp[j]
            xyz_new[index][0]=xyz_tmp[j][0]
            xyz_new[index][1]=xyz_tmp[j][1]
            xyz_new[index][2]=xyz_tmp[j][2] 
            index=index+1 
      for j in range(natoms_insert):
         if names_insert[i] == elements[i]:
            names[index]=names_insert[j] 
            xyz_new[index][0]=xyz_insert_new[j][0]
            xyz_new[index][1]=xyz_insert_new[j][1]
            xyz_new[index][2]=xyz_insert_new[j][2]  
            index=index+1
      elem_num[i]=elem_num_tmp[i]
   cartesian=True

# D: TRANSLATE FROM DIRECT TO CARTESIAN ###################################
if frac2cart:
   print(" Translate struture from fractional/direct to cartesian coordinates...")   
   if cartesian:
      print(" The POSCAR file already has cartesian coordinates!")
      exit(1)
#  If one of the previous jobs were done, overwrite the xyz array
   if shift_job or multiply_job or map2unit:
      xyz=xyz_new 

   xyz_new=trans_frac2cart(xyz,natoms,a_vec,b_vec,c_vec)
   cartesian=True
   select_new=[]
   if selective:
      for i in range(natoms):
         select_new.append(coord_select[i])


   print(" done!\n")
  
# E: TRANSLATE FROM CARTESIAN TO DIRECT ###################################
if cart2frac:
   print(" Translate struture from cartesian to fractional/direct coordinates...")
   if not cartesian:
      print(" The POSCAR file already has direct coordinates!")
      exit(1)    
#  If one of the previous jobs were done, overwrite the xyz array
   if shift_job or multiply_job or frac2cart or map2unit:
      xyz=xyz_new


   xyz_new=trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec)

   cartesian=False
   select_new=[]
   if selective:
      for i in range(natoms):
         select_new.append(coord_select[i])


   print(" done!\n")

# F: WRITE XYZ FILE #######################################################
if writexyz:
   print(" Write structure to xyz file poscar_mod.xyz...")
#  If no other job has been done before, copy xyz directly to xyz_new
   if not shift_job and not multiply_job and not frac2cart and not cart2frac and not map2unit:
      xyz_new=xyz 
#  If the structure is in direct coordinates, translate it first to cartesian!
   if (not cartesian):
      xyz_new=trans_frac2cart(xyz_new,natoms,a_vec,b_vec,c_vec) 

#  Translate it to direct coordinates and move all atoms into unit cell
   xyz_new=trans_cart2frac(xyz_new,natoms,a_vec,b_vec,c_vec)
   
   for i in range(natoms):
      for j in range(3):
         while xyz_new[i][j] < 0.0:
            xyz_new[i][j] = xyz_new[i][j] + 1.0
         while xyz_new[i][j] > 1.0:
            xyz_new[i][j] = xyz_new[i][j] - 1.0
   
#  Translate to cartesian coordinates for writeout
   xyz_new=trans_frac2cart(xyz_new,natoms,a_vec,b_vec,c_vec)

   original_stdout=sys.stdout
   with open("poscar_mod.xyz","w") as f:
      sys.stdout = f
      print(natoms)
      print(" System converted to xyz by modify_poscar.py")
      for i in range(natoms):
         print(names[i],"  ",str(xyz_new[i][0]),"  ",str(xyz_new[i][1]),"  ",
                    str(xyz_new[i][2]))
#  Translate the structure back to direct coordinates, for final POSCAR printout          
   if (not cartesian):
      xyz_new=trans_cart2frac(xyz_new,natoms,a_vec,b_vec,c_vec)

   sys.stdout=original_stdout
   print(" done!\n")
    


# If the bottom option is activated, turn on the selective switch
if bottom:
   selective = True
   if ((not multiply_job) and (not shift_job) and (not frac2cart) and (not cart2frac)):
      xyz_new=xyz
      select_new=[]


# Finally, write out the modified POSCAR file POSCAR_mod

print(" New file POSCAR_mod will be written...")

original_stdout=sys.stdout
with open("POSCAR_mod","w") as f:
   sys.stdout = f

   # rescale the z-component of the unit cell vector such that the desired 
   # total length in z-direction of the system is retained

   print("Modified POSCAR written by modify_poscar.py")
   print(sys_scale)
   print("{:18.11f}".format(a_vec[0]) + "  " + "{:18.11f}".format(a_vec[1]) +
                  "  " + "{:18.11f}".format(a_vec[2]))
   print("{:18.11f}".format(b_vec[0]) + "  " + "{:18.11f}".format(b_vec[1]) + 
                  "  " + "{:18.11f}".format(b_vec[2]))
   print("{:18.11f}".format(c_vec[0]) + "  " + "{:18.11f}".format(c_vec[1]) + 
                  "  " + "{:18.11f}".format(c_vec[2]))
   for i in range(nelem):
      f.write(elements[i] + " ")
   print(" ")
# Combine element sections of both molecules  
   elem_num_all=[]
   elements_all =[]
   for i in range(nelem):
      elem_num_all.append(elem_num[i])
   for elem in elem_num_all:
      f.write(str(elem) + " ")
   
   print(" ")
# Write combined element numbers section
   if selective == 1:
      print("Selective Dynamics")
   if cartesian:
      print("Cartesian")
   else:
      print("Direct")

   if selective:
      for i in range(natoms):
# If the bottom option is activated, freeze the atoms below a certain z-value          
         if bottom:
            if xyz_new[i][2] <= freeze_bottom:
               select_new.append(" F F F ")
            else:   
               select_new.append(" T T T ")
# For direct coordinates, move all atoms into the central unit cell!               
         if (not cartesian):
            for j in range(3):
               while xyz_new[i][j] < 0.0:
                  xyz_new[i][j] = xyz_new[i][j] + 1.0
               while xyz_new[i][j] > 1.0:   
                  xyz_new[i][j] = xyz_new[i][j] - 1.0 

         print("{:20.11f}".format(xyz_new[i][0]) + " " + "{:20.11f}".format(xyz_new[i][1]) + " " +
              "{:20.11f}".format(xyz_new[i][2]) + "     " + select_new[i])
   else:       
      for i in range(natoms):
         if (not cartesian):
            for j in range(3):
               while xyz_new[i][j] < 0.0:
                  xyz_new[i][j] = xyz_new[i][j] + 1.0 
               while xyz_new[i][j] > 1.0:
                  xyz_new[i][j] = xyz_new[i][j] - 1.0          
         print("{:20.11f}".format(xyz_new[i][0]) + " " + "{:20.11f}".format(xyz_new[i][1]) + " " +
              "{:20.11f}".format(xyz_new[i][2])) 

sys.stdout=original_stdout
print(" done!")

print ('''
 Execution of modify_poscar.py successfully finished!
 File POSCAR_mod has been written!
   ''')
