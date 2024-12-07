#!/usr/bin/env python3
#   
#    analyze_poscar: Analyze the structure of a system in a 
#      POSCAR file. Different options are possible: measurement
#      of inclination angles of molecular groups towards coordinate
#      angles (e.g., for the rotation of benzene rings)
#    Part of VASP4CLINT
#     Julien Steffen, 2024 (julien.steffen@fau.de)
#

import sys
import os
import re
import numpy as np
from numpy import linalg as LA

print('''
 This script takes a POSCAR file in direct or cartesian coordinates 
 and analyzes/measures its structure, based on given keywords. 
 The list of possible options:
  -density : Calculates the density of the POSCAR file (g/cm^3)
  -atoms=a,b,c : Selects a list of atoms to be used for an evaluation
     option, for example to define a plane for angle measurement.
     Example: atoms=3,6,17,18,35
  -plane_angle=[coord-plane]: Activates the calculation of the angle
     of the plane defined by the given atoms and a coordinate plane
     given by the [coord-plane], which can be xy, xz or yz.
     Image flags are corrected automatically.
     Example: plane_angle=xy
 More than one job can be done at once, the ordering of operation is the 
 same as the ordering of keywords above 
''')

atoms_defined=False
angle_job=False
dens_job=False

# Read in the command line arguments
for arg in sys.argv:
   if arg == "-density":
      dens_job=True
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-atoms":
         atom_list=actval.split(",")
         atoms_defined=True
      if param == "-plane_angle":
         coord_plane=actval
         angle_job=True

if dens_job:
   print("The density of the unit cell in POSCAR will be calculated.")
else:  

   if ((not angle_job)):
      print("Please give a least one valid evaluation job!")
      sys.exit(1)

   if ((not atoms_defined)):
      print("Please give a number of atoms with -atoms keyword!")
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

# A dictionary for the atomic masses
elements_dict = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}

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


#
#  This function fits a 3D plane to a number of given points (atom
#  coordinates in 3D space. 
#
def fit_plane(points):
    # Convert the list of points to a NumPy array
    # Extract the x, y, z coordinates from the points
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    # Create the design matrix A
    A = np.c_[x, y, np.ones(points.shape[0])]

    # Solve for the plane coefficients using least squares
    # The plane equation is Ax + By + Cz + D = 0, with the normal vector being (A, B, C)
    coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)

    # The normal vector is (A, B, -1) because we solve for Ax + By + z = D
    normal_vector = np.array([coeffs[0], coeffs[1], -1])

    # Normalize the normal vector
    normal_vector /= np.linalg.norm(normal_vector)

    return normal_vector

#  Simple function to determine the angle between two 3D vectors (used for 
#   calculation of angle between atom-defined plane and chosen coordinate 
#   plane 

def vector_angle(v1, v2):
    # Normalize the vectors
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)

    # Compute the dot product
    dot_product = np.dot(v1_norm, v2_norm)

    # Ensure the dot product to be in the range [-1, 1] to avoid numerical issues
    dot_product = np.clip(dot_product, -1.0, 1.0)

    # Compute the angle in radians
    angle = np.arccos(dot_product)

    # Convert the angle to degrees
    angle = angle*180/(np.pi)

    if angle < 0:
       angle = -angle

    # Always take an angle between 0 and 90 degrees
    if angle >= 90:
       angle = 180 - angle

    angle = float(angle[0])
 
    return angle 

#
#  Calculate the density of the unit cell in the POSCAR file
#
if dens_job:
#
#    First, calculate the volume of the unit cell
#  
   a_np=np.array(a_vec)
   b_np=np.array(b_vec)
   c_np=np.array(c_vec) 

   cross_product = np.cross(b_np, c_np)

   triple_product = np.dot(a_np, cross_product)

   volume=abs(triple_product)

   print("vol",volume)

   mass_tot=0.0
   for i in range(nelem):
      el_act=elements[i]
      el_act=elements[i].upper()
      mass_act=elements_dict[elements[i]]
      mass_tot=mass_tot+mass_act*elem_num[i]

   density=mass_tot/volume*1.66053906660
   print(density)

#  Obtain the other set of structures, respectively (fractional if 
#   cartesian is given, cartesian if fractional is given)

if not cartesian:
   xyz_frac=xyz
   xyz=trans_frac2cart(xyz_frac,natoms,a_vec,b_vec,c_vec)

if cartesian:
   xyz_frac=trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec)


if angle_job:

   #  Determine, which coordinate plane (xy, xz or yz) shall be taken
   #  as reference for coordinate measurement
   ref_vector=np.zeros(3)
   if coord_plane == "xy":
      ref_vector=[(0,0,1)]
  

   #  Determine the list of atom coordinates of the atoms to be fitted   
   points=np.zeros((len(atom_list),3))

   if (len(atom_list) < 3):
       print("Please give at least three atoms to define a plane!")
       sys.exit(1)

   #  Correct for unit cell problems: move all atoms together if their 
   #  coordinate (x, y or z) differ by more than 0.5 in direct 
   #  coordinates (correct for image flags)

   pos_ref=np.zeros(3)
   pos_ref[0]=xyz_frac[int(atom_list[0])-1][0]
   pos_ref[1]=xyz_frac[int(atom_list[0])-1][1]
   pos_ref[2]=xyz_frac[int(atom_list[0])-1][2]

   for index in atom_list:
      for i in range(3):
         while ((xyz_frac[int(index)][i] - pos_ref[i]) > 0.5):
             xyz_frac[int(index)][i] = xyz_frac[int(index)][i] - 1.0
         while ((xyz_frac[int(index)][i] - pos_ref[i]) < -0.5):
             xyz_frac[int(index)][i] = xyz_frac[int(index)][i] + 1.0 

   xyz=trans_frac2cart(xyz_frac,natoms,a_vec,b_vec,c_vec)
 
   #  Store the coordinates of the points for the transmission to the 
   #  angle calculation function

   i=0
   for index in atom_list:
      index=int(index)
      points[i][0]=xyz[index-1][0]
      points[i][1]=xyz[index-1][1]
      points[i][2]=xyz[index-1][2]
      i=i+1

   normal_vector = fit_plane(points)

   #  Now calculate the angle between the normal vector of the fitted 
   #  plane and the normal vector of the chosen coordinate plane

   angle = vector_angle(ref_vector,normal_vector)   

   print(" The plane spanned up by the chosen points and the " + str(coord_plane) + "-plane") 
   print("   have an angle of " + str(angle) + "°")


print ('''
 Execution of analyze_poscar.py successfully finished!
   ''')
