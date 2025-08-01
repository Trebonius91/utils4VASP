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
 This script takes a POSCAR file in direct or cartesian coordinates 
 and analyzes/measures its structure, based on given keywords. 
 The list of possible options:
  -overview : Print overview of utils4VASP scripts/programs
  -gen_kpoints : Generates a KPOINTS file for a given sampling density
     If vacuum of more than 7 Angs. is introduced in one axis, the 
     number of k-points in this direction will be set to 1.
  -k_density=[value] : The desired reciprocal space sampling density
     in 2pi/Angstrom (default: 0.03)
  -gen_potcar : Generates a POTCAR file for the given POSCAR
  -potcar_path=[path] : The path with the VASP POTCAR collection for all
     elements (default: ~/ (home-directory)) 
  -count_elec : Determine the total number of valence electrons of the 
     system (NELECT keyword), by assuming that all elements are modeled
     by their standard PAW PBE potentials.
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
kpoints_job=False
potcar_job=False
elec_job=False
k_dens=0.04
pot_path="~/"  # The default POTCAR path, change if needed!
vac_min=7.0

# Read in the command line arguments
for arg in sys.argv:
   if arg == "-density":
      dens_job=True
   if arg == "-gen_kpoints":
      kpoints_job=True   
   if arg == "-gen_potcar":
      potcar_job=True
   if arg == "-count_elec":
      elec_job=True
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-atoms":
         atom_list=actval.split(",")
         atoms_defined=True
      if param == "-plane_angle":
         coord_plane=actval
         angle_job=True
      if param == "-k_density":
         k_dens=actval
      if param == "-potcar_path":
         pot_path=actval

if dens_job:
   print(" The density of the unit cell in POSCAR will be calculated.")
elif angle_job: 
   print(" The angle between a plane and a coord. plane will be calculated.")
elif kpoints_job:
   print(" A KPOINTS file for the POSCAR file will be written.")
elif potcar_job: 
   print(" A POTCAR file for the POSCAR file will be written.")
elif elec_job: 
   print(" The total number of valence electrons will be calculated.")
else:
   print(" Please give a valid job!\n")
   sys.exit(1)

if angle_job:
   if ((not atoms_defined)):
      print(" Please give a number of atoms with -atoms keyword!")
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

# A dictionary for the number of valence electrons (standard PBE PAW)
element_valence_dict = {"H": 1.0, "He": 2.0, "Li": 1.0, "BE": 2.0, "B": 3.0,\
                 "C": 4.0, "N": 5.0, "O": 6.0, "F": 7.0, "NE": 8.0, "NA": 1.0,\
                 "MG": 2.0, "AL": 3.0, "SI": 4.0, "P": 5.0, "S": 6.0, "CL": 7.0,\
                 "AR": 8.0, "SC": 3.0, "TI": 4.0, "V": 5.0, "CR": 6.0,\
                 "MN": 7.0, "FE": 8.0, "CO": 9.0, "NI": 10.0, "CU": 11.0,\
                 "ZN": 12.0, "GA": 3.0, "GE": 4.0, "AS": 5.0, "SE": 6.0,\
                 "BR": 7.0, "KR": 8.0, "MO": 6.0, "TC": 7.0, "RU": 8.0,\
                 "RH": 9.0, "PD": 10.0, "AG": 11.0, "CD": 12.0, "IN": 3.0,\
                 "SN": 4.0, "SB": 5.0, "TE": 6.0, "I": 7.0, "XE": 8.0,\
                 "LA": 11.0, "CE": 12.0, "PR": 13.0, "ND": 14.0, "PM": 16.0,\
                 "SM": 16.0, "EU": 17.0, "GD": 18.0, "TB": 19.0, "DY": 20.0,\
                 "HO": 21.0, "ER": 22.0, "TM": 23.0, "YB": 24.0, "LU": 25.0,\
                 "HF": 4.0, "TA": 5.0, "W": 6.0, "RE": 7.0, "OS": 8.0, "IR": 9.0,\
                 "PT": 10.0, "AU": 11.0, "HG": 12.0, "TL": 3.0, "PB": 4.0,\
                 "BI": 5.0, "PO": 6.0, "AT": 7.0, "RN": 8.0, "AC": 11.0,\
                 "TH": 12.0, "PA": 13.0, "U": 14.0, "NP": 15.0, "PU": 16.0,\
                 "AM": 17.0, "CM": 18.0}
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
#  Generate a POTCAR file for the current POSCAR file, with a given path
#
if potcar_job:
   if os.path.isfile("POTCAR"):
      os.system("rm POTCAR")
   for el in elements:
      os.system("cat "+pot_path+"potcar/PAW_PBE.52/"+el+"/POTCAR >> POTCAR")

   print(" The POTCAR file was written.")
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
      mass_act=elements_dict[el_act]
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

#
#  Generate a KPOINTS file for the current POSCAR file, with a given sampling
#   density in the Brillouin zone
#
if kpoints_job:
#
#    First, calculate the volume of the unit cell
#  
   a_np=np.array(a_vec)
   b_np=np.array(b_vec)
   c_np=np.array(c_vec)

   cross_product = np.cross(b_np, c_np)

   triple_product = np.dot(a_np, cross_product)

   volume=abs(triple_product)


   kpoint_a_ex=np.linalg.norm(np.cross(b_np,c_np))/(volume*k_dens)
   kpoint_a = round(kpoint_a_ex)
   if kpoint_a < 1:
      kpoint_a = 1
   print(" k-points along a-axis:",kpoint_a_ex," (rounded:",kpoint_a,")")
   kpoint_b_ex=np.linalg.norm(np.cross(a_np,c_np))/(volume*k_dens)
   kpoint_b = round(kpoint_b_ex)
   if kpoint_b < 1:
      kpoint_b = 1
   print(" k-points along b-axis:",kpoint_b_ex," (rounded:",kpoint_b,")")
   kpoint_c_ex=np.linalg.norm(np.cross(a_np,b_np))/(volume*k_dens)
   kpoint_c = round(kpoint_c_ex)
   if kpoint_c < 1:
      kpoint_c = 1
   print(" k-points along c-axis:",kpoint_c_ex," (rounded:",kpoint_c,")")
   kpoint_c_ex=np.linalg.norm(np.cross(a_np,b_np))/(volume*k_dens)

#
#    Check if vacuum exists in any of the three directions in space
#    If yes, set the number of k-points to 1
#    
   for i in range(natoms):
      for j in range(3):
         while xyz_frac[i][j] > 1.0:
            xyz_frac[i][j]=xyz_frac[i][j]-1.0
         while xyz_frac[i][j] < 0.0:
            xyz_frac[i][j]=xyz_frac[i][j]+1.0


   vacuum=np.zeros((3))
   for i in range(3):
      xyz_axis=np.zeros((natoms))
      for j in range(natoms):
         xyz_axis[j]=xyz_frac[j][i]
      xyz_axis=sorted(xyz_axis)
      distances=np.zeros((natoms))
      distances[0]=xyz_axis[0]-xyz_axis[natoms-1]+1.0
      for j in range(natoms-1):
         distances[j+1]=xyz_axis[j+1]-xyz_axis[j]
      if i == 0:
         vacuum[i]=distances.max()*np.linalg.norm(a_np)    
      if i == 1:
         vacuum[i]=distances.max()*np.linalg.norm(b_np)
      if i == 2:
         vacuum[i]=distances.max()*np.linalg.norm(c_np)
   
   if (vacuum[0] > vac_min):
      print(" ")
      print(" A vacuum of more than 7 Ang. (",vacuum[0],") exists in a direction!")
      print(" The number of k-points along a will be set to 1.")
      k_point_a=1
   if (vacuum[1] > vac_min):
      print(" ")
      print(" A vacuum of more than 7 Ang. (",vacuum[1],") exists in b direction!")
      print(" The number of k-points along b will be set to 1.")
      k_point_b=1
   if (vacuum[2] > vac_min):
      print(" ")
      print(" A vacuum of more than 7 Ang. (",vacuum[2],") exists in c direction!")
      print(" The number of k-points along c will be set to 1.")
      k_point_c=1



   original_stdout=sys.stdout
   with open("KPOINTS","w") as f:
      sys.stdout = f
      print("k-points, density = ",k_dens,"2pi/Ang.")
      print("0")
      print("Gamma")
      print(kpoint_a,"  ",kpoint_b,"  ",kpoint_c)
      print("0 0 0 ")
   sys.stdout = original_stdout
   print(" ")
   print(" The KPOINTS file was written.")


#   Determine the total number of electrons in the system

if elec_job:
   elec_tot=0.0
   print(" Analysis of the valence electrons... :")
   for i in range(nelem):
      el_act_print=elements[i]
      el_act=elements[i].upper()
      elec_act=element_valence_dict[el_act]
      elec_tot=elec_tot+elec_act*elem_num[i]
      print("  *",el_act_print,":  ",elec_act," valence electrons  (",elem_num[i]," atoms)")
   print(" ")
   print(" Total number of valence electrons: ",elec_tot)   
   original_stdout=sys.stdout
   with open("NELECT.dat","w") as f:
      sys.stdout = f
      print(elec_tot)
   sys.stdout = original_stdout

   print(" File NELECT.dat with valence electron number written.")   

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
