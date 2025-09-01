#!/usr/bin/env python3
#   
#    manage_mlff_md: Performs a MD trajectory calcuation with a 
#      potentially unstable VASP ML-FF, with the current focus 
#      on slab calculations. Whenever an atom has moved too far 
#      from the initial slab geometry away (z-coordinate), the
#      current calculation is canceled and a previous structure 
#      from the XDATCAR is taken as new POSCAR
#    Part of VASP4CLINT
#     Julien Steffen, 2023 (julien.steffen@fau.de)
#

import sys
import os
import re
import time
import math
import numpy as np

print('''
 This script performs a MD trajectory calcuation with a
 potentially unstable VASP ML-FF, with the current focus
 on slab calculations. Whenever an atom has moved too far
 from the initial slab geometry away (z-coordinate), the
 current calculation is canceled and a previous still resonable
 structure from the XDATCAR is taken as new POSCAR.
 The atoms in the POSCAR file should have a slab structure 
 and at z=0 should be vacuum!
 During the calculations, the finished part of the trajectory 
 are collected together in a XDATCAR_all file.
 The current status of the script is written to file 
 'manage_mlff_md.log'.
 If the calculation shall be terminated, touch a file named 
 'kill' in the current folder.
 To start the script, give the following command:
 manage_mlff_md.py -tolerance=[value] -buffer=[value]
   -tolerance=[value] gives the maximum distance in Angstrom
     an atom can depart from the surface before the calculation
     is restarted (resonable choice: 2-3 Angstroms(?))
   -buffer=[value] gives for the reset a distance any atom 
     should below the maximum tolerance such that no direct
     restart for the next trajectory is needed 
     (resonable choice: 0.2-0.5 Angstroms(?))
   -check_range=[number1],[number2]: Checks only a part of the 
     atoms in the system and resets them if too far away. Can 
     be useful if, e.g., a surface slab is simulated together
     with a gas phase. Pick the atom numbers of the slab into 
     this range if no atom of the slab shall go into the gas 
     phase durinthy dynamics.
''')

#
#     Maximum elongation of atoms in both sides of the 
#     slab, in Angstroms
#
z_tolerance=-100.0
#
#     Tolerance buffer for newly written POSCAR, such that it 
#     is not directly canceled
#
tol_buffer=-100.0
#
#     First and last atom index to be used
#
check_range=False
ind_first=0
ind_last=0
#
#     Read in needed parameters from command line 
#
for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-tolerance":
         z_tolerance=float(actval)
      if param == "-buffer":
         tol_buffer=float(actval)
      if param == "-check_range":
         actval1,actval2=actval.split(",")
         ind_first=int(actval1)
         ind_last=int(actval2)
         check_range = True

if z_tolerance < -10.0:
   print("Please give a value with the -tolerance command!")
   exit(1)
if tol_buffer < -10.0:  
   print("Please give a value with the -buffer command!")
   exit(1)    

#
#     From now on, print everything into the logfile!
#
            

with open("manage_mlff_md.log", "w") as logfile:
   get_start=False

#
#    Open POSCAR file to find allowed minimum and maximum z-values 
#
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
         print(" The POSCAR file has selective dynamics.", file=logfile)
         print(" ", file=logfile)
         logfile.flush()
      if selective:
         line = infile.readline().rstrip("\n")
      line_parts = line.split()
      if (line_parts[0] != "Direct" and line_parts[0] != "direct"
                and line_parts[0] != " Direct" and line_parts[0] != " direct"):
         cartesian=True
         print(" The POSCAR has cartesian coordinates.", file=logfile)
         logfile.flush()
      else:
         print(" The POSCAR has direct coordinates.", file=logfile)
         logfile.flush()

      print(" ")
      for i in range(natoms):
         line = infile.readline().rstrip("\n")
         xyz_read = line.rstrip().split()[0:3]
         if selective:
            select_read = line.rstrip().split()[3:6]
            coord_select.append(select_read[0] + " " + select_read[1] + " " + select_read[2])
         for j in range(3):
            xyz[i][j]=float(xyz_read[j])

   poscar_in.close()
#
#    Set checked range to all atoms if not defined otherwise
#    Else, set it to python range (atom1 has index 0)
#
   if not check_range:
      ind_first=0
      int_last=natoms-1
   else:
      ind_first=ind_first-1
      ind_last=ind_last-1

#
#    Determine minimum and maximum values 
#

   if cartesian:
      zmax=abs((np.amax(xyz[ind_first:ind_last][:],axis=0))[2])
      zmin=abs((np.amin(xyz[ind_first:ind_last][:],axis=0))[2])
   else:
      zmax=abs((np.amax(xyz[ind_first:ind_last][:],axis=0))[2]*c_vec[2])
      zmin=abs((np.amin(xyz[ind_first:ind_last][:],axis=0))[2]*c_vec[2])

#
#    Read in INCAR file for number of MD steps 
#

   incar_name="INCAR"

   incar_act = open(incar_name,"r")

   with incar_act as infile:
      incar_array=infile.readlines()
   incar_act.close()

   step_freq = 1  # defaut value for output frequency
   for line in incar_array:
      if re.search("NSW",line):
         split_act=line.split()
         steps_all=int(split_act[2])
      if re.search("ML_OUTBLOCK",line):
         split_act=line.split()
         step_freq=int(split_act[2])
 
   xdat_steps=int(steps_all/step_freq)

#
#     Remove old CONTCAR, XDATCAR and XDATCAR_old files
#
   if os.path.isfile("CONTCAR"):
      os.system("rm CONTCAR")
   if os.path.isfile("XDATCAR"):
      os.system("rm XDATCAR")
   if os.path.isfile("XDATCAR_all"):
      os.system("rm XDATCAR_all")

#
#     Write header of new unified XDATCAR file
#
   original_stdout=sys.stdout
   with open("XDATCAR_all","w") as f:
      sys.stdout = f
      print("Combined XDATCAR written by manage_mlff_md.py")
      print(sys_scale)
      print(str(a_vec[0]) + "  " + str(a_vec[1]) + "  " + str(a_vec[2]))
      print(str(b_vec[0]) + "  " + str(b_vec[1]) + "  " + str(b_vec[2]))
      print(str(c_vec[0]) + "  " + str(c_vec[1]) + "  " + str(c_vec[2]))
      for i in range(nelem):
         f.write(elements[i] + " ")
      print(" ")
      elem_num_all=[]
      elements_all =[]
      for i in range(nelem):
         elem_num_all.append(elem_num[i])
      for elem in elem_num_all:
         f.write(str(elem) + " ")
      print(" ")
      if selective == 1:
         print("Selective Dynamics")
   sys.stdout=original_stdout


#
#    Save initial INCAR and POSCAR files 
#

   os.system("cp INCAR INCAR_initial")
   os.system("cp POSCAR POSCAR_initial")
#
#    Start first calculation
#
   os.system("sbatch slurm_script > slurm_log")

   slurm_logfile=open("slurm_log","r")
   with slurm_logfile as infile:
      line = infile.readline()
      line_split = line.rstrip().split()
   
   slurm_logfile.close()   
   job_num=line_split[3]


   print("First job (No.",job_num,") started!", file=logfile)
   logfile.flush()
   time.sleep(3)

#
#     Go in big loop, check status of CONTCAR file every 30 seconds
#
   unfinished=True
   reset=False
   finished_steps=0
   while unfinished:    
      time.sleep(10)
#
#     Check if kill file is present, if yes, 
#
      if os.path.isfile("kill"):
         print("The file kill has been found!", file=logfile)
         print("The calculations will be aborted...", file=logfile)
         logfile.flush()
         os.system("scancel " + str(job_num))
         sys.exit(0)
#
#     Check status of job
#
      os.system("squeue --job " + str(job_num) + " > slurm_log")
      slurm_logfile=open("slurm_log","r")
      with slurm_logfile as infile:
         line = infile.readline()
         line = infile.readline()
         line_split = line.rstrip().split()
      try:   
         job_status=line_split[4]
      except Exception as inst:
         print("Job has been canceled!", file=logfile)
         logfile.flush()
         try:
            xdat_name="XDATCAR"
            xdat_act = open(xdat_name,"r")

            with xdat_act as infile:
               xdat_array=infile.readlines()
#
#     Determine number of XDATCAR frames 
#
            nlines=len(xdat_array)
            nframes=math.floor((nlines-7)/(natoms+1))
         except Exception as inst2:
            nframes=0 
#
#     If the total number of frames has been reached, exit 
#
         if nframes + finished_steps >= int(steps_all/step_freq):
            logfile.flush()
            xdat_all=open("XDATCAR_all","a")
            with xdat_all as infile:
                for i in range (7,7+(nframes)*(natoms+1)):
                   infile.write(xdat_array[i])
            xdat_all.close()
            print("HURRAY! The ML-FF trajectory is finished!", file=logfile)
            os.system("mv XDATCAR_all XDATCAR")
            sys.exit(0)

         else:
            print("The calculation is not finished yet..", file=logfile)
            print("It will be restarted.", file=logfile)
            logfile.flush()
            os.system("sbatch slurm_script > slurm_log")

            slurm_logfile=open("slurm_log","r")
            with slurm_logfile as infile:
               line = infile.readline()
               line_split = line.rstrip().split()

            slurm_logfile.close()
            job_num=line_split[3]

      if job_status == "PD":
         print("Job still waits...", file=logfile)
         logfile.flush()
#
#     If no job is active at all, check if all steps already were calculated
#     If yes, exit the script with success
#     If no, simply restart the job (possibly a VASP error)
#
      elif job_status == "job":
         print("Job has been canceled!", file=logfile)
         logfile.flush()

         xdat_name="XDATCAR"
         xdat_act = open(xdat_name,"r")

         with xdat_act as infile:
            xdat_array=infile.readlines()
#
#     Determine number of XDATCAR frames 
#
         nlines=len(xdat_array)
         nframes=math.floor((nlines-7)/(natoms+1))
#
#     If the total number of frames has been reached, exit 
#
         if nframes + finished_steps >= int(steps_all/step_freq):

            xdat_all=open("XDATCAR_all","a")
            with xdat_all as infile:
                for i in range (7,7+(nframes)*(natoms+1)):
                  infile.write(xdat_array[i])
            xdat_all.close()

            print("HURRAY! The ML-FF trajectory is finished!", file=logfile)
            logfile.flush()
            os.system("mv XDATCAR_all XDATCAR")
            sys.exit(0)
         else:    
            print("The calculation is not finished yet..", file=logfile)
            print("It will be restarted.", file=logfile)
            logfile.flush()
            os.system("sbatch slurm_script > slurm_log")

            slurm_logfile=open("slurm_log","r")
            with slurm_logfile as infile:
               line = infile.readline()
               line_split = line.rstrip().split()

            slurm_logfile.close()
            job_num=line_split[3]

      elif job_status == "R":
         print("Job is currently active", file=logfile)
         logfile.flush()
      else:
         print("Job status is unknown! (",job_status,")", file=logfile) 
         logfile.flush()
      if job_status == "R":
         loop_exit=False 
         if os.path.isfile("CONTCAR"):
            contcar_name="CONTCAR"

            contcar_act = open(contcar_name,"r")

            with contcar_act as infile:
               line = infile.readline()
               line = infile.readline()
               line = infile.readline()
               line = infile.readline()
               line = infile.readline()
               line = infile.readline()
               line = infile.readline()
               line = infile.readline()
               if selective:
                  line = infile.readline()
 
               for i in range(natoms):
                  try:
                     line = infile.readline().rstrip("\n")
                     xyz_read = line.rstrip().split()[0:3]
                  except Exception as inst:
                      loop_exit=True 
                  if selective:
                     select_read = line.rstrip().split()[3:6]
                     coord_select.append(select_read[0] + " " + select_read[1] + " " + select_read[2])
                  for j in range(3):
                     try: 
                        xyz[i][j]=float(xyz_read[j])
                     except Exception as inst:
                        loop_exit=True
                        break
                  if loop_exit:
                     break
            if loop_exit:
               print("CONTCAR file was not complete!", file=logfile)
               logfile.flush()
               zmax_act=zmax
               zmin_act=zmin
#
#     Determine current z minimum and maximum values 
#
            else:
               zmax_act=abs((np.amax(xyz[ind_first:ind_last][:],axis=0))[2])*c_vec[2]
               zmin_act=abs((np.amin(xyz[ind_first:ind_last][:],axis=0))[2])*c_vec[2]
#
#     Determine z-elongations and check if they are too large
#
               z_min_elong=zmin-zmin_act
               z_max_elong=zmax_act-zmax
               print("The current CONTCAR appears to be OK.", file=logfile)
               logfile.flush()
               if z_min_elong > z_tolerance:
                  print("One atom is too far below the slab!", file=logfile) 
                  logfile.flush()
                  reset=True
               if z_max_elong > z_tolerance:
                  print("One atom is too far above the slab!", file=logfile)
                  logfile.flush()
                  reset=True
#
#     Reset calculation if z-elongation is too large
#
               if reset:  
#
#     Cancel old calculation
#
                  os.system("scancel " + str(job_num))
#
#     Remove actual CONTCAR file 
#
                  os.system("rm CONTCAR")
#
#     Check actual XDATCAR file
#

                  print("The calcuation will be reset!", file=logfile)
                  logfile.flush()
                  xdat_name="XDATCAR"

                  xdat_act = open(xdat_name,"r")

                  with xdat_act as infile:
                     xdat_array=infile.readlines()
#
#     Determine number of XDATCAR frames 
#
                  nlines=len(xdat_array)
                  nframes=math.floor((nlines-7)/(natoms+1))
                  print("number of frames:",str(nframes), file=logfile)
                  logfile.flush()
#
#     Search for first frame that is below the threshold..
#
                  frame_act=nframes-100
                  if frame_act <= 0:
                     frame_act = 1
                  while reset:
                     xyz_act=np.zeros((natoms,3))
                     for i in range(natoms):
                        line = xdat_array[7+(frame_act-1)*(natoms+1)+i+1]
                        xyz_read = line.rstrip().split()[0:3]
                        for j in range(3):
                           xyz_act[i][j]=float(xyz_read[j])


                     zmax_check=abs((np.amax(xyz_act[ind_first:ind_last][:],axis=0))[2])*c_vec[2]
                     zmin_check=abs((np.amin(xyz_act[ind_first:ind_last][:],axis=0))[2])*c_vec[2]
                     z_min_elong=zmin-zmin_check
                     z_max_elong=zmax_check-zmax
                     if (z_min_elong > (z_tolerance-tol_buffer)) or (z_max_elong > (z_tolerance-tol_buffer)):
                        frame_act=frame_act-5
#
#     In order to avoid an endless loop: reset the calculation even if the threshold is not 
#      met if no calculation in the trajectory has been found that obeys the criteria
#
                        if frame_act <= 0:
                           frame_act = 1
                           reset = False
                     else:
#
#     If a good frame was found, write it to a POSCAR file
#
                        reset=False
                        original_stdout=sys.stdout
                        with open("POSCAR","w") as f:
                           sys.stdout = f
  
# rescale the z-component of the unit cell vector such that the desired 
# total length in z-direction of the system is retained

                           print("Modified POSCAR written by manage_mlff_md.py")
                           print(sys_scale)
                           print(str(a_vec[0]) + "  " + str(a_vec[1]) + "  " + str(a_vec[2]))
                           print(str(b_vec[0]) + "  " + str(b_vec[1]) + "  " + str(b_vec[2]))
                           print(str(c_vec[0]) + "  " + str(c_vec[1]) + "  " + str(c_vec[2]))
                           for i in range(nelem):
                              f.write(elements[i] + " ")
                           print(" ")
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
                           print("Direct")
                           if selective:
                              for i in range(natoms):
                                  print(str(xyz_act[i][0]) + " " + str(xyz_act[i][1]) + " " +
                                       str(xyz_act[i][2]) + " " + select_new[i])
                           else:
                              for i in range(natoms):
                                  print(str(xyz_act[i][0]) + " " + str(xyz_act[i][1]) + " " +
                                       str(xyz_act[i][2]))
                        sys.stdout=original_stdout
#
#   Update total number of frames
#
                        xdat_steps=xdat_steps-frame_act
                        finished_steps=finished_steps+frame_act
#
#   Write new INCAR file 
#
                        original_stdout=sys.stdout
                        with open("INCAR","w") as f:
                           sys.stdout = f
                           for line in incar_array:
                              if re.search("NSW",line):
                                 print("NSW = ", str(xdat_steps*step_freq))
                              else: 
                                 print(line.rstrip())
                        sys.stdout=original_stdout           
#
#   Write all useful frames to XDATCAR_all
#
 
                        print("First ",frame_act,"frames of actual XDATCAR are appended to XDATCAR_all", file=logfile)
                        logfile.flush()
                        xdat_all=open("XDATCAR_all","a")
                        with xdat_all as infile:
                            for i in range (7,7+(frame_act-1)*(natoms+1)):
                              infile.write(xdat_array[i]) 
                        xdat_all.close()
                        print("Already ",str(finished_steps*step_freq)," of ",str(steps_all)," MD steps finished.", 
                                  file=logfile)
                        logfile.flush()
#
#   Start new calculation
#
                        os.system("sbatch slurm_script > slurm_log")

                        slurm_logfile=open("slurm_log","r")
                        with slurm_logfile as infile:
                           line = infile.readline()
                           line_split = line.rstrip().split()

                        slurm_logfile.close()
                        job_num=line_split[3]
                        print("New job (No.",job_num,") started!", file=logfile)
                        logfile.flush()

         else: 
            print("Calculation still in ML-FF initialization mode.", file=logfile) 
            logfile.flush()

