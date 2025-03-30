from ase import units
from ase import build
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.md.npt import NPT
import numpy as np
import time
from timeit import default_timer as timer
import datetime
from datetime import timedelta

from mace.calculators import mace_mp

#
#    Begin the time measurement
#
start = timer()
#
#   Define the MACE model (local model file)
#
macemp = mace_mp(model="mace.model") 

#
#   General input variables
#
#   The desired temperature in K
temperature = 300.0
#   The desired external pressure in bar
sigma = 1.0
#   The desired total number of MD steps
nsteps = 30000
#   The desired time step in fs 
deltat = 1.0
#   The write out frequency (every N steps, energy and structure are written)
nfreq = 100
#   The Parrinello-Rahman barostat time constant in GPa 
pfactor = 2e6
#   The Nose-Hoover thermostat time constant in fs
ttime = 20.0



#
#    Read in the geometry from the POSCAR file
#
atoms = read('POSCAR')

natoms = len(atoms)
#
#    Setup the MACE foundation model
#
atoms.calc = macemp

#
#    Initialize velocities.
#
MaxwellBoltzmannDistribution(atoms, temperature * units.kB)

#
#    Initialize the NVT dynamics object
#    Cell is only variable in x,y and z axis lengths 
#    but keeps rectangular shape
#
dyn = NPT(atoms,
      deltat*units.fs,
      temperature_K = temperature,
      externalstress = sigma*units.bar,
      ttime = ttime*units.fs,
      pfactor = pfactor*units.GPa*(units.fs**2),
      logfile = "dummy",
      trajectory = "traj_dum",
      loginterval=nfreq,
      mask=[[0,0,0],[0,0,0],[0,0,0]]
      )

#
#    The logfile for the current status and energy printout
#
logfile=open("mace.log","w")
gradfile=open("gradients.dat","w")
#
#    Open the XDATCAR file during the whole dynamics
#

with open("XDATCAR","w") as traj:

   def write_xyz():
      write(traj,atoms,format="vasp-xdatcar")

#  Print status and CONTCAR XDATCAR during MD routine
   def print_dyn():

#  Write new frame of XDATCAR
      write(traj,atoms,format="vasp-xdatcar")

#  Update the CONTCAR file       
      with open("CONTCAR","w") as contcar:
         write(contcar,atoms,format="vasp")

#  Write information and new frame of XDATCAR
      traj.flush()
      logfile.flush()
      gradfile.flush()
      step_act=dyn.get_number_of_steps()
      epot = atoms.get_potential_energy()
      ekin = atoms.get_kinetic_energy()
      temp = ekin / (1.5 * units.kB) / natoms
      etot = epot + ekin
      volume=atoms.get_volume()
      mass_tot=sum(atoms.get_masses())
      dens=mass_tot/volume*1.66053906660
      force_vec=atoms.get_forces()
      print('Step No.  %i  ' % (step_act),file=gradfile)
      for i in range(natoms):
         print('   %.8f    %.8f    %.8f  ' % (force_vec[i][0],force_vec[i][1],force_vec[i][2]),file=gradfile) 

      print('    %i      %.5f        %.5f       %.5f       %.5f       %.5f        %.5f ' % ((step_act),epot,ekin,etot,temp,volume,dens),file=logfile)


   dyn.attach(print_dyn, interval=nfreq)

#
#    The logfile for the current status and energy printout
#
   logfile=open("mace.log","w")
   gradfile=open("gradients.dat","w")

   print("Starting the MACE foundation model calculation for this system.",file=logfile)
   current_time = datetime.datetime.now()
   print("The current date and time is: ",current_time,file=logfile)
   print(" - Number of MD steps: ",nsteps,file=logfile)
   print(" - MD timestep (fs): ",deltat,file=logfile)
   print(" - Writing frequency: ",nfreq,file=logfile)
   print(" - Desired temperature (K): ",temperature,file=logfile)
   print(" - Parrinello-Rahman barostat parameter (GPa): ",pfactor,file=logfile)
   print(" - Nose-Hoover thermostat time constant (fs): ",ttime,file=logfile)
   print(" ",file=logfile)

#
#    Enter the big MD loop!
#

#
#    Write initial geometry to CONTCAR file
#
   
   with open("CONTCAR","w") as contcar:
      write(contcar,atoms,format="vasp")


   print("  Step     pot.energy (eV)    kin.energy (eV)   tot.energy (eV)   temperature (K)  volume(A^3)  density(g/cm^3)",file=logfile)
   epot = atoms.get_potential_energy()
   ekin = atoms.get_kinetic_energy()
   temp = ekin / (1.5 * units.kB) / natoms
   etot = epot + ekin
   volume=atoms.get_volume()
   mass_tot=sum(atoms.get_masses())
   force_vec=atoms.get_forces()
   dens=mass_tot/volume*1.66053906660
   print('Step No.  %i  ' % (0),file=gradfile)
   for i in range(natoms):
      print('   %.8f    %.8f    %.8f  ' % (force_vec[i][0],force_vec[i][1],force_vec[i][2]),file=gradfile)

   print('    %i       %.5f        %.5f       %.5f        %.5f        %.5f        %.5f ' % (0,epot,ekin,etot,temp,volume,dens),file=logfile)

   dyn.run(nsteps) 
       

#
#    End the time measurement
#
end = timer()
print("  ",file=logfile)
print("The calculation has been finished successfully!",file=logfile)
print("Total time spent (h:m:s): ", (timedelta(seconds=end-start)),file=logfile)
print("Exiting normally...",file=logfile)
