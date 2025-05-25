!
!    eval_bader: evaluates the result of a Bader charge calculation
!      after preprocessing the VASP calculation results with the 
!      chgsum.pl and bader programs by the Henkelman group.
!      The files POSCAR and ACF.dat need to be present
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!

program eval_bader 
implicit none 
integer::i,j,k
integer::readstat,el_num,idum,natoms
integer::readstat1,readstat2
character(len=120)::adum,arg
character(len=2),allocatable::el_init(:),el_names(:)
character(len=2),allocatable::el_atoms(:)
character(len=2),allocatable::el_mods(:)
character(len=5),allocatable::val_vec(:)
integer::mod_num
integer,allocatable::el_numbers(:)
real(kind=8),allocatable::charge_ref(:),charge_list(:)
real(kind=8),allocatable::charge_avg(:)  ! the final average charge per element
real(kind=8),allocatable::charge_mod(:)
real(kind=8)::charge_act,rdum,a_vec(3),b_vec(3),c_vec(3)
real(kind=8)::cell_scale
logical::coord_direct
real(kind=8),allocatable::xyz(:,:),xyz_print(:,:)


!
!    only print the overview of utils4VASP scripts/programs and stop
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-overview") then
      write(*,*)
      write(*,*) "utils4VASP: Setup and Evaluation of DFT and MLP simulations with VASP"
      write(*,*) "The following scripts and programs are contained:"
      write(*,*) "Setup:"
      write(*,*) " - gen_incar.py: Generate INCAR file templates for several job types"
      write(*,*) " - analyze_poscar.py: Analyze POSCAR, generate KPOINTS and POTCAR"
      write(*,*) " - build_alloy.py: Build regular alloy structures of various shapes"
      write(*,*) " - modify_poscar.py: Modify POSCAR: multiply, transform, shift, insert etc."
      write(*,*) " - cut_unitcell: Generate surface slab unit cells of arbitrary shape"
      write(*,*) " - build_adsorb.py: Place adsorbates on surfaces, set translation, rotation"
      write(*,*) " - split_freq: Divide frequency calculations for large molecules"
      write(*,*) "Evaluation:"
      write(*,*) " - modify_xdatcar: Modify trajectory files: multiply, shift, pick etc."
      write(*,*) " - analyze_bulk: Analyze bulk MD trajectories for RDFs, diffusion etc."
      write(*,*) " - analyze_slab: Analyze slab MD trajectories for RDFs, density etc."
      write(*,*) " - check_geoopt.py: Monitor geometry optimizations with selective dynamics"
      write(*,*) " - manage_cls: Prepare, evaluate core level energy calculations"
      write(*,*) " - eval_bader: Evaluate and visualize Bader charge calculations"
      write(*,*) " - eval_stm: Generate STM images with different settings"
      write(*,*) " - partial_dos: Plot atom and orbital resolved density of states"
      write(*,*) " - manage_neb.py: Setup, monitor and restart NEB calculations"
      write(*,*) "ML-FF:"
      write(*,*) " - mlff_select: Heuristic selection of local reference configurations"
      write(*,*) " - eval_vasp_ml.py: Visualize results of VASP ML-FF on the fly learnings"
      write(*,*) " - vasp2trainset: Generate ML-FF training sets from VASP calculations"
      write(*,*) " - mlp_quality: Determine quality of MLPs for VASP validation set"
      write(*,*) "Management:"
      write(*,*) " - md_long.sh: Automated restart of MD calculations with slurm"
      write(*,*) " - opt_long.sh: Automated restart of geometry optimizations with slurm"
      write(*,*) " - ml_long.sh: Automated restart of VASP ML-FF on the fly learnings"
      write(*,*) " - mace_long.sh: Automated restart of MACE MD trajectories with ASE "
      write(*,*)
      stop
   end if
end do

!
!    Print general information and all possible keywords of the program    
!
write(*,*)
write(*,*) "PROGRAM eval_bader: evaluation of Bader charge calculations"
write(*,*) " from preprocessed VASP output for arbitrary systems."
write(*,*) "Before starting this program, the files AECCAR0, AECCAR2 and "
write(*,*) " CHGCAR must be present as output and preprocessed with the "
write(*,*) " tools by the Henkelman group, findable at:"
write(*,*) " http://theory.cm.utexas.edu/henkelman/code/bader/"
write(*,*) "Step 1: chgsum.pl AECCAR0 AECCAR2"
write(*,*) "Step 2: bader CHGCAR -ref CHGCAR_sum"
write(*,*) "After this, the file ACF.dat should be present, besides POSCAR"
write(*,*) "Now, start this program. "
write(*,*) "NOTE: Currently, this program always assumes that the standard "
write(*,*) " POTCAR files were used for the calculation (no _h, _s,_pv etc)"
write(*,*) " If other POTCARs are used for certain elements, the number of "
write(*,*) " valence electrons must be given for the list of elements being"
write(*,*) " being affected with the command line argument:"
write(*,*) " -valence:el1=val1,el2=val2,..., example: -valence:Ga=13,Pt=16"
write(*,*) "To give an overview about utils4VASP, add the -overview command."
write(*,*) 
!
!    Read in custom element valences if needed from -valence argument
!
allocate(val_vec(10))
val_vec="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-valence:") then
!
!    Ignore bad readstat at this point...
!
      read(arg(10:),*,iostat=readstat) val_vec
      if (readstat .ne. 0) then
      end if
   end if
end do
mod_num=0
do i=1,10
   if (val_vec(i) .eq. "XX") exit
   mod_num=mod_num+1
end do
allocate(el_mods(mod_num))
allocate(charge_mod(mod_num))
do i=1,mod_num
   if (val_vec(i)(2:2) .eq. "=") then
      read(val_vec(i)(1:1),*,iostat=readstat1) el_mods(i) 
      read(val_vec(i)(3:5),*,iostat=readstat2) charge_mod(i)
      if (readstat1 .ne. 0 .or. readstat2 .ne. 0) then
         write(*,*) "Something went wrong in the -valence command for element",i,"!"
         stop
      end if        
   else if (val_vec(i)(3:3) .eq. "=") then    
      read(val_vec(i)(1:2),*,iostat=readstat1) el_mods(i)
      read(val_vec(i)(3:5),*,iostat=readstat2) charge_mod(i)
      if (readstat1 .ne. 0 .or. readstat2 .ne. 0) then
         write(*,*) "Something went wrong in the -valence command for element",i,"!"
         stop
      end if         
   end if        
end do

!
!    Open the POSCAR file, check if its there
!
open(unit=56,file="POSCAR",status="old",iostat=readstat)

if (readstat .ne. 0) then
   write(*,*) "ERROR! The file POSCAR is not there!"
   stop
end if
       
read(56,*)
read(56,*) cell_scale
read(56,*) a_vec(:)
read(56,*) b_vec(:)
read(56,*) c_vec(:)
read(56,'(a)') adum

!
!    Determine number of different elements in file 
!
do i=10,1,-1
   allocate(el_init(i))
   read(adum,*,iostat=readstat) el_init(:)
   if (readstat .eq. 0) then
      el_num=i
      allocate(el_names(el_num))
      el_names(1:el_num)=el_init(1:el_num)
      exit
   end if        
   deallocate(el_init)
   
end do

allocate(el_numbers(el_num))
read(56,*) el_numbers
natoms=sum(el_numbers)
read(56,*) adum
!
!    Read coordinate section for final PDB printout of structure
!
coord_direct=.false.
if (adum .eq. "Selective" .or. adum .eq. "selective") then
   read(56,*) adum
   if (adum .eq. "Direct" .or. adum .eq. "direct") then
      coord_direct=.true.
   end if   
else if (adum .eq. "Direct" .or. adum .eq. "direct") then
   coord_direct=.true.
end if        
allocate(xyz(3,natoms))
allocate(xyz_print(3,natoms))
do i=1,natoms
   read(56,*) xyz(:,i) 
end do
close(56)
!
!    If needed, convert coordinates to cartesian format
!
if (coord_direct) then
   do i=1,natoms
      xyz_print(:,i)=xyz(1,i)*a_vec(:)+xyz(2,i)*b_vec(:)+xyz(3,i)*c_vec(:)
   end do
else 
   xyz_print(:,:)=xyz(:,:)
end if        
write(*,'(a,i7)') " Number of elements in POSCAR:",el_num
!
!    Determine reference charges of elements 
!    Always assume the simple POTCAR files, no _sv,_pv etc
!
allocate(charge_ref(el_num))
charge_ref=0.d0
!
!    Overwrite valencies with -valence entries if given
!

do i=1,el_num
   do j=1,mod_num
      if (el_names(i) .eq. el_mods(j)) then
         charge_ref(i)=charge_mod(j)
      end if
   end do
end do

do i=1,el_num
!
!     Cycle if the charge of this element already has been defined
!
   if (charge_ref(i) .gt. 0.001d0) cycle
   if (el_names(i) .eq. "H") then
      charge_ref(i)=1.d0
   else if (el_names(i) .eq. "He") then
      charge_ref(i)=2.0d0
   else if (el_names(i) .eq. "Li") then
      charge_ref(i)=1.0d0
   else if (el_names(i) .eq. "Be") then
      charge_ref(i)=2.0d0
   else if (el_names(i) .eq. "B") then
      charge_ref(i)=3.0d0
   else if (el_names(i) .eq. "C") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "N") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "O") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "F") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Ne") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Na") then
      charge_ref(i)=1.0d0
   else if (el_names(i) .eq. "Mg") then
      charge_ref(i)=2.0d0
   else if (el_names(i) .eq. "Al") then
      charge_ref(i)=3.0d0
   else if (el_names(i) .eq. "Si") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "P") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "S") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "Cl") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Ar") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Sc") then
      charge_ref(i)=3.0d0
   else if (el_names(i) .eq. "Ti") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "V") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "Cr") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "Mn") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Fe") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Co") then
      charge_ref(i)=9.0d0
   else if (el_names(i) .eq. "Ni") then
      charge_ref(i)=10.0d0
   else if (el_names(i) .eq. "Cu") then
      charge_ref(i)=11.0d0
   else if (el_names(i) .eq. "Zn") then
      charge_ref(i)=12.0d0
   else if (el_names(i) .eq. "Ga") then
      charge_ref(i)=3.0d0
   else if (el_names(i) .eq. "Ge") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "As") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "Se") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "Br") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Kr") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Mo") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "Tc") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Ru") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Rh") then
      charge_ref(i)=9.0d0
   else if (el_names(i) .eq. "Pd") then
      charge_ref(i)=10.0d0
   else if (el_names(i) .eq. "Ag") then
      charge_ref(i)=11.0d0
   else if (el_names(i) .eq. "Cd") then
      charge_ref(i)=12.0d0
   else if (el_names(i) .eq. "In") then
      charge_ref(i)=3.0d0
   else if (el_names(i) .eq. "Sn") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "Sb") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "Te") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "I") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Xe") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "La") then
      charge_ref(i)=11.0d0
   else if (el_names(i) .eq. "Ce") then
      charge_ref(i)=12.0d0
   else if (el_names(i) .eq. "Pr") then
      charge_ref(i)=13.0d0
   else if (el_names(i) .eq. "Nd") then
      charge_ref(i)=14.0d0
   else if (el_names(i) .eq. "Pm") then
      charge_ref(i)=16.0d0
   else if (el_names(i) .eq. "Sm") then
      charge_ref(i)=16.0d0
   else if (el_names(i) .eq. "Eu") then
      charge_ref(i)=17.0d0
   else if (el_names(i) .eq. "Gd") then
      charge_ref(i)=18.0d0
   else if (el_names(i) .eq. "Tb") then
      charge_ref(i)=19.0d0
   else if (el_names(i) .eq. "Dy") then
      charge_ref(i)=20.0d0
   else if (el_names(i) .eq. "Ho") then
      charge_ref(i)=21.0d0
   else if (el_names(i) .eq. "Er") then
      charge_ref(i)=22.0d0
   else if (el_names(i) .eq. "Tm") then
      charge_ref(i)=23.0d0
   else if (el_names(i) .eq. "Yb") then
      charge_ref(i)=24.0d0
   else if (el_names(i) .eq. "Lu") then
      charge_ref(i)=25.0d0
   else if (el_names(i) .eq. "Hf") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "Ta") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "W") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "Re") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Os") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Ir") then
      charge_ref(i)=9.0d0
   else if (el_names(i) .eq. "Pt") then
      charge_ref(i)=10.0d0
   else if (el_names(i) .eq. "Au") then
      charge_ref(i)=11.0d0
   else if (el_names(i) .eq. "Hg") then
      charge_ref(i)=12.0d0
   else if (el_names(i) .eq. "Tl") then
      charge_ref(i)=3.0d0
   else if (el_names(i) .eq. "Pb") then
      charge_ref(i)=4.0d0
   else if (el_names(i) .eq. "Bi") then
      charge_ref(i)=5.0d0
   else if (el_names(i) .eq. "Po") then
      charge_ref(i)=6.0d0
   else if (el_names(i) .eq. "At") then
      charge_ref(i)=7.0d0
   else if (el_names(i) .eq. "Rn") then
      charge_ref(i)=8.0d0
   else if (el_names(i) .eq. "Ac") then
      charge_ref(i)=11.0d0
   else if (el_names(i) .eq. "Th") then
      charge_ref(i)=12.0d0
   else if (el_names(i) .eq. "Pa") then
      charge_ref(i)=13.0d0
   else if (el_names(i) .eq. "U") then
      charge_ref(i)=14.0d0
   else if (el_names(i) .eq. "Np") then
      charge_ref(i)=15.0d0
   else if (el_names(i) .eq. "Pu") then
      charge_ref(i)=16.0d0
   else if (el_names(i) .eq. "Am") then
      charge_ref(i)=17.0d0
   else if (el_names(i) .eq. "Cm") then
      charge_ref(i)=18.0d0
   else 
      write(*,*) "ERROR! No reference charge for element ",el_names(i)," available!"
      write(*,*) "Give it via the -valence command line argument!"
      stop
   end if        
end do
write(*,*) " Used number of valence electrons for the included elements:"
do i=1,el_num
   write(*,'(3a,f12.6)') "   -",el_names(i),":",charge_ref(i)
end do

write(*,*) 
write(*,*) "Calculate the charges ..."
write(*,*)
!
!    Now open ACF.dat file 
!

allocate(charge_avg(el_num))
charge_avg=0.d0
open(unit=70,file="ACF.dat",status="old")
read(70,*)
read(70,*)
allocate(charge_list(natoms))
allocate(el_atoms(natoms))
k=1
do i=1,el_num
   do j=1,el_numbers(i)
      read(70,*) idum,rdum,rdum,rdum,charge_act
      charge_avg(i)=charge_avg(i)+(charge_ref(i)-charge_act)     
      charge_list(k)=(charge_ref(i)-charge_act) 
      el_atoms(k)=el_names(i)
      k=k+1
   end do   
   charge_avg(i)=charge_avg(i)/el_numbers(i)
end do

close(70)

!
!    Print out individual Bader partial charges to separate file
!
open(unit=19,file="bader_charges.dat")
write(19,*) "# Bader charges of all atoms, written by eval_bader"
write(19,*) "# index      charge(e)"
do i=1,natoms
   write(19,'(i9,f15.8)') i, charge_list(i)
end do

close(19)
write(*,*) "Bader charges of atoms written to 'bader_charges.dat'."

!
!    Print out structure and charges together to file charges.pdb
!    (can be used for VMD visualization)
!
open(unit=20,file="charges.pdb")
write(20,'(a)') "COMPND    FINAL HEAT OF FORMATION =     0.000000"
write(20,'(a)') "AUTHOR    GENERATED BY PROGRAM EVAL BADER"
do i=1,natoms
   write(20,'(a,i5,a,a,a,3f8.3,f7.3,f5.2,a,a)') "HETATM",i," ",el_atoms(i), &
                 & "   UNL     1    ",xyz_print(:,i),charge_list(i),0d0,"          ",el_atoms(i)
end do
write(20,*) "MASTER        0    0    0    0    0    0    0    0  180    0  180    0"
write(20,*) "END"
close(20)
write(*,*) "File with coordinates and charges written to 'charges.pdb'"
write(*,*) " Open this file with VMD and select 'coloring method: occupancy'"

open(unit=21,file="POSCAR_charge")
write(21,*) "POSCAR file with charges, written by eval_bader"
write(21,*) cell_scale
write(21,*) a_vec(:)
write(21,*) b_vec(:)
write(21,*) c_vec(:)
do i=1,el_num
   write(21,'(a,a)',advance="no") " ",el_names(i)
end do
write(21,*)
write(21,*) el_numbers(:)
write(21,*) "Cartesian"
do i=1,natoms
   write(21,'(4f23.15)') xyz_print(:,i),charge_list(i)
end do
close(21)
write(*,*) "File with coordinates and charges written to 'POSCAR_charge'"
write(*,*)
!
!    Print out resulting average charges 
!

write(*,*) "The resulting average charges are:"
do i=1,el_num
   write(*,'(3a,f12.6)') "   -",el_names(i),":  ",charge_avg(i)
end do

write(*,*)
write(*,*) "eval_bader terminated normally."
write(*,*)
end program eval_bader      
