!
!    analyze_dft: Analyzes, evaluates and manages VASP
!      DFT single point calculations. Currenty supported
!      are the following calculation types:
!      - Bader partial charges
!      - Core level shift calculations
!      - (partial) density of state calculations
!      - STM image calculations
!    Part of VASP4CLINT
!     Julien Steffen, 2025 (julien.steffen@fau.de)
!
program analyze_dft
implicit none 
integer::i
character(len=120)::cdum,arg
logical::switch_bader,switch_stm,switch_pdos,switch_cls
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
      write(*,*) " - gen_poscar.py: Generate POSCAR of alloy and surface structures"
      write(*,*) " - modify_poscar.py: Modify POSCAR: multiply, transform, shift, insert etc."
      write(*,*) " - build_adsorb.py: Place adsorbates on surfaces, set translation, rotation"
      write(*,*) " - split_freq: Divide frequency calculations for large molecules"
      write(*,*) "Evaluation:"
      write(*,*) " - modify_xdatcar: Modify trajectory files: multiply, shift, pick etc."
      write(*,*) " - analyze_md: Analyze MD trajectories for RDFs, diffusion, density etc."
      write(*,*) " - analyze_dft: Analyze DFT calculations (Bader charges, STM, CLS, pDOS)"
      write(*,*) " - check_geoopt.py: Monitor geometry optimizations with selective dynamics"
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
write(*,*)
write(*,*) "PROGRAM analyze_dft: Analyzes evaluates and manages VASP"
write(*,*) " DFT single point calculations."
write(*,*) "Currently, four different calculation types can be treated."
write(*,*) "Select one of them by adding the respective command line "
write(*,*) " argument to access the respective subroutine with further"
write(*,*) " details and more detailed commands."
write(*,*) "The following command line arguments can be given:"
write(*,*) "  -bader: Evaluates a Bader charge calculation "
write(*,*) "  -stm: generates STM images from surface calculations."
write(*,*) "  -pdos: partial density of states calculations are evaluated"
write(*,*) "  -cls: Setup and evaluate core level shift calculations for"
write(*,*) "    several atoms in a structure."
write(*,*) "To give an overview about utils4VASP, add the -overview command."


!
!    If a Bader charge calculation shall be evaluated
!
switch_bader = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-bader") then
      switch_bader = .true.
   end if
end do

!
!    If a STM calculation shall be evaluated
!
switch_stm = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:4))  .eq. "-stm") then
      switch_stm = .true.
   end if
end do

!
!    If a partial density of states calculation shall be evaluated
!
switch_pdos = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:5))  .eq. "-pdos") then
      switch_pdos = .true.
   end if
end do

!
!    If a core level energy/shift calculation shall be evaluated
!
switch_cls = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:4))  .eq. "-cls") then
      switch_cls = .true.
   end if
end do

if ((.not. switch_bader) .and. (.not. switch_stm) .and. &
    & (.not. switch_pdos) .and. (.not. switch_cls)) then
   write(*,*)
   write(*,*) "Please give one of the calculation types!"
   write(*,*)
   stop
end if

!
!    Now call the respective subroutines
!
if (switch_bader) then
   call calc_bader
end if

if (switch_stm) then
   call calc_stm
end if

if (switch_pdos) then
   call calc_pdos
end if

if (switch_cls) then
   call calc_cls
end if

!
!   Print ending message
!
write(*,*)
write(*,*) "analyze_dft has finished all tasks, goodbye!"
write(*,*)

end program analyze_dft


!
!    calc_bader: evaluates the result of a Bader charge calculation
!      after preprocessing the VASP calculation results with the
!      chgsum.pl and bader programs by the Henkelman group.
!      The files POSCAR and ACF.dat need to be present
!

subroutine calc_bader
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
!    Print general information and all possible keywords of the program
!
write(*,*)
write(*,*) "Option -bader: evaluation of Bader charge calculations"
write(*,*) " from preprocessed VASP output for arbitrary systems."
write(*,*) "Before starting the analysis, the files AECCAR0, AECCAR2 and "
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

end subroutine calc_bader


!
!    calc_stm: evaluates the results of a STM calculation
!      (partial charge density, to be evaluated by the Tersoff-Hamann
!      approach). A PARCHG file is read in and depending on the
!      user keywords, a constant current or constant height STM
!      picture is evaluated. Both a bitmap file and a plottable
!      data file with a gnuplot file to plot it are given as
!      results.
!

subroutine calc_stm
implicit none
integer::readstat
integer::i,j,k,l,m,n
integer::k_act,l_act,m_act
integer::natoms,nelems
integer::grida,gridb,gridc
integer::nx,ny,nz
integer::stm_gridx,stm_gridy
integer::near_x,near_y,near_z
integer::repeat_x,repeat_y
integer::include_x,include_y,include_z
real(kind=8)::scale
real(kind=8)::a_fac,cutoff,g_width
real(kind=8)::x_len,y_len,z_len
real(kind=8)::stm_height,isos_dens
real(kind=8)::x_plot_min,x_plot_max
real(kind=8)::y_plot_min,y_plot_max
real(kind=8)::dist_act,xtoy
real(kind=8)::weight_sum,weight
real(kind=8)::prev_z,prev_z_init
real(kind=8)::def_z_shift,z_step
real(kind=8)::stm_act
real(kind=8)::coord_mat(3,3)
real(kind=8)::act_pos(3)
real(kind=8)::vec_act(3)
real(kind=8)::diff_vec(3)
real(kind=8),allocatable::xyz(:,:)
real(kind=8),allocatable::cdens(:,:,:)
real(kind=8),allocatable::c_coord(:,:,:,:)
real(kind=8),allocatable::stm_dat(:,:)
real(kind=8),allocatable::stm_pos(:,:,:)
logical::eval_stat(10)
character(len=120)::cdum,arg
character(len=30)::stm_mode
character(len=2),allocatable::el_names_read(:),el_names(:)
character(len=2),allocatable::at_names(:)
integer,allocatable::el_nums(:)


!
!    Print general information and all possible keywords of the program
!
write(*,*)
write(*,*) "Option -stm: Plots a STM image from a VASP"
write(*,*) " PARCHG file based on the Tersoff-Hamann approach."
write(*,*) "Constant height and constant current images are possible."
write(*,*) "The only needed input is a PARCHG file for a certain"
write(*,*) " energy range below the Fermi level, corresponding to an"
write(*,*) " experimental STM tunneling voltage."
write(*,*) "The following command line arguments can be given:"
write(*,*) "  -mode=[height or current] : Determines if the constant"
write(*,*) "    height or constant current STM mode shall be used."
write(*,*) "  -pos=[value] : For constant height STMs: the position of"
write(*,*) "    imaginary tip along the z-axis of the unit cell."
write(*,*) "    The value needs to be between 0 and the height of the"
write(*,*) "    unit cell! (in Angstroms)"
write(*,*) "  -dens=[value] : For constant current STMs : The charge "
write(*,*) "    density at which the tip shall be located at each grid"
write(*,*) "    point along x and y, resulting in a z-value (height)."
write(*,*) "    Reasonable values are between 0.01 and 0.1"
write(*,*) "  -repeat_x=[number] : For the generated picture: the number"
write(*,*) "    of unit cell repetitions along x-axis to generate a larger  "
write(*,*) "    picture which shows the periodicity. (DEFAULT: 1)"
write(*,*) "  -repeat_y=[number] : For the generated picture: the number"
write(*,*) "    of unit cell repetitions along y-axis to generate a larger  "
write(*,*) "    picture which shows the periodicity. (DEFAULT: 1)"
write(*,*) "  -grid_x=[number] : The number of grid points along x-axis "
write(*,*) "    in the unit cell at which the STM picture shall be "
write(*,*) "    calculated. (DEFAULT: 100)"
write(*,*) "  -grid_y=[number] : The number of grid points along y-axis "
write(*,*) "    in the unit cell at which the STM picture shall be "
write(*,*) "    calculated. (DEFAULT: 100)"
write(*,*) "  -gauss_width=[value] : The full width of half maximum of the"
write(*,*) "    Gaussian that is used to smooth the image by including the "
write(*,*) "    charge densities around an evaluation point by weighting them"
write(*,*) "    with the Gaussian function value. The larger the value, the"
write(*,*) "    smoother/smeared the picture will be. (DEFAULT: 0.3 Angs.)"
write(*,*)
!
!    For the constant current mode: the moving up of the tip for each
!    new row position (in Angstroms)
!
def_z_shift=2.d0
!
!    Read command line arguments
!
!    If the constant height or the constant current mode
!      shall be used
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-mode=") then
      read(arg(7:),*,iostat=readstat) stm_mode
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -mode=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "height") then
   write(*,*) "The constant height mode will be evaluated."
else if (trim(stm_mode) .eq. "current") then
   write(*,*) "The constant current mode will be evaluated."
else
   write(*,*) "Please give either the constant height mode (-mode=height) "
   write(*,*) " or the constant current mode (-mode=current)!"
   stop
end if
!
!    If constant height mode: read in the height (z/c coordinate)
!
stm_height=-1.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:5))  .eq. "-pos=") then
      read(arg(6:),*,iostat=readstat) stm_height
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -pos=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "height") then
   if (stm_height .lt. 0.d0) then
      write(*,*) "Please give the height (z-coordinate) at which the STM"
      write(*,*) " tip shall scan the surface!"
      stop
   else
      write(*,'(a,f13.6,a)') " The STM tip will scan the surface at z=", &
                      & stm_height," Angstroms"
   end if
end if
!
!    If constant current mode: read in the isos density
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-dens=") then
      read(arg(7:),*,iostat=readstat) isos_dens
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -dens=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "current") then
   if (isos_dens .lt. 0.d0) then
      write(*,*) "Please give a positive density (e/Ang^3) to which the "
      write(*,*) " shall be moved down in each grid point!"
      stop
   else
      write(*,'(a,f10.6,a)') " The local density to which the STM tip will be dropped is ", &
                      & isos_dens," e/Ang.^3"
   end if
end if

!
!     The number of STM image repetitions along the x axis
!
!     Default (only actual unit cell)
!
repeat_x=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-repeat_x=") then
      read(arg(11:),*,iostat=readstat) repeat_x
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -repeat_x=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i3)') " Number of unit cell repetitions along x: ",repeat_x

!
!     The number of STM image repetitions along the y axis
!
!     Default (only actual unit cell)
!
repeat_y=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-repeat_y=") then
      read(arg(11:),*,iostat=readstat) repeat_y
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -repeat_y=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i3)') " Number of unit cell repetitions along y: ",repeat_y
!
!     The number of grid points for STM printout along the x axis
!
!     Default
!
stm_gridx=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-grid_x=") then
      read(arg(9:),*,iostat=readstat) stm_gridx
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -grid_x=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i6)') " Number of plotted STM grid points (unit cell) along x: ",stm_gridx

!
!     The number of grid points for STM printout along the y axis
!
!     Default
!
stm_gridy=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-grid_y=") then
      read(arg(9:),*,iostat=readstat) stm_gridy
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -grid_y=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i6)') " Number of plotted STM grid points (unit cell) along y: ",stm_gridy

!
!     The width at half maximum (FWHM) of the Gaussian used for
!       smoothing the picture
!
!     Default
!
g_width=0.3d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-gauss_width=") then
      read(arg(14:),*,iostat=readstat) g_width
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -gauss_width=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,f12.7)') " Width of the smooting Gaussian at half height (Angs.): ",g_width
!
!    Determine the properties and the function of the smoothing Gaussian
!
a_fac=dlog(2.d0)/((g_width/2d0)**2)
cutoff=sqrt(dlog(100d0)/a_fac)

write(*,*)
write(*,*) "Read in the PARCHG file ..."
!
!    Open file PARCHG and read in relevant data
!
open(unit=15,file="PARCHG",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file PARCHG is not there!"
   stop
end if
!
!    Read the unit cell definition
!
read(15,*)
read(15,*) scale
read(15,*) coord_mat(1,:)
read(15,*) coord_mat(2,:)
read(15,*) coord_mat(3,:)
!
!    The lengths of coordinate axes
!
x_len=sqrt(dot_product(coord_mat(1,:),coord_mat(1,:)))
y_len=sqrt(dot_product(coord_mat(2,:),coord_mat(2,:)))
z_len=sqrt(dot_product(coord_mat(3,:),coord_mat(3,:)))
!
!    Read in the elements
!
allocate(el_names_read(10))
el_names_read="XX"
read(15,'(a)') cdum
read(cdum,*,iostat=readstat) el_names_read
nelems=0
do i=1,10
   if (el_names_read(i) .eq. "XX") exit
   nelems=nelems+1
end do
allocate(el_names(nelems),el_nums(nelems))

do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(15,*) el_nums

!
!    Define the element symbols for all atoms in the system
!
natoms = sum(el_nums)
allocate(at_names(natoms))
k=0
do i=1,nelems
   do j=1,el_nums(i)
      k=k+1
      at_names(k)=el_names(i)
   end do
end do
!
!    Read in the atomic coordinates
!
allocate(xyz(3,natoms))
read(15,*)
do i=1,natoms
   read(15,*) xyz(:,i)
end do
!
!    Read in the number of grid points
!
read(15,*)
read(15,*) grida,gridb,gridc
allocate(cdens(grida,gridb,gridc))
!
!    Read in the partial charge density
!    Use compactified read command for array
!
read(15,*) (((cdens(nx,ny,nz),nx=1,grida),ny=1,gridb),nz=1,gridc)
close(15)
write(*,*) "  done!"
write(*,'(a)') " Length of coordinate axes (Angstroms):"
write(*,'(a,f12.6,a,f12.6,a,f12.6,a)') "   x:",x_len,", y:",y_len,", z:",z_len
write(*,'(a,i6,a,i6,a,i6,a,i9,a)') " Number of partial density grid points:"
write(*,'(a,i6,a,i6,a,i6,a,i9,a)') "   x:",grida, &
              & ", y:",gridb,", z:",gridc, "   (total:",grida*gridb*gridc,")"
write(*,*)
!
!    Calculate x,y,z cartesian coordinates for all grid points, for the
!    subsequent determination of distances to the current sampling point
!
allocate(c_coord(3,grida,gridb,gridc))
do i=1,gridc
   do j=1,gridb
      do k=1,grida
         vec_act(1)=real((k-1))/real(grida)
         vec_act(2)=real((j-1))/real(gridb)
         vec_act(3)=real((i-1))/real(gridc)
         c_coord(:,k,j,i)=vec_act(:)!matmul(vec_act,coord_mat)
      end do
   end do
end do
!
!    Allocate data array for 2D STM picture
!
allocate(stm_dat(stm_gridy,stm_gridx))
allocate(stm_pos(2,stm_gridy,stm_gridx))
stm_dat=0.d0
!
!    Determine number of grid points to be included into the Gaussian
!    smoothin region
!    Along z only half the length! (since we do not probe into depth)
!
include_x=nint(real(cutoff)/real(x_len)*real(grida))
include_y=nint(real(cutoff)/real(y_len)*real(gridb))
include_z=nint(real(cutoff)/real(z_len)*real(gridc)/2)
write(*,*) "Generate the STM picture at all plot grid points ...  "
!
!    A: CONSTANT HEIGHT MODE
!
if (trim(stm_mode) .eq. "height") then
!
!    Loop over all picture gridpoints in the x-y (or a-b) plane and calculate
!    the Gaussian-weighted interpolation of nearby grid points
!
   eval_stat = .false.
   do i=1,stm_gridx
!
!    Every 10% of the read in, give a status update
!
      do j=1,10
         if (real(i)/real(stm_gridx) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      do j=1,stm_gridy
         act_pos(1)=real((i-1))/real(stm_gridx)
         act_pos(2)=real((j-1))/real(stm_gridy)
         act_pos(3)=stm_height/z_len
!
!    Determine nearest point in density grid
!
         near_x=nint(act_pos(1)*grida)
         near_y=nint(act_pos(2)*gridb)
         near_z=nint(act_pos(3)*gridc)
       !  act_pos(:)=matmul(act_pos,coord_mat)

!
!    Now loop over all other grid points within the proposed cutoff of the
!    Gaussian smearing function
!
         weight_sum=0
         do k=near_x-include_x,near_x+include_x
            if (k .lt. 1) then
               k_act=grida+k
            else if (k .gt. grida) then
               k_act=k-grida
            else
               k_act=k
            end if
            do l=near_y-include_y,near_y+include_y
               if (l .lt. 1) then
                  l_act=gridb+l
               else if (l .gt. gridb) then
                  l_act=l-gridb
               else
                  l_act=l
               end if
               do m=near_z-include_z,near_z+include_z
                  if (m .lt. 1) then
                     m_act=gridc+m
                  else if (m .gt. gridc) then
                     m_act=m-gridc
                  else
                     m_act=m
                  end if
!
!    Now sum over all grid points within the chosen cutoff for consideration
!    and add the components together
!
                !  write(*,*) k_act,l_act,m_act
                  diff_vec(:)=act_pos(:)-c_coord(:,k_act,l_act,m_act)

                  do n=1,3
                     if (diff_vec(n) .gt. 0.5d0) then
                        diff_vec(n) = diff_vec(n)-0.5d0
                     end if
                     if (diff_vec(n) .lt. 0.d0) then
                        diff_vec(n) = diff_vec(n)+0.5d0
                     end if
                  end do

                  diff_vec=matmul(diff_vec,coord_mat)
                  dist_act=sqrt(dot_product(diff_vec,diff_vec))
                  weight=exp(-g_width*dist_act*dist_act)
                  weight_sum=weight_sum+weight
                  stm_dat(j,i)=stm_dat(j,i)+weight*cdens(k_act,l_act,m_act)
               end do
            end do
         end do
         stm_dat(j,i)=stm_dat(j,i)/weight_sum
         act_pos=matmul(act_pos,coord_mat)
         stm_pos(:,j,i)=act_pos(1:2)

      end do
   end do
end if
!
!    B: CONSTANT CURRENT MODE
!
if (trim(stm_mode) .eq. "current") then
!
!    Loop over all picture gridpoints in the x-y (or a-b) plane and calculate
!    the Gaussian-weighted interpolation of nearby grid points
!    For constant current, start for the first point of each row at 0.8 times the
!    maximum z-value of the cell and scan below until a value equal or higher the
!    desired current (local charge density) is reached, similar to a STM tip that
!    is moved down to a surface
!    Within each row, from one point to the next, only scan the next +/- 1 Angstrom
!    relative to the previous point, assuming that the charge density won't change too
!    much between two gridpoints
!
!    The default movement along z for scanning the correct height, a fraction of the
!    z-grid density within the given PARCHG file.
!
   eval_stat = .false.
   z_step=1.d0/(gridc*3d0)
   do i=1,stm_gridx
!
!    Every 10% of the read in, give a status update
!
      do j=1,10
         if (real(i)/real(stm_gridx) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      if (i .eq. 1) then
         prev_z=0.8d0-def_z_shift/z_len
      end if
!
!    Now loop trough the remaining points of the current row and start with the init z value
!    of this row
!
!      prev_z=prev_z_init
      do j=1,stm_gridy
         act_pos(1)=real((i-1))/real(stm_gridx)
         act_pos(2)=real((j-1))/real(stm_gridy)
         act_pos(3)=prev_z+def_z_shift/z_len
!
!    Scan the starting position for a new grid point
!
         do
            if (act_pos(3) .ge. 1.d0) then
               write(*,*) "The tip would be located outside the cell!"
               write(*,*) "Increase the -dens parameter if possible!"
               stop
            end if
!
!    Determine nearest point in density grid
!
            near_x=nint(act_pos(1)*grida)
            near_y=nint(act_pos(2)*gridb)
            near_z=nint(act_pos(3)*gridc)
!
!    Now loop over all other grid points within the proposed cutoff of the
!    Gaussian smearing function
!
            weight_sum=0.d0
            stm_act=0.d0
            do k=near_x-include_x,near_x+include_x
               if (k .lt. 1) then
                  k_act=grida+k
               else if (k .gt. grida) then
                  k_act=k-grida
               else
                  k_act=k
               end if
               do l=near_y-include_y,near_y+include_y
                  if (l .lt. 1) then
                     l_act=gridb+l
                  else if (l .gt. gridb) then
                     l_act=l-gridb
                  else
                     l_act=l
                  end if
                  do m=near_z-include_z,near_z+include_z
                     if (m .lt. 1) then
                        m_act=gridc+m
                     else if (m .gt. gridc) then
                        m_act=m-gridc
                     else
                        m_act=m
                     end if
!
!    Now sum over all grid points within the chosen cutoff for consideration
!    and add the components together
!
                     diff_vec(:)=act_pos(:)-c_coord(:,k_act,l_act,m_act)

                     do n=1,3
                        do while (abs(diff_vec(n)) .gt. 0.5d0)
                           diff_vec(n) = diff_vec(n) -sign(1.0d0,diff_vec(n))
                        end do
                     end do
                     diff_vec=matmul(diff_vec,coord_mat)
                     dist_act=sqrt(dot_product(diff_vec,diff_vec))
                     weight=exp(-g_width*dist_act*dist_act)
                     weight_sum=weight_sum+weight
                     stm_act=stm_act+cdens(k_act,l_act,m_act)
                  end do
               end do
            end do
!
!    Leave the scanning loop if either the corrent density has been found or
!    the tip reached the ground of the unit cell
!
            stm_act=stm_act/weight_sum
            if ((stm_act .ge. isos_dens) .or. (act_pos(3) .lt. 0.01d0)) then
               act_pos=matmul(act_pos,coord_mat)
               stm_pos(:,j,i)=act_pos(1:2)
               stm_dat(j,i)=act_pos(3)
               prev_z=act_pos(3)/z_len
               exit
            end if
            act_pos(3)=act_pos(3)-z_step
         end do
      end do
   end do
end if

write(*,*) " done!"
!
!    Write plot file with STM data
!
write(*,*) "Write file with 2D plot data for STM image (stm_plot.dat) ..."
open(unit=37,file="stm_plot.dat",status="replace")
write(37,*) "# This data has been generated by the program eval_stm."
write(37,*) "# x-coord     y-coord     intenstiy  "
x_plot_min=100.d0
x_plot_max=-100.d0
y_plot_min=100.d0
y_plot_max=-100.d0
do k=1,repeat_x
   do i=1,stm_gridx
      do l=1,repeat_y
         do j=1,stm_gridy
            vec_act(1)=real(k-1)
            vec_act(2)=real(l-1)
            vec_act(3)=1.d0
            vec_act=matmul(vec_act,coord_mat)
            vec_act(1)=vec_act(1)+stm_pos(1,j,i)
            vec_act(2)=vec_act(2)+stm_pos(2,j,i)
            if (vec_act(1) .lt. x_plot_min) then
               x_plot_min=vec_act(1)
            else if (vec_act(1) .gt. x_plot_max) then
               x_plot_max=vec_act(1)
            end if
            if (vec_act(2) .lt. y_plot_min) then
               y_plot_min=vec_act(2)
            else if (vec_act(2) .gt. y_plot_max) then
               y_plot_max=vec_act(2)
            end if
            write(37,*) vec_act(1), vec_act(2), stm_dat(j,i)
         end do
      end do
      write(37,*)
   end do
end do
close(37)
write(*,*) " done!"
!
!    Write gnuplot file for the plot
!    Automatically align the plotted png picture to the proportion of
!     x and y axes lengths
!
write(*,*) "Write gnuplot file for 2D plot of STM image (stm_plot.gnu) ..."
xtoy=(x_plot_max-x_plot_min)/(y_plot_max-y_plot_min)
open(unit=38,file="stm_plot.gnu",status="replace")
write(38,*) "# This file generates a png picture from 'stm_plot.dat'."
write(38,*) "set terminal png size 2400,",nint(2000/xtoy)," lw 3.5 font 'Helvetica,46'"
write(38,*) "set output 'stm_plot.png'"
write(38,*) "set encoding iso_8859_1"
write(38,*) "set lmargin screen 0.10"
write(38,*) "set rmargin screen 0.80"
write(38,*) "set tmargin screen 0.13"
write(38,*) "set bmargin screen 0.95"
write(38,*) "set xrange [",x_plot_min,":",x_plot_max,"]"
write(38,*) "set yrange [",y_plot_min,":",y_plot_max,"]"
write(38,*) "set xlabel 'x-coordinate ({\305})"
write(38,*) "set ylabel 'y-coordinate ({\305})"
if (trim(stm_mode) .eq. "height") then
   write(38,*) "set cblabel 'local density (e/{\305}^3)' offset -.8,0"
end if
if (trim(stm_mode) .eq. "current") then
   write(38,*) "set cblabel 'height ({\305})' offset -.8,0"
end if
write(38,*) "set cblabel offset 0.8,0"
write(38,*) "set palette gray"
write(38,*) "set pm3d map interpolate 0,0"
write(38,*) "splot 'stm_plot.dat' with pm3d notitle"
close(38)

write(*,*) " done!"
write(*,*) "Execute gnuplot to obtain image (stm_plot.png) ..."
call system("gnuplot stm_plot.gnu")
write(*,*) " done!"
write(*,*)
!
!    Write second gnuplot file for plot without axis labels or tics
!
write(*,*) "Write gnuplot file for 2D plot wihtout axis labels (stm_plot_blank.gnu) ..."
xtoy=(x_plot_max-x_plot_min)/(y_plot_max-y_plot_min)
open(unit=38,file="stm_plot_blank.gnu",status="replace")
write(38,*) "# This file generates a png picture from 'stm_plot.dat'."
write(38,*) "set terminal png size 2400,",nint(2400/xtoy)," lw 3.5 font 'Helvetica,46'"
write(38,*) "set output 'stm_plot_blank.png'"
write(38,*) "set encoding iso_8859_1"
write(38,*) "set lmargin screen 0.0"
write(38,*) "set rmargin screen 1.0"
write(38,*) "set tmargin screen 0.0"
write(38,*) "set bmargin screen 1.0"
write(38,*) "set xrange [",x_plot_min,":",x_plot_max,"]"
write(38,*) "set yrange [",y_plot_min,":",y_plot_max,"]"
write(38,*) "unset xtics"
write(38,*) "unset ytics"
write(38,*) "unset cbtics"
write(38,*) "unset colorbox"

write(38,*) "set palette gray"
write(38,*) "set pm3d map interpolate 0,0"
write(38,*) "splot 'stm_plot.dat' with pm3d notitle"
close(38)

write(*,*) " done!"
write(*,*) "Execute gnuplot to obtain image without labels (stm_plot_blank.png) ..."
call system("gnuplot stm_plot_blank.gnu")
write(*,*) " done!"


end subroutine calc_stm


!
!    calc_pdos: evaluates DOSCAR files and prints out the DOS
!      plots for one or several elements, atom indices or orbitals
!
subroutine calc_pdos
implicit none
integer::i,j,k
integer::readstat,counter,ispin
integer::nelems,natoms,natoms_test,nelems_choice
character(len=120)::cdum,adum
character(len=2),allocatable::el_list(:),el_list_read(:)
character(len=5),allocatable::orb_list_read(:),orb_list(:)
character(len=100)::arg
character(len=2),allocatable::el_names_read(:),el_names(:),at_names(:)
integer,allocatable::el_nums(:)
real(kind=8)::enmax,enmin,efermi,fdum
real(kind=8),allocatable::e_dos(:),dos_tot(:),dos_tot_up(:),dos_tot_down(:)
real(kind=8),allocatable::dos_partial(:,:,:)
real(kind=8),allocatable::dos_part_out(:),dos_part_out_up(:),dos_part_out_down(:)
integer::npoints,ninds_choice,norb_choice
integer,allocatable::ind_list(:),ind_list_read(:)
logical::dos_elements,dos_indices,act_read,el_all,orb_all,use_orbital
logical,allocatable::orb_used(:)


!
!    Print general information and all possible keywords of the program
!
write(*,*)
write(*,*) "Option -pdos: Evaluation of VASP density of states (DOS)"
write(*,*) " calculations. Element- or atom-resolved DOS can be printed out,"
write(*,*) " further, one or several certain orbitals can be chosen."
write(*,*) "The POSCAR and DOSCAR files of a VASP calculation need to be "
write(*,*) " present, the calculation must be done with the LORBIT=11 command."
write(*,*) "The following command line arguments can/must be given:"
write(*,*) " -element=[list of elements] : which element shall be analyzed."
write(*,*) "    one or several can be given. State 'all' if all atoms shall be "
write(*,*) "    incorporated into the DOS evaluation. Examples: -element=Pt or "
write(*,*) "    -element=Pt,Ga,In or -element=all"
write(*,*) " -index=[one or more indices] : If certain atoms shall be analyzed."
write(*,*) "    One or more indices can be given. Examples: -index=103 or "
write(*,*) "    -index=13,103,201"
write(*,*) " -orbital=[name(s)] : Which orbitals shall be analyzed for the given "
write(*,*) "    elements or atom indices. Either 'all' for all orbitals or "
write(*,*) "    one or several descriptors need to be given, possible are:"
write(*,*) "    s,px,py,pz,dxy,dyz,dz2,dxz,dx2y2. As shortcuts can further be "
write(*,*) "    used: p (all p-orbitals), d (all d-orbitals).  "
write(*,*) "    Examples: -orbital=all or -orbital=s or -orbital=p,dxy,dxz"
write(*,*) "If both -element=all and -orbital=all are given, only the global DOS"
write(*,*) " will be written, since then the whole wavefunction is chosen."
write(*,*) "The Fermi energy will always be subtracted to build the x-axis. If"
write(*,*) " you want to plot with without the correction, use the Fermi energy"
write(*,*) " printed in the headers of the files."
write(*,*)
!
!    Read in and process the command line arguments
!

!
!    The elements whose DOS shall be calculated
!
dos_elements=.false.
el_all=.false.
allocate(el_list_read(10))
el_list_read="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-element=") then
      read(arg(10:),*,iostat=readstat) el_list_read
      dos_elements=.true.
   end if
end do

if (dos_elements) then
   nelems_choice=0
   do i=1,10
      if (el_list_read(i) .eq. "XX") exit
!
!    If all elements are chosen (all)
!
      if (el_list_read(i) .eq. "al") then
         el_all = .true.
         exit
      end if
      nelems_choice=nelems_choice+1
   end do

   allocate(el_list(nelems_choice))
   do i=1,nelems_choice
      el_list(i)=el_list_read(i)
   end do
   if (el_all) then
      write(*,*) "All atoms in the system were chosen for the DOS evalulation."
   else
      write(*,*) "The following elements were chosen for the DOS evaluation:"
      do i=1,nelems_choice
         write(*,*) " - ",el_list(i)
      end do
   end if
end if


!
!     The atom indices whose DOS shall be calculated
!

dos_indices=.false.
allocate(ind_list_read(1000))
ind_list_read=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-index=") then
      read(arg(8:),*,iostat=readstat) ind_list_read
      if (dos_elements) then
         write(*,*) "Both -element=... and -index=... flags were given."
         write(*,*) " Please use only one of them at a time!"
         stop
      end if
      dos_indices=.true.
   end if
end do

if (dos_indices) then
   ninds_choice=0

   do i=1,1000
      if (ind_list_read(i) .eq. 0) exit
      ninds_choice=ninds_choice+1
   end do

   allocate(ind_list(ninds_choice))
   do i=1,ninds_choice
      ind_list(i)=ind_list_read(i)
   end do

   write(*,*) "The following atoms were chosen for the DOS evaluation:"
   do i=1,ninds_choice
      write(*,*) " - No. ",ind_list(i)
   end do

end if

if ((.not. dos_indices) .and. (.not. dos_elements)) then
   write(*,*) "Please give either a -element=... or -index=... selection!"
   stop
end if

!
!    The orbital angular momentum specifiers
!
orb_all=.false.
allocate(orb_list_read(10))
orb_list_read="XXXXX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-orbital=") then
      read(arg(10:),*,iostat=readstat) orb_list_read
      use_orbital=.true.
   end if
end do
if (.not. use_orbital) then
   write(*,*) "Please give a -orbital=... selection!"
   stop
end if

norb_choice=0
do i=1,10
   if (orb_list_read(i) .eq. "XXXXX") exit
!
!    If all elements are chosen (all)
!
   if (orb_list_read(i) .eq. "all") then
      orb_all = .true.
      exit
   end if
   norb_choice=norb_choice+1
end do

allocate(orb_list(norb_choice))
do i=1,norb_choice
   orb_list(i)=orb_list_read(i)
end do

allocate(orb_used(9))
orb_used = .false.
if (orb_all) orb_used = .true.
do i=1,norb_choice
   if (orb_list(i) .eq. "s") orb_used(1) = .true.
   if (orb_list(i) .eq. "p") orb_used(2:4) = .true.
   if (orb_list(i) .eq. "d") orb_used(5:9) = .true.
   if (orb_list(i) .eq. "px") orb_used(2) = .true.
   if (orb_list(i) .eq. "py") orb_used(3) = .true.
   if (orb_list(i) .eq. "pz") orb_used(4) = .true.
   if (orb_list(i) .eq. "dxy") orb_used(5) = .true.
   if (orb_list(i) .eq. "dyz") orb_used(6) = .true.
   if (orb_list(i) .eq. "dz2") orb_used(7) = .true.
   if (orb_list(i) .eq. "dxz") orb_used(8) = .true.
   if (orb_list(i) .eq. "dx2y2") orb_used(9) = .true.
end do
write(*,*) "The following orbitals were chosen for the DOS evaluation:"
write(*,*) "  (T: used, F: not used)"
write(*,*) " - s     : ",orb_used(1)
write(*,*) " - px    : ",orb_used(2)
write(*,*) " - py    : ",orb_used(3)
write(*,*) " - pz    : ",orb_used(4)
write(*,*) " - dxy   : ",orb_used(5)
write(*,*) " - dyz   : ",orb_used(6)
write(*,*) " - dz2   : ",orb_used(7)
write(*,*) " - dxz   : ",orb_used(8)
write(*,*) " - dx2y2 : ",orb_used(9)
write(*,*)


!
!    First, open the POSCAR file for elements and atom indices
!

open(unit=28,file="POSCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   stop "The POSCAR file is not there!"
end if

read(28,*)
read(28,*)
read(28,*)
read(28,*)
read(28,*)
!
!    Determine the element symbols and their numbers
!
allocate(el_names_read(10))
el_names_read="XX"
read(28,'(a)') cdum

read(cdum,*,iostat=readstat) el_names_read
nelems=0

do i=1,10
   if (el_names_read(i) .eq. "XX") exit
   nelems=nelems+1
end do
!
!    Define permanent element symbol and number arrays and the total
!    number of atoms in the system
!
allocate(el_names(nelems),el_nums(nelems))
do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(28,*) el_nums

natoms = sum(el_nums)

!
!    Define individual elements for each atom
!

allocate(at_names(natoms))
counter = 1
do i=1,nelems
   do j=1,el_nums(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do

close(28)

!
!    Second, read in the DOSCAR completely
!

open(unit=29,file="DOSCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   stop "The DOSCAR file is not there!"
end if
read(29,*) natoms_test
if (natoms_test .eq. natoms) then
   write(*,*) "Check compatibility of POSCAR and DOSCAR ..."
   write(*,*) " ...  success!"
else
   write(*,*) "Check compatibility of POSCAR and DOSCAR ..."
   write(*,*) "POSCAR has ", natoms, " atoms, DOSCAR has ",natoms_test," atoms!"
   stop "Check POSCAR and DOSCAR!"
end if
read(29,*)
read(29,*)
read(29,*)
read(29,*)
read(29,*) enmax,enmin,npoints,efermi
write(*,*) "General information:"
write(*,'(a,f16.8,a)') "  - Minimum DOS energy: ",enmin," eV"
write(*,'(a,f16.8,a)') "  - Maximum DOS energy: ",enmax," eV"
write(*,'(a,f16.8,a)') "  - Fermi energy: ",efermi, " eV"
write(*,*)
allocate(e_dos(npoints))

!
!     Determine if the calculation was spin-polarized or not and read in the
!     total DOS (always written at the beginning of DOSCAR)
!
write(*,*) "Read in the total DOS ..."
ispin=1
do i=1,npoints
   if (i .eq. 1) then
      read(29,'(a)') adum
      read(adum,*,iostat=readstat) fdum,fdum,fdum,fdum,fdum
      if (readstat .eq. 0) then
         write(*,*) "You have done a spin-polarized calculation (ISPIN=2)!"
         write(*,*) "Two columns will be written for the DOS: up and down."
         ispin=2
         allocate(dos_tot_up(npoints),dos_tot_down(npoints))
         read(adum,*) e_dos(1),dos_tot_up(1),fdum,dos_tot_down(1)
      else
         write(*,*) "You have done no spin-polarized calculation (ISPIN=1)!"
         write(*,*) "The DOS will be written into one column."
         ispin=1
         allocate(dos_tot(npoints))
         read(adum,*) e_dos(1),dos_tot(1)
      end if
      cycle
   end if
   if (ispin .eq. 2) then
      read(29,*) e_dos(i),dos_tot_up(i),fdum,dos_tot_down(i)
   else if (ispin .eq. 1) then
      read(29,*) e_dos(i),dos_tot(i)
   end if
end do
write(*,*) " ...  success!"
!
!     Write the total DOS to a file
!
open(unit=30,file="dos_total.dat",status="replace")
write(30,*) "# DOS calculated by partial_dos "
write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"

if (ispin .eq. 1) then
   write(30,*) "#   energy(eV)             DOS"
else if (ispin .eq. 2) then
   write(30,*) "#   energy(eV)             DOS(up)            DOS(down)"
end if
do i=1,npoints
   if (ispin .eq. 1) then
      write(30,*) e_dos(i)-efermi,dos_tot(i)
   else if (ispin .eq. 2) then
      write(30,*) e_dos(i)-efermi,dos_tot_up(i),dos_tot_down(i)
   end if
end do
close(30)
write(*,*) "Total DOS of the system written to file 'dos_total.dat'."

!
!    If spin polarization is present, write an averaged DOS
!     (e.g., for effective calculation of DOS overlaps)
!
if (ispin .eq. 2) then
   open(unit=30,file="dos_tot_sum.dat",status="replace")
   write(30,*) "# DOS calculated by partial_dos (alpha and beta summed up) "
   write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"
   write(30,*) "#   energy(eV)             DOS"
   do i=1,npoints
      write(30,*) e_dos(i)-efermi,dos_tot_up(i)+dos_tot_down(i)
   end do
   close(30)
   write(*,*) "Summed total DOS (alpha+beta) written to file 'dos_tot_sum.dat'."
end if

if (el_all .and. orb_all) then
   write(*,*)
   write(*,*) "All atoms and all orbitals are chosen, the partial DOS will be skipped."
   goto 22
end if
!
!     Now read the individual DOS of all atoms and their respective orbitals
!
write(*,*)
write(*,*) "Read in the partial DOS ..."
if (ispin .eq. 1) then
   allocate(dos_partial(natoms,npoints,9))
else if (ispin .eq. 2) then
   allocate(dos_partial(natoms,npoints,18))
end if
do i=1,natoms
   read(29,*,iostat=readstat) fdum
   if (readstat .ne. 0) then
      stop "Missing data in DOSCAR? Did you add the LORBIT=11 keyword?"
   end if
   do j=1,npoints
      read(29,*,iostat=readstat) fdum,dos_partial(i,j,:)
      dos_partial(i,j,:)=abs(dos_partial(i,j,:))
      if (readstat .ne. 0) then
         stop "Missing data in DOSCAR? Did you add the LORBIT=11 keyword?"
      end if
   end do
end do
close(29)
write(*,*) " ...  success!"

!
!     Now select the DOS data according to element,index and orbital
!
write(*,*)
write(*,*) "Sum up partial DOS according to your element/index/orbital choice."
if (ispin .eq. 1) then
   allocate(dos_part_out(npoints))
   dos_part_out=0.d0
else if (ispin .eq. 2) then
   allocate(dos_part_out_up(npoints),dos_part_out_down(npoints))
   dos_part_out_up=0.d0
   dos_part_out_down=0.d0
end if

do i=1,natoms
   act_read=.false.
   if (dos_indices) then
      do j=1,ninds_choice
         if (i .eq. ind_list(j)) then
            act_read=.true.
         end if
      end do
   else if (dos_elements) then
      if (el_all) act_read=.true.
      do j=1,nelems_choice
         if (at_names(i) .eq. el_list(j)) then
            act_read=.true.
         end if
      end do
   end if
   if (act_read) then
      do j=1,npoints
         if (ispin .eq. 1) then
            if (orb_all) then
               dos_part_out(j)=dos_part_out(j)+sum(dos_partial(i,j,:))
            else
               do k=1,9
                  if (orb_used(k)) then
                     dos_part_out(j)=dos_part_out(j)+dos_partial(i,j,k)
                  end if
               end do
            end if
         else if (ispin .eq. 2) then
            if (orb_all) then
               do k=1,9
                  dos_part_out_up(j)=dos_part_out_up(j)+dos_partial(i,j,k*2-1)
                  dos_part_out_down(j)=dos_part_out_down(j)+dos_partial(i,j,k*2)
               end do
            else
!
!     The upspin component
!
               do k=1,9
                  if (orb_used(k)) then
                     dos_part_out_up(j)=dos_part_out_up(j)+dos_partial(i,j,k*2-1)
                  end if
               end do
!
!     The downspin component
!
               do k=1,9
                  if (orb_used(k)) then
                     dos_part_out_down(j)=dos_part_out_down(j)+dos_partial(i,j,k*2)
                  end if
               end do
            end if
         end if
      end do
   end if
end do

write(*,*) " ...  success!"
if (ispin .eq. 1) then
   if (sum(dos_part_out) .lt. 0.0001d0) then
      write(*,*)
      write(*,*) "The partial DOS is zero! Maybe no valid element or index chosen?"
      write(*,*)
   end if
else if (ispin .eq. 2) then
   if (sum(dos_part_out_up) .lt. 0.0001d0) then
      write(*,*)
      write(*,*) "The partial DOS is zero! Maybe no valid element or index chosen?"
      write(*,*)
   end if
end if
!
!     Write the partial DOS to a file
!
open(unit=30,file="dos_partial.dat",status="replace")
write(30,*) "# DOS calculated by partial_dos "
write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"

if (ispin .eq. 1) then
   write(30,*) "#   energy(eV)             DOS  "
else if (ispin .eq. 2) then
   write(30,*) "#   energy(eV)             DOS(up)            DOS(down) "
end if
do i=1,npoints
   if (ispin .eq. 1) then
      write(30,*) e_dos(i)-efermi,dos_part_out(i)
   else if (ispin .eq. 2) then
      write(30,*) e_dos(i)-efermi,dos_part_out_up(i),dos_part_out_down(i)
   end if
end do
close(30)
write(*,*) "Partial DOS of the system written to file 'dos_partial.dat'."

!
!    If spin polarization is present, write an averaged DOS
!     (e.g., for effective calculation of DOS overlaps)
!
if (ispin .eq. 2) then
   open(unit=30,file="dos_part_sum.dat",status="replace")
   write(30,*) "# DOS calculated by partial_dos (alpha and beta summed up) "
   write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"
   write(30,*) "#   energy(eV)             DOS"
   do i=1,npoints
      write(30,*) e_dos(i)-efermi,dos_part_out_up(i)+dos_part_out_down(i)
   end do
   close(30)
   write(*,*) "Summed partial DOS (alpha+beta) written to file 'dos_part_sum.dat'."
end if

22 continue

end subroutine calc_pdos

!
!    calc_cls: set up and evaluate core level shift (CLS)
!      calculations for structures with several atoms
!
subroutine calc_cls
implicit none
integer::i,j,k
character(len=2),allocatable::el_list_read(:)
character(len=2),allocatable::el_list(:)
character(len=80)::arg,adum,select_string
character(len=80)::coord_string
character(len=40)::foldername
character(len=120)::a120
character(len=120),allocatable::potcar_block(:)
integer::readstat,el_num
integer::nelems_choice,counter,block_len
logical::cls_elements,el_all
logical::mode_eval,mode_setup
logical::use_slurm
logical::coord_direct
logical::no_fermi
integer::quantum_n,quantum_l
integer,allocatable::list_active(:)
character(len=2),allocatable::el_active(:)
real(kind=8)::cell_scale,a_vec(3),b_vec(3),c_vec(3)
real(kind=8)::e_fermi
real(kind=8)::gauss_width
character(len=2),allocatable::el_init(:)
character(len=2),allocatable::el_names(:)
character(len=2),allocatable::at_names(:)
character(len=20)::spec_name
real(kind=8),allocatable::fs_val(:),is_val(:),fs_is_val(:)
real(kind=8),allocatable::fs_all(:),is_all(:),fs_is_all(:)
real(kind=8)::dft_is_ref,dft_fs_ref,exp_ref
integer,allocatable::el_numbers(:)
integer::natoms,num_active
integer::npoints
real(kind=8)::x_lo,x_hi,deltax,x_act,y_act
real(kind=8),allocatable::xyz(:,:),xyz_print(:,:)

!
!    Print general information and all possible keywords of the option
!
write(*,*)
write(*,*) "Option -cls: Management of VASP core level shift (CLS)"
write(*,*) " calculations of structures with several atoms."
write(*,*) "In the setup mode, the input files for a new calculation are"
write(*,*) " prepared, in choosing the atoms whose final state energies "
write(*,*) " shall be calculated. A POSCAR file with the structure as well"
write(*,*) " as a INCAR (for final state calc.), KPOINTS and POTCAR file"
write(*,*) " need to be located in the folder."
write(*,*) " In the INCAR file, the usual keywords for a final state calc."
write(*,*) " need to be present: CLN, CLL, CLZ, in order to specify the "
write(*,*) " orbital and the species. The species (CLN) always needs to be"
write(*,*) " 'number of the species in the system' + 1"
write(*,*) "In the evaluation mode, the prepared and finished calculations"
write(*,*) " are evaluated and the results collected such that one plot or"
write(*,*) " picture can be produced for the whole structure."
write(*,*) " The evaluation must be called in the same folder as the setup."
write(*,*) "The following command line arguments can/must be given:"
write(*,*) " -setup : chooses the setup mode for a new calculation."
write(*,*) " -eval : chooses the evaluation mode for a done calculation."
write(*,*) " -element=[element]: which element shall be analzed."
write(*,*) "   For each atom of the chosen element in the system, a "
write(*,*) "   separate final state energy calculation will be done!"
write(*,*) " -slurm : assumes that slurm is used as queue manager. A file"
write(*,*) "   named slurm_script will be copied into each folder."
write(*,*) "During the evaluation mode, Gaussian-broadened plots of the "
write(*,*) " FS and IS and FS-IS (final state effect) spectra are made."
write(*,*) " Optionally, they can be modified by the following keywords:"
write(*,*) " -plot_points=[number] : Number of plot points in the spectrum"
write(*,*) "   (default: 1000)"
write(*,*) " -gauss_width=[value]: Exponential prefactor of the Gaussian "
write(*,*) "   applied to broaden the lines of the spectrum (default 40.0)"
write(*,*) " -no_fermi : Deactivates the default correction of core level"
write(*,*) "   energies by the Fermi level of the respective calculation."


mode_setup=.false.
mode_eval=.false.
npoints=1000
gauss_width=40.d0
!
!     Determine calculation mode
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-setup") then
      mode_setup = .true.
      write(*,*)
      write(*,*) "Mode A was chosen, the setup will be done."
      write(*,*)
   end if
end do

do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-eval") then
      if (mode_setup) then
         write(*,*)
         write(*,*) "Please choose either Mode -eval or -setup, not both!"
         write(*,*)
         stop
      else
         mode_eval = .true.
         write(*,*)
         write(*,*) "Mode B was chosen, the evaluation will be done."
         write(*,*)
      end if
   end if
end do

!
!    Abort if no mode has been chosen
!
if ((.not. mode_setup) .and. (.not. mode_eval)) then
   write(*,*)
   write(*,*) "Please choose either Mode -setup or Mode -eval!"
   write(*,*)
   stop
end if


!
!    The elements whose core level energies shall be calculated
!
nelems_choice=0
cls_elements=.false.
el_all=.false.
allocate(el_list_read(10))
el_list_read="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-element=") then
      read(arg(10:),*,iostat=readstat) el_list_read
      cls_elements=.true.
   end if
end do

if (cls_elements) then
   nelems_choice=0
   do i=1,10
      if (el_list_read(i) .eq. "XX") exit
!
!    If all elements are chosen (all)
!
      if (el_list_read(i) .eq. "al") then
         el_all = .true.
         exit
      end if
      nelems_choice=nelems_choice+1
   end do

   allocate(el_list(nelems_choice))
   do i=1,nelems_choice
      el_list(i)=el_list_read(i)
   end do
   if (el_all) then
      write(*,*) "All atoms in the system were chosen for the DOS evalulation."
   else
      write(*,*) "The following elements were chosen for the DOS evaluation:"
      do i=1,nelems_choice
         write(*,*) " - ",el_list(i)
      end do
   end if
end if
if (nelems_choice .lt. 1 .and. mode_setup) then
   write(*,*) "Please give one element with the -element command!"
   stop
end if

use_slurm=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-slurm") then
      write(*,*) "The calculations will be started with slurm_scripts."
      use_slurm=.true.
   end if
end do

npoints=1000
gauss_width=40.0
if (mode_eval) then
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:13))  .eq. "-plot_points=") then
         read(arg(14:),*,iostat=readstat) npoints
         if (readstat .ne. 0) then
            write(*,*)
            stop "Check the command -plot_points=..., something went wrong!"
            write(*,*)
         end if
      end if
   end do
   if (npoints .ne. 1000) then
      write(*,*) " - The number of plot points has been set to ",npoints
   end if


   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:13))  .eq. "-gauss_width=") then
         read(arg(14:),*,iostat=readstat) gauss_width
         if (readstat .ne. 0) then
            write(*,*)
            stop "Check the command -gauss_width=..., something went wrong!"
            write(*,*)
         end if
      end if
   end do
   if ((gauss_width - 40d0) .gt. 0.0001d0) then
      write(*,*) " - The width of line broadening Gaussians has been set to ",gauss_width
   end if
end if
!
!    If the correction by Fermi level shall be deactivated
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-no_fermi") then
      no_fermi = .true.
      write(*,*)
      write(*,*) "The core level energies will not be corrected for the Fermi level!"
      write(*,*)
   end if
end do
!
!    Open the POSCAR file and read in its content
!    This is needed for both modes!
!
open(unit=56,file="POSCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The POSCAR file could not been found!"
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
!    For direct coordinates: convert coordinates to cartesians
!
if (coord_direct) then
   do i=1,natoms
      xyz_print(:,i)=xyz(1,i)*a_vec(:)+xyz(2,i)*b_vec(:)+xyz(3,i)*c_vec(:)
   end do
end if

!
!    Define individual elements for each atom
!

allocate(at_names(natoms))
counter = 1
do i=1,el_num
   do j=1,el_numbers(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do



!
!    MODE A (setup)
!

if (mode_setup) then
   allocate(list_active(natoms))
   allocate(el_active(natoms))
!
!    Determine the number and indices of the atoms that
!    shall be calculated
!
   counter=1
   do i=1,nelems_choice
      do j=1,natoms
         if (at_names(j) .eq. el_list(i)) then
            list_active(counter)=j
            el_active(counter)=el_list(i)
            counter=counter+1
         end if
      end do
   end do
   num_active=counter-1

   open(unit=57,file="active_list.dat",status="replace")
   write(57,*) "# This file contains all atoms chosen for final state CLS calculations"
   do i=1,num_active
      write(57,*) list_active(i),el_active(i)
   end do
   write(*,*) "File active_list.dat written. (used for the eval mode!)"
!
!    Setup the IS calculation folder
!    One calculation is sufficient for the whole system!
!
   call system("mkdir IS")
   call chdir("IS")
   call system("cp ../KPOINTS .")
   call system("cp ../INCAR .")
   call system("cp ../POTCAR .")
   call system("cp ../POSCAR .")
   if (use_slurm) then
      call system("cp ../slurm_script .")
   end if
!
!    Order a initial state (ICORELEVEL = 1) calculation
!
   call system("sed -i '/ICORELEVEL/c\ICORELEVEL = 1' INCAR")

   write(*,*)
   write(*,*) "Global initial state calculation prepared in folder IS"
!
!     Start the calculation if the slurm_script usage is activated
!
   if (use_slurm) then
      call system("sbatch slurm_script")
      call system("sleep 0.5")
      write(*,*) "IS calculation has started!"
   end if

   call chdir("..")



!
!    Now, setup a FS calculation folder for each selected atom!
!    In the respective POSCAR, the atom is printed in the last
!    line
!    Further, a new species is added to the POTCAR file, resembling
!    the element of the atom in the last line
!

   write(*,*)
   write(*,*) "Final state calculations for all active atoms of the "
   write(*,*) " chosen element prepared in folders:"

   allocate(potcar_block(8000))
   do i=1,num_active
      write(foldername,*) list_active(i)
      call system("mkdir " // adjustl(trim(foldername)))
      write(*,'(a,i4,a,a)') "   calculation ",i,":     ", trim(foldername)
      call chdir(adjustl(trim(foldername)))
      call system("cp ../KPOINTS .")
      call system("cp ../INCAR .")
      call system("cp ../POTCAR .")
      if (use_slurm) then
         call system("cp ../slurm_script .")
      end if
      open(unit=45,file="POSCAR",status="replace")
      write(45,*) "Final state calculation for atom ",list_active(i)
      write(45,*) cell_scale
      write(45,*) a_vec(:)
      write(45,*) b_vec(:)
      write(45,*) c_vec(:)
      do j=1,el_num
         write(45,'(a3,a2)',advance="no") "   ",el_names(j)
      end do
      write(45,'(a3,a2)',advance="no") "   ",el_active(i)
      write(45,*)
      do j=1,el_num
         if (el_names(j) .eq. el_active(i)) then
            write(45,'(a3,i8)',advance="no") "   ",el_numbers(j)-1
         else
            write(45,'(a3,i8)',advance="no") "   ",el_numbers(j)
         end if
      end do
      write(45,'(a3,i8)',advance="no") "   ",1
      write(45,*)
      write(45,*) "Direct"
!
!     Write the coordinate section, move the line of the active atom
!     to the last line
!
      do j=1,natoms
         if (j .ne. list_active(i)) then
            write(45,*) xyz(:,j)
         end if
      end do
      write(45,*) xyz(:,list_active(i))
      close(45)
!
!     Write the entry for the POTCAR file at the last line
!
      block_len=0
      open(unit=67,file="POTCAR",status="old")
      do
         read(67,'(a)',iostat=readstat) a120
         if (readstat .ne. 0) exit
         if (index(a120,'PAW_PBE '//trim(el_active(i))) .ne. 0) then

            counter=1
            potcar_block(counter)=a120
            do
               counter=counter+1
               read(67,'(a120)') a120
               potcar_block(counter)=a120
               if (index(a120,'End of Dataset') .ne. 0) exit
            end do
            block_len=counter
         end if
         if (block_len .ne. 0) exit
      end do
      close(67)
      open(unit=68,file="POTCAR",status="old",access="append")
      do j=1,block_len
         write(68,'(a120)') potcar_block(j)
      end do
      close(68)
!
!     Start the calculation if the slurm_script usage is activated
!
      if (use_slurm) then
         call system("sbatch slurm_script")
         call system("sleep 0.5")
      end if

      call chdir("..")

   end do
   if (use_slurm) then
      write(*,*) "Final state calculations started!"
   else
      write(*,*) "Please go into all generated folders and start the calculations"
      write(*,*) " manually, or add the -slurm option and a slurm_script in the "
      write(*,*) " main folder."
   end if

   write(*,*)
   write(*,*) "Wait until all calculations are finished, then, run this program"
   write(*,*) " again, now in the -eval mode."
end if

!
!    MODE B: (evaluation)
!

if (mode_eval) then
!
!     First, read in all atoms that shall be considered
!

   open(unit=57,file="active_list.dat",status="old")
   counter=-1
   do
      read(57,*,iostat=readstat)
      if (readstat .ne. 0) exit
      counter=counter+1
   end do
   close(57)

   if (counter .le. 0) then
      write(*,*) "The file 'active_list.dat' contains to atom indices!"
      stop
   end if

   num_active=counter
   allocate(list_active(num_active))
   allocate(el_active(num_active))
!
!     The array with the initial and final state energies
!
   allocate(fs_val(num_active))
   allocate(is_val(num_active))
   allocate(fs_is_val(num_active))
   allocate(fs_all(natoms))
   allocate(is_all(natoms))
   allocate(fs_is_all(natoms))


   write(*,*) "List of atoms that were calculated:"
   open(unit=57,file="active_list.dat",status="old")
   read(57,*)
   write(*,*) "(number    atom index     element)"
   do i=1,num_active
      read(57,*) list_active(i),el_active(i)
      write(*,'(a,i5,a,i5,a,a,a)') "  -",i,":     ",list_active(i),"    &
     &       (",trim(el_active(i)),")"
   end do
   close(57)
!
!     Read in the FS results for each atom, separately
!     Loop through folders of calculations
!
   write(*,*)
   do i=1,num_active
      write(foldername,*) list_active(i)
      call chdir(adjustl(trim(foldername)))
      write(*,'(a,a,a)') "  Evaluating folder ",trim(foldername)," ..."
!
!     For first folder: open INCAR and determine orbital to be analyzed
!
      if (i .eq. 1) then
         quantum_n=0
         quantum_l=0
         open(unit=58,file="INCAR",status="old")
         do
            read(58,'(a)',iostat=readstat) a120
            if (readstat .ne. 0) exit
            if (index(a120,"CLN") .ne. 0) then
               read(a120,*) adum,adum,quantum_n
            else if (index(a120,"CLL") .ne. 0) then
               read(a120,*) adum,adum,quantum_l
            end if
         end do
         close(58)
         if (quantum_n .lt. 1) then
            write(*,*) "The n quantum number has not been given in the INCAR file!"
            stop
         end if
         if (quantum_l .lt. 1) then
            write(*,*) "The l quantum number has not been given in the INCAR file!"
            stop
         end if
         write(*,*) "The chosen quantum numbers are: "
         write(*,'(a,i1)') "  - N = ",quantum_n
         if (quantum_l .eq. 0) then
            write(*,'(a,i1,a)') "  - L = ",quantum_l," (s)"
         else if (quantum_l .eq. 1) then
            write(*,'(a,i1,a)') "  - L = ",quantum_l," (p)"
         else if (quantum_l .eq. 2) then
            write(*,'(a,i1,a)') "  - L = ",quantum_l," (d)"
         else if (quantum_l .eq. 3) then
            write(*,'(a,i1,a)') "  - L = ",quantum_l," (f)"
         end if
      end if
!
!     Open OUTCAR file and read in core level energies for current index
!     It is assumed that the "active" atom is always the last one!
!
      write(*,*) natoms
      if (natoms .lt. 10) then
         write(spec_name,'(i1,a)') natoms,"-"
      else if (natoms .lt. 100) then
         write(spec_name,'(i2,a)') natoms,"-"
      else if (natoms .lt. 1000) then
         write(spec_name,'(i3,a)') natoms,"-"
      end if
      open(unit=47,file="OUTCAR",status="old")
      do
         read(47,'(a)',iostat=readstat) a120
         if (readstat .ne. 0) exit
         if (index(a120," "//trim(spec_name)//" ") .ne. 0) then
!
!     Depending on the N and L quantum numbers of the core level, determine
!      the position of the number to be read in
!
!     The 1s orbital
            if (quantum_n .eq. 1 .and. quantum_l .eq. 0) then
               read(a120,*) adum,adum,fs_val(i)
!     The 2s orbital
            else if (quantum_n .eq. 2 .and. quantum_l .eq. 0) then
               read(a120,*) adum,adum,adum,adum,fs_val(i)
!     The 2p orbital
            else if (quantum_n .eq. 2 .and. quantum_l .eq. 1) then
               read(a120,*) adum,adum,adum,adum,adum,adum,fs_val(i)
!     The 3s orbital
            else if (quantum_n .eq. 3 .and. quantum_l .eq. 0) then
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,adum,fs_val(i)
!     The 3p orbital
            else if (quantum_n .eq. 3 .and. quantum_l .eq. 1) then
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,fs_val(i)
!     The 3d orbital
            else if (quantum_n .eq. 3 .and. quantum_l .eq. 2) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,fs_val(i)
!     The 4s orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 0) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,fs_val(i)
!     The 4p orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 1) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,fs_val(i)
!     The 4d orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 2) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,fs_val(i)
!     The 4f orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 3) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,fs_val(i)
!     The 5s orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 0) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,fs_val(i)
!     The 5p orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 1) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,fs_val(i)
!     The 5d orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 2) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,fs_val(i)
!     The 5f orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 3) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,fs_val(i)

            end if
         end if
!
!     Read in the Fermi energy for the correct shift in energy
!
         if (index(a120,"Fermi energy:") .ne. 0) then
            read(a120,*) adum,adum,e_fermi
         end if
      end do
      close(47)
!
!     Calculate the current final state energy, in the same scale as experiment
!

      if (no_fermi) then
         fs_val(i)=-fs_val(i)
      else
         fs_val(i)=-(fs_val(i)-e_fermi)
      end if

      call chdir("..")
   end do

!
!     Now read in the IS results of all atoms, for the same orbitals as the
!     FS calculations
!
   call chdir("IS")
   do i=1,num_active
      if (list_active(i) .lt. 10) then
         write(spec_name,'(i1,a)') list_active(i),"-"
      else if (list_active(i) .lt. 100) then
         write(spec_name,'(i2,a)') list_active(i),"-"
      else if (list_active(i) .lt. 1000) then
         write(spec_name,'(i3,a)') list_active(i),"-"
      end if
      open(unit=47,file="OUTCAR",status="old")
      do
         read(47,'(a)',iostat=readstat) a120
         if (readstat .ne. 0) exit
         if (index(a120," "//trim(spec_name)//" ") .ne. 0) then
!
!     Depending on the N and L quantum numbers of the core level, determine
!      the position of the number to be read in
!
!     The 1s orbital
            if (quantum_n .eq. 1 .and. quantum_l .eq. 0) then
               read(a120,*) adum,adum,is_val(i)
!     The 2s orbital
            else if (quantum_n .eq. 2 .and. quantum_l .eq. 0) then
               read(a120,*) adum,adum,adum,adum,is_val(i)
!     The 2p orbital
            else if (quantum_n .eq. 2 .and. quantum_l .eq. 1) then
               read(a120,*) adum,adum,adum,adum,adum,adum,is_val(i)
!     The 3s orbital
            else if (quantum_n .eq. 3 .and. quantum_l .eq. 0) then
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,adum,is_val(i)
!     The 3p orbital
            else if (quantum_n .eq. 3 .and. quantum_l .eq. 1) then
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,is_val(i)
!     The 3d orbital
            else if (quantum_n .eq. 3 .and. quantum_l .eq. 2) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,is_val(i)
!     The 4s orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 0) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,is_val(i)
!     The 4p orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 1) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,is_val(i)
!     The 4d orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 2) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,is_val(i)
!     The 4f orbital
            else if (quantum_n .eq. 4 .and. quantum_l .eq. 3) then
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,is_val(i)
!     The 5s orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 0) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,is_val(i)
!     The 5p orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 1) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,is_val(i)
!     The 5d orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 2) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,is_val(i)
!     The 5f orbital
            else if (quantum_n .eq. 5 .and. quantum_l .eq. 3) then
               read(47,'(a)',iostat=readstat) a120
               read(47,'(a)',iostat=readstat) a120
               read(a120,*) adum,adum,adum,adum,adum,adum,adum,is_val(i)

            end if
         end if
!
!     Read in the Fermi energy for the correct shift in energy
!
         if (index(a120,"Fermi energy:") .ne. 0) then
            read(a120,*) adum,adum,e_fermi
         end if
      end do
!
!     Calculate the current initial state energy, in the same scale as experiment
!
      if (no_fermi) then
         is_val(i)=-is_val(i)
      else
         is_val(i)=-(is_val(i)-e_fermi)
      end if

      close(47)
   end do

   call chdir("..")

!
!     Reference the CLS values to the experimental XPS data
!
   do i=1,num_active
!     For Pt 4f
      if (trim(el_active(i)) .eq. "Pt") then
         dft_is_ref=67.0709
         dft_fs_ref=76.6047
         exp_ref=71.0d0
!     For Ga 3d
      else if (trim(el_active(i)) .eq. "Ga") then
         dft_is_ref=14.4171
         dft_fs_ref=21.9421
         exp_ref=18.7d0
      end if

      fs_val(i)=fs_val(i)-dft_fs_ref+exp_ref
      is_val(i)=is_val(i)-dft_is_ref+exp_ref
      fs_is_val(i)=fs_val(i)-is_val(i)
   end do

!
!     Produce Gaussian-broadened spectrum of core levels for comparison with experiment
!
!     The final state:
!
!     First, determine plot-limits
!
   x_lo=real(floor(minval(fs_val)))-1.d0
   x_hi=real(ceiling(maxval(fs_val)))+1.d0
   deltax=(x_hi-x_lo)/real(npoints)

   open(unit=57,file="plot_fs.dat",status="replace")
   write(57,*) "# This file contains a final state CLS spectrum for the current system"
   write(57,*) "# Produced by the script manage_cls, part of VASP4CLINT"
   do i=1,npoints
      x_act=x_lo+(i-1)*deltax
      y_act=0
      do j=1,num_active
         y_act=y_act+exp(-(x_act-fs_val(j))**2*gauss_width)
      end do
      write(57,*) x_act,y_act
   end do
   close(57)
!
!     The initial state:
!
   x_lo=real(floor(minval(is_val)))-1.d0
   x_hi=real(ceiling(maxval(is_val)))+1.d0
   deltax=(x_hi-x_lo)/real(npoints)

   open(unit=57,file="plot_is.dat",status="replace")
   write(57,*) "# This file contains a initial state CLS spectrum for the current system"
   write(57,*) "# Produced by the script manage_cls, part of VASP4CLINT"
   do i=1,npoints
      x_act=x_lo+(i-1)*deltax
      y_act=0
      do j=1,num_active
         y_act=y_act+exp(-(x_act-is_val(j))**2*gauss_width)
      end do
      write(57,*) x_act,y_act
   end do
   close(57)

!
!     The final state effect (FS-IS):
!
   x_lo=real(floor(minval(fs_is_val)))-1.d0
   x_hi=real(ceiling(maxval(fs_is_val)))+1.d0
   deltax=(x_hi-x_lo)/real(npoints)

   open(unit=57,file="plot_fs-is.dat",status="replace")
   write(57,*) "# This file contains a final state effect (FS-IS) spectrum for the current system"
   write(57,*) "# Produced by the script manage_cls, part of VASP4CLINT"
   do i=1,npoints
      x_act=x_lo+(i-1)*deltax
      y_act=0
      do j=1,num_active
         y_act=y_act+exp(-(x_act-fs_is_val(j))**2*gauss_width)
      end do
      write(57,*) x_act,y_act
   end do
   close(57)

   write(*,*)
   write(*,*) "Files for Gaussian-broadened spectra of CLS values written:"
   write(*,*) " - plot_fs.dat : Final state values "
   write(*,*) " - plot_is.dat : Initial state values "
   write(*,*) " - plot_fs-is.dat : Final state effect values "
!
!     Fill arrays with CLS values, spanning all atoms of the system
!
   fs_all=0.d0
   is_all=0.d0
   fs_is_all=0.d0
   do i=1,num_active
      fs_all(list_active(i))=fs_val(i)
      is_all(list_active(i))=is_val(i)
      fs_is_all(list_active(i))=fs_is_val(i)
   end do

!
!     Print out PDB file for coloring of atoms by their CLS
!     (can be visualized with VMD)
!
!     The FS values
!
   open(unit=20,file="show_fs.pdb")
   write(20,'(a)') "COMPND    FINAL HEAT OF FORMATION =     0.000000"
   write(20,'(a)') "AUTHOR    GENERATED BY PROGRAM MANAGE CLS"
   do i=1,natoms
      write(20,'(a,i5,a,a,a,3f8.3,f7.3,f5.2,a,a)') "HETATM",i," ",at_names(i), &
                 & "   UNL     1    ",xyz_print(:,i),fs_all(i),0d0,"          ",at_names(i)
   end do
   write(20,*) "MASTER        0    0    0    0    0    0    0    0  180    0  180    0"
   write(20,*) "END"
   close(20)
!
!     The IS values
!
   open(unit=20,file="show_is.pdb")
   write(20,'(a)') "COMPND    FINAL HEAT OF FORMATION =     0.000000"
   write(20,'(a)') "AUTHOR    GENERATED BY PROGRAM MANAGE CLS"
   do i=1,natoms
      write(20,'(a,i5,a,a,a,3f8.3,f7.3,f5.2,a,a)') "HETATM",i," ",at_names(i), &
                 & "   UNL     1    ",xyz_print(:,i),is_all(i),0d0,"          ",at_names(i)
   end do
   write(20,*) "MASTER        0    0    0    0    0    0    0    0  180    0  180    0"
   write(20,*) "END"
   close(20)
!
!     The FS-IS values (final state effect)
!
   open(unit=20,file="show_fs-is.pdb")
   write(20,'(a)') "COMPND    FINAL HEAT OF FORMATION =     0.000000"
   write(20,'(a)') "AUTHOR    GENERATED BY PROGRAM MANAGE CLS"
   do i=1,natoms
      write(20,'(a,i5,a,a,a,3f8.3,f7.3,f5.2,a,a)') "HETATM",i," ",at_names(i), &
                 & "   UNL     1    ",xyz_print(:,i),fs_is_all(i),0d0,"          ",at_names(i)
   end do
   write(20,*) "MASTER        0    0    0    0    0    0    0    0  180    0  180    0"
   write(20,*) "END"
   close(20)

   write(*,*)
   write(*,*) "Files for spatial visualization of CLS values written:"
   write(*,'(a,f10.5,a,f10.5,a)') " - show_fs.pdb : Final state values (range: ", &
             & minval(fs_val),"  to ",maxval(fs_val),")"
   write(*,'(a,f10.5,a,f10.5,a)') " - show_is.pdb : Initial state values (range: ", &
             & minval(is_val),"  to ",maxval(is_val),")"
   write(*,'(a,f10.5,a,f10.5,a)') " - show_fs-is.pdb : Final state effect values (range: ", &
             & minval(fs_is_val),"  to ",maxval(fs_is_val),")"
   write(*,*) " Open these files with VMD and select 'coloring method: occupancy'"

end if


end subroutine calc_cls

