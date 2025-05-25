!
!    manage_cls: set up and evaluate core level shift (CLS)
!      calculations for structures with several atoms
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!
program manage_cls
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
write(*,*) "PROGRAM manage_cls: Management of VASP core level shift (CLS)"
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
write(*,*) " -overview:  print an overview of all scripts and programs "
write(*,*) "   in utils4VASP"
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

write(*,*)
write(*,*) "manage_cls terminated normally."
write(*,*)


end program manage_cls        
