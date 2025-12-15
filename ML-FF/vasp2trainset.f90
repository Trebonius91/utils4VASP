!
!    vasp2trainset: This program takes the output of either a nudged 
!      elastic band (NEB) calculation, a VASP on-the-fly machine
!      learning force field (ML-FF) calculation, a VASP ML-FF production run
!      or a MACE foundation model calculation and generates input files
!      containing structures, energies and gradients that can be used to fit
!      an artificial neural network (ANN) potential with the aenet program from
!      atomisticnet on github or a MACE finetuning model for the system of interest
!      (see there for details)      
!
!    Part of utils4vasp
!
!     Julien Steffen,     2023 (julien.steffen@fau.de) - main program
!     Maximilian Bechtel, 2025 (maxi.bechtel@fau.de)   - extension to VASP ML-FF
!                                                      - production runs
!

program vasp2trainset

implicit none

integer::i,j,k,inc,inc2,inc3
logical::eval_neb,eval_mlff,does_exist,eval_md,eval_aimd
character(len=120)::arg
character(len=160)::a160
character(len=80)::a80,md_mode
character(len=50)::basename,out_name,foldername
character(len=3)::neb_paths(100)
real(kind=8)::energy
real(kind=8)::cell_vecs(3,3)
real(kind=8),allocatable::xyz(:,:)
real(kind=8),allocatable::grad(:,:)
character(len=1),allocatable::selective(:,:)
logical::struc_select
integer::at_numbers(10)
real(kind=8)::fdum
character(len=50)::adum
character(len=2),allocatable::at_names(:)
character(len=150),allocatable::xdat_content(:)
integer::readstat
integer::xdat_lines
integer::nfolders,nelems,natoms,nstrucs,nframes
integer::read_freq,frame_act,frames_read
! real(kind=8)::stress(3,3)
integer::el_nums(10)
logical::write_aenet,write_mace
logical::header
logical::skip_large
logical,allocatable::skip_frame(:)
character(len=2)::el_syms(10)
character(len=150)::cdum
character(len=80)::outcar_name
real(kind=8)::factor,scale_dum
integer::openstat,alloc_stat
character(len=100)::alloc_err
integer::traj_freq
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
      write(*,*) " - gen_poscar.py: Generate POSCAR of alloy and surface structures"
      write(*,*) " - analyze_poscar.py: Analyze POSCAR, generate KPOINTS and POTCAR"
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
write(*,*) " *** utils4VASP -- utility scripts and programs for VASP ***"
write(*,*) "PROGRAM vasp2trainset: Transforms the output of a nudged elastic band (NEB)"
write(*,*) " or an on-the-fly VASP ML-FF generation trajectory to the training set"
write(*,*) " format needed for either (a) the aenet program from atomisticnet for"
write(*,*) " an artificial neural network (ANN) potential or (b) the MACE suite to"
write(*,*) " fit one of the MACE foundation models for the system of interest."
write(*,*) "For NEB calculations, start the program in the main folder where "
write(*,*) " the frames are stored in subfolders 00,01,02, ..."
write(*,*) "For ML-FF calculations, start the program where the ML_AB file generated"
write(*,*) " by the on-the-fly learning of interest is located."
write(*,*) "The command line argument decides which kind of calculation will be done:"
write(*,*) " -neb : evaluates a NEB calculation."
write(*,*) " -ml_ff : evaluates a ML-FF learning calculation."
write(*,*) " -aimd : evaluates a VASP AIMD calculation (OUTCAR file)"
write(*,*) " -outcar_name=[name] : Name of the OUTCAR file (default: OUTCAR)"
write(*,*) " -read_freq : Only each Nth AIMD frame will be converted to training set."
write(*,*) " -md_traj=[mode] : picks frames from a VASP trajectory (XDATCAR)"
write(*,*) " Ensemble (NVT or NpT) of the trajectory is detected automatically."
write(*,*) "   If mode = setup, a folder with a POSCAR file will be generated for each"
write(*,*) "   XDATCAR frame, a VASP single point calculation needs to be done in each"
write(*,*) "   of the folders (IBRION=-1,NSW=0)."
write(*,*) "   If mode = eval, the folders with the VASP calculations will be evaluated"
write(*,*) "   and the training set file(s) will be written with the energy/gradients."  
write(*,*) " -traj_freq=[number] : For the md_traj (mode setup), only each Nth trajectory"
write(*,*) "   frame will be written to file for VASP calculation."
write(*,*) " -name=[speficier] : give a unique string to specify the current calculation"
write(*,*) " -aenet : The output will be for the aenet program (atomisticnet)"
write(*,*) " -mace : The output will be for the MACE program (fit foundation model)"
write(*,*) " -skip_large : If structures with large gradient components shall be skipped"
write(*,*) "   during the collection of the results."
write(*,*) "For the aenet option, after execution, a number of XSF files (one for each"
write(*,*) " structure) is written into a new folder named 'xsf_files'. Further, a list"
write(*,*) " of all file names is written into the file 'xsf_list.dat'. The file names "
write(*,*) " begin with the basename defined with -name. For both NEB and ML-FF training"
write(*,*) " sets, the selective dynamics of POSCAR files can be read in, such that in "
write(*,*) " the evaluation with mlp_quality, only activated degrees of freedom are "
write(*,*) " considered (else, gradient errors can occur, since aenet only trains the"
write(*,*) " energies, and no information about fixed atoms will be gathered."
write(*,*) "For the mace option, a single file named [base].xyz containing the training"
write(*,*) " set files is written, where base is defined by the -name keyword. Further,"
write(*,*) " header lines are added where the atomic energies can be added later."
write(*,*) "For an overview of utils4VASP, give the -overview command."

!
!     If a NEB shall be evaluated
!
eval_neb = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-neb") then
      eval_neb = .true.
   end if
end do

!
!     If a ML-FF (ML_AB) shall be evaluated
!
eval_mlff = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-ml_ff") then
      eval_mlff = .true.
   end if
end do

!
!     If a AIMD trajectory (OUTCAR) shall be evaluated
!
eval_aimd = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-aimd") then
      eval_aimd = .true.
   end if
end do
!
!     The frequency in which the AIMD frames will be written 
!     to the training set file
!
read_freq=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-read_freq=") then
      read(arg(12:),*) read_freq
   end if
end do

!
!     The name of the OUTCAR file, if not the default (OUTCAR)
!
outcar_name="OUTCAR"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-outcar_name=") then
      read(arg(14:),*) outcar_name
   end if
end do

!
!     If a MD trajectory shall be evaluated
!
md_mode="xxxx"
eval_md = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-md_traj=") then
      eval_md = .true.
      read(arg(10:),*) md_mode
   end if
end do
!
!     The frequency in which the MD frames will be written 
!     out to VASP single point calculation input folders
!
traj_freq=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-traj_freq=") then
      read(arg(12:),*) traj_freq
   end if
end do

!
!     Which basename the written XSF files shall have
!
basename="xxxxx"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-name=") then
      read(arg(7:),*) basename
   end if
end do

!
!     If the training set in the aenet format shall be written
!
write_aenet = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-aenet") then
      write_aenet = .true.
   end if
end do

!
!     If the training set in the MACE format shall be written
!
write_mace = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-mace") then
      write_mace = .true.
   end if
end do

!
!     If structures with large gradient components shall be skipped 
!     during the collection of the results
!
skip_large = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-skip_large") then
      skip_large = .true.
   end if
end do



if (.not. eval_neb .and. .not. eval_mlff .and. .not. eval_aimd .and. .not. eval_md) then
   write(*,*)
   write(*,*) "Please choose either the -neb or the -ml_ff or the -aimd or the -md_traj option!"
   write(*,*)
   stop
end if
if (eval_neb .and. eval_mlff .or. eval_neb .and. eval_md .or. eval_mlff .and. eval_md) then
   write(*,*)
   write(*,*) "Please choose either the -neb or the -ml_ff or the -md_traj option!"
   write(*,*)
   stop
end if
if (eval_md) then
   if (trim(md_mode) .ne. "setup" .and. trim(md_mode) .ne. "eval") then
      write(*,*)
      write(*,*) "You have chosen the -eval_md option, but no valid mode (setup or eval)"
      write(*,*)
      stop
   end if
end if
if (trim(md_mode) .ne. "setup") then
   if (.not. write_aenet .and. .not. write_mace) then
      write(*,*)
      write(*,*) "Please choose either the -aenet or -mace output option!"
      write(*,*)
      stop
   end if 
end if
if (write_aenet .and. write_mace) then
   write(*,*)
   write(*,*) "Please choose either the -aenet or -mace output option!"
   write(*,*)
   stop
end if
if (trim(md_mode) .ne. "setup") then
   if (trim(basename) .eq. "xxxxx") then
      if (write_aenet) then
         write(*,*)
         write(*,*) "Please give a basename for the XSF files with the -name keyword!"
         write(*,*)
         stop
      end if
      if (write_mace) then
         write(*,*)
         write(*,*) "Please give a basename for the trainset file with the -name keyword!"
         write(*,*)
         stop
      end if
   end if
end if         
!
!     Generate the folder containing the XFS files with the output
!     information (for aenet)
!
if (write_aenet) then
   inquire(file="xsf_files/",exist=does_exist)
   if (.not. does_exist) then
      call system("mkdir xsf_files")
   end if
!
!     Open the file listing the names of all XFS files (for aenet)
!
   open(unit=45,file="xsf_list.dat",status="replace")
   write(45,*) "# This file contains the names of all XSF files, in the formate"
   write(45,*) "# needed to setup the aenet calculation."
   write(45,*) "# Copy the following lines into the ... input file."
   write(45,*)
end if
!
!     Open the training set file (for MACE)
!
if (write_mace) then
   open(unit=48,file=trim(basename)//".xyz")   
end if
!
!     Evaluate a NEB calculation
!
if (eval_neb) then
!
!     Determine the number of NEB frames
!   
   i=1
   write(neb_paths(1),'(i1,i1,a)') 0,0,"/"
   do
      if (i .lt. 10) then
         write(neb_paths(i+1),'(i1,i1,a)') 0,i,"/"    
      else     
         write(neb_paths(i+1),'(i2,a)') i,"/"
      end if  
      inquire(file=neb_paths(i+1), exist=does_exist)
      if (.not. does_exist) exit
      i=i+1
   end do
   nfolders=i
   write(*,*) "There are ",nfolders-1," NEB frames present."
!
!     Now go into all folders and open the OUTCAR files
!     For frame 00 and N+1, the OUTCARs are viewed as optional
!
   do i=1,nfolders
      el_nums=0
      inc3=1
      call chdir(neb_paths(i))
      inquire(file="OUTCAR", exist=does_exist)
      if (.not. does_exist) then
         if (i .eq. 1 .or. i .eq. nfolders) then
            write(*,*) "Folder ",neb_paths(i)," with no OUTCAR skipped."
            call chdir("..")
            cycle
         else 
            write(*,*) "No OUTCAR file could be found in the folder",neb_paths(i),"!"
            stop
         end if        
      end if
!
!     Read in the relevant informations and write them into the XSF folder 
!
      open(unit=11,file="OUTCAR",status="old")
      do 
         read(11,'(a)',iostat=readstat) a160
         if (readstat .ne. 0) exit
!
!     Number of atoms and element symbols of them
!
         if (index(a160,"ions per type =") .ne. 0) then
            read(a160,*,iostat=readstat) a80,a80,a80,a80,el_nums
            do j=1,10
               if (el_nums(j) .ne. 0) then
                  nelems=j
               end if     
            end do
            natoms=sum(el_nums(1:nelems))
            if (allocated(xyz)) deallocate(xyz)
            allocate(xyz(3,natoms))
            if (allocated(grad)) deallocate(grad)
            allocate(grad(3,natoms))
            if (allocated(at_names)) deallocate(at_names)
            allocate(at_names(natoms))
!
!     Determine element symbols for all atoms
!        
            at_names="XX"
            inc2=1
            do j=1,nelems
               do k=1,el_nums(j)
                  at_names(inc2)=el_syms(j)
                  inc2=inc2+1
               end do
            end do

         end if   
         if (index(a160,"TITEL") .ne. 0) then
            read(a160,*) a80,a80,a80,el_syms(inc3)
            inc3=inc3+1
         end if        
!
!     Read energy, basis vectors, positions and gradients for current frame
!

         if (index(a160,'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)').ne.0) then
            read(11,*)
            read(11,*)
            read(11,*)
            read(11,*) a80,a80,a80,a80,a80,a80,energy
         end if
         if(index(a160,'VOLUME and BASIS-vectors').ne.0)then
            read(11,*)
            read(11,*)
            read(11,*)
            read(11,*)
            read(11,*) cell_vecs(1,:)
            read(11,*) cell_vecs(2,:)
            read(11,*) cell_vecs(3,:)
         end if
         if (index(a160,'POSITION                                       TOTAL-FORCE (eV/Angst)') &
                         & .ne. 0) then
            read(11,*)
            do j=1,natoms
               read(11,*) xyz(:,j),grad(:,j)
            end do
         end if        

      end do
      close(11)

!
!     Read the POSCAR file for information of selective dynamics
!
      if (allocated(selective)) deallocate(selective)
      allocate(selective(3,natoms))
      selective="T"
      struc_select=.true.
      open(unit=12,file="POSCAR",status="old",iostat=readstat)
      if (readstat .ne. 0) then
         write(*,*)
         write(*,*) "No POSCAR found, all atoms are assumed to be movable."
         write(*,*)
         struc_select=.false.
      else
         do j=1,7
            read(12,*)
         end do 
         read(12,*) a160
         if ((index(a160,'Selective').ne.0) .or. (index(a160,'selective').ne.0)) then
            read(12,*)
            do j=1,natoms
               read(12,*) fdum,fdum,fdum,selective(:,j)
            end do 
         end if        
      end if         


      call chdir("..")
!
!     Write output for current frame (aenet)
!
      if (write_aenet) then
         if (i-1 .lt. 10) then
            write(out_name,'(a,a,i1)') trim(basename),"_frame",i-1
         else
            write(out_name,'(a,a,i2)') trim(basename),"_frame",i-1
         end if
         open(unit=30,file="xsf_files/" // trim(out_name) // ".xsf")
         write(30,'(a,f14.8,a)') "# total energy = ",energy," eV"
         write(30,*)
         write(30,'(a)') "CRYSTAL"
         write(30,'(a)') "PRIMVEC"
         do j=1,3
            write(30,'(3f14.8)') cell_vecs(j,:)
         end do
         write(30,'(a)') "PRIMCOORD"
         write(30,'(i6,i1)') natoms,1
         do j=1,natoms
            write(30,'(a,a,6f14.8,6a)') at_names(j)," ",xyz(:,j),grad(:,j)," ",selective(1,j), &
                         & " ",selective(2,j)," ",selective(3,j)
         end do

         close(30)
!
!     Write name of output file to list of filenames
!
         write(45,'(a,a,a)') "./xsf_files/", trim(out_name),".xsf"
      end if
!
!     Write output for current frame (MACE)
!
      if (write_mace) then
         write(48,*) natoms
         write(48,'(a,9f13.8,a,a,f15.8,a)') 'Lattice="',cell_vecs(1,:),cell_vecs(2,:),& 
                 & cell_vecs(3,:),'" Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3',&
                 & ' Nmols=1 REF_energy=',energy,' pbc="T T T"'
         do j=1,natoms
            write(48,'(a,a,3f14.8,a,3f14.8)') at_names(j),"  ",xyz(:,j),"  0  ",grad(:,j)
         end do
      end if
   end do

end if        

!
!     Evaluate a ML-FF calculation (ML_AB file)
!
if (eval_mlff) then


!
!     Read the POSCAR file for information of selective dynamics
!     If it is not there, assume that all atoms are movable
!  
   open(unit=12,file="POSCAR",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*)      
      write(*,*) "No POSCAR found, all atoms are assumed to be movable."
      write(*,*)
      struc_select=.false.
   else
      do j=1,6
         read(12,*)
      end do
      at_numbers=0
      read(12,'(a)') a160
      read(a160,*,iostat=readstat) at_numbers
      natoms=sum(at_numbers)
      allocate(selective(3,natoms))
      selective="T"
      struc_select=.true.

      read(12,*) a160
      if ((index(a160,'Selective').ne.0) .or. (index(a160,'selective').ne.0)) then
         read(12,*)
         do j=1,natoms
            read(12,*) fdum,fdum,fdum,selective(:,j)
         end do
      end if
   end if
   
!
!     Try to open the ML_AB file, abort if its not there
!
   inquire(file="ML_AB", exist=does_exist)
   if (.not. does_exist) then
      write(*,*) "The file ML_AB is not there!"
      stop
   end if         
   open(unit=27,file="ML_AB",status="old")
   read(27,*)
   read(27,*)
   read(27,*)
   read(27,*)
   read(27,*) nstrucs
   write(*,*) "Number of structures in the ML_AB fille: ",nstrucs
   inc=0
   do
      read(27,'(a)',iostat=readstat) a160
      if (readstat .ne. 0) exit
!
!     Read in all information for the current configuration
!
      if (index(a160,"Configuration num.") .ne. 0) then
         do i=1,7
            read(27,*)
         end do
         read(27,*) nelems
         do i=1,3
            read(27,*) 
         end do
         read(27,*) natoms
         if (.not. struc_select) then
            if (allocated(selective)) deallocate(selective)
            allocate(selective(3,natoms))
            selective="T"
         end if        
         do i=1,3
            read(27,*)
         end do
         do i=1,nelems
            read(27,*) el_syms(i),el_nums(i)
         end do      
         do i=1,7
            read(27,*)
         end do
         read(27,*) cell_vecs(1,:)
         read(27,*) cell_vecs(2,:)
         read(27,*) cell_vecs(3,:)
         do i=1,3
            read(27,*)
         end do
         if (allocated(xyz)) deallocate(xyz)
         allocate(xyz(3,natoms),stat=alloc_stat,errmsg=alloc_err)
            if (alloc_stat .ne. 0) then
            write(*,*) "Allocation of xyz array failed:",trim(alloc_err)
            write(*,*) "Please read in a smaller ML_AB or increase memory!"
            stop
         end if
         if (allocated(grad)) deallocate(grad)
         allocate(grad(3,natoms),stat=alloc_stat,errmsg=alloc_err)
            if (alloc_stat .ne. 0) then
            write(*,*) "Allocation of grad array failed:",trim(alloc_err)
            write(*,*) "Please read in a smaller ML_AB or increase memory!"
            stop
         end if
         if (allocated(at_names)) deallocate(at_names)
         allocate(at_names(natoms),stat=alloc_stat,errmsg=alloc_err)
            if (alloc_stat .ne. 0) then
            write(*,*) "Allocation of at_names array failed:",trim(alloc_err)
            write(*,*) "Please read in a smaller ML_AB or increase memory!"
            stop
         end if

!
!     Determine element symbols for all atoms
!
         at_names="XX"
         inc2=1
         do j=1,nelems
            do k=1,el_nums(j)
               at_names(inc2)=el_syms(j)
               inc2=inc2+1
            end do
         end do

         do i=1,natoms
            read(27,*) xyz(:,i)
         end do
         do i=1,3
            read(27,*)
         end do
         read(27,*) energy
         do i=1,3
            read(27,*)
         end do
         do i=1,natoms
            read(27,*) grad(:,i)
         end do
         inc=inc+1
!
!     Write the main output file for current structure (aenet)
!
         if (write_aenet) then
            if (inc .lt. 10) then
               write(out_name,'(a,a,i1)') trim(basename),"_struc",inc
            else if (inc .lt. 100) then 
               write(out_name,'(a,a,i2)') trim(basename),"_struc",inc
            else if (inc .lt. 1000) then
               write(out_name,'(a,a,i3)') trim(basename),"_struc",inc
            else if (inc .lt. 10000) then
               write(out_name,'(a,a,i4)') trim(basename),"_struc",inc
            else
               write(out_name,'(a,a,i5)') trim(basename),"_struc",inc
            end if        
            open(unit=30,file="xsf_files/" // trim(out_name) // ".xsf")
            write(30,'(a,f14.8,a)') "# total energy = ",energy," eV"
            write(30,*)
            write(30,'(a)') "CRYSTAL"
            write(30,'(a)') "PRIMVEC"
            do j=1,3
               write(30,'(3f14.8)') cell_vecs(j,:)
            end do
            write(30,'(a)') "PRIMCOORD"
            write(30,'(i6,i1)') natoms,1
            do j=1,natoms
               write(30,'(a,a,6f14.8,6a)') at_names(j)," ",xyz(:,j),grad(:,j)," ",selective(1,j), &
                         & " ",selective(2,j)," ",selective(3,j)
            end do
            close(30)
!
!     Write name of output file to list of filenames
!
            write(45,'(a,a,a)') "./xsf_files/", trim(out_name),".xsf"
         end if
!
!     Write the current frame to the training set file (MACE)
!
         if (write_mace) then
            write(48,*) natoms
            write(48,'(a,9f13.8,a,a,f15.8,a)') 'Lattice="',cell_vecs(1,:),cell_vecs(2,:),&
                 & cell_vecs(3,:),'" Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3',&
                 & ' Nmols=1 REF_energy=',energy,' pbc="T T T"'
            do j=1,natoms
               write(48,'(a,a,3f14.8,a,3f14.8)') at_names(j),"  ",xyz(:,j),"  0  ",grad(:,j)
            end do
         end if

      end if  
         
      if (inc .gt. nstrucs) exit
   end do

   close(27)

end if       

!
!     Evaluate an AIMD trajectory: read in positions and forces from OUTCAR file
!
if (eval_aimd) then
   write(*,*)
   write(*,*) "Chosen mode: evaluation of AIMD trajectory."
   write(*,*) "OUTCAR file ",trim(outcar_name)," is read ..."
   write(*,*)
!
!     First, read in number of atoms and element symbols
!
   inc3=1
   el_nums=0
   frames_read=0
   open(unit=11,file=outcar_name,status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The ",trim(outcar_name)," file of the AIMD simulation is missing!"
      stop
   end if 
   do
      read(11,'(a)',iostat=readstat) a160
      if (readstat .ne. 0) exit
!
!     Number of atoms and element symbols of them
!
      if (index(a160,"ions per type =") .ne. 0) then
         read(a160,*,iostat=readstat) a80,a80,a80,a80,el_nums
         do j=1,10
            if (el_nums(j) .ne. 0) then
               nelems=j
            end if
         end do
         natoms=sum(el_nums(1:nelems))
         if (allocated(xyz)) deallocate(xyz)
         allocate(xyz(3,natoms))
         if (allocated(grad)) deallocate(grad)
         allocate(grad(3,natoms))
         if (allocated(at_names)) deallocate(at_names)
         allocate(at_names(natoms))
!
!     Determine element symbols for all atoms
!        
         at_names="XX"
         inc2=1
         do j=1,nelems
            do k=1,el_nums(j)
               at_names(inc2)=el_syms(j)
               inc2=inc2+1
            end do
         end do

      end if
      if (index(a160,"TITEL") .ne. 0) then
         read(a160,*) a80,a80,a80,el_syms(inc3)
         inc3=inc3+1
         if (inc3 .eq. nelems) exit
      end if
   end do
   close(11)
   write(*,*) "Number of atoms in the system: ",natoms
!
!    Second, read in the positions and forces for all frames
!
   frame_act=0
   open(unit=11,file=outcar_name,status="old")
   do
      read(11,'(a)',iostat=readstat) a160
      if (readstat .ne. 0) exit
!
!     Number of atoms and element symbols of them
!
      if (index(a160,"direct lattice vectors") .ne. 0) then
         frame_act = frame_act + 1
!
!     Only read in the frame if its in the correct reading frequency
!
         if (modulo(frame_act,read_freq) .eq. 0) then
            read(11,*) cell_vecs(1,:)
            read(11,*) cell_vecs(2,:)
            read(11,*) cell_vecs(3,:)
            frames_read=frames_read+1
            do
               read(11,'(a)',iostat=readstat) a160
               if (readstat .ne. 0) exit
               if (index(a160,"TOTAL-FORCE (eV/Angst)") .ne. 0) then
                  read(11,*)
                  do j=1,natoms
                     read(11,*) xyz(:,j),grad(:,j)
                  end do
!
!     Read the energy after the force section
!
                  do j=1,12
                     read(11,*)
                  end do
                  read(11,*) adum,adum,adum,adum,adum,adum,energy
                  exit
               end if
            end do
!
!     Write output for current frame (aenet)
!
            if (write_aenet) then
                  if (i .lt. 10) then
                     write(out_name,'(a,a,i1)') trim(basename),"_frame",i
                  else if (i .lt. 100) then
                     write(out_name,'(a,a,i2)') trim(basename),"_frame",i
                  else if (i .lt. 1000) then
                     write(out_name,'(a,a,i3)') trim(basename),"_frame",i
                  else
                     write(out_name,'(a,a,i4)') trim(basename),"_frame",i
                  end if
                  open(unit=30,file="xsf_files/" // trim(out_name) // ".xsf")
                  write(30,'(a,f14.8,a)') "# total energy = ",energy," eV"
                  write(30,*)
                  write(30,'(a)') "CRYSTAL"
                  write(30,'(a)') "PRIMVEC"
                  do j=1,3
                     write(30,'(3f14.8)') cell_vecs(j,:)
                  end do
                  write(30,'(a)') "PRIMCOORD"
                  write(30,'(i6,i1)') natoms,1
                  do j=1,natoms
                     write(30,'(a,a,6f14.8)') at_names(j)," ",xyz(:,j),grad(:,j)
                  end do
                  close(30)
!
!     Write name of output file to list of filenames
!
                  write(45,'(a,a,a)') "./xsf_files/", trim(out_name),".xsf"
            end if

!
!     Write output for current frame (MACE)
!
            if (write_mace) then
                  write(48,*) natoms
                  write(48,'(a,9f13.8,a,a,f20.10,a)') ' Lattice="',cell_vecs(1,:),cell_vecs(2,:),&
                          & cell_vecs(3,:),'" Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3',&
                          & ' Nmols=1 REF_energy=',energy,' pbc="T T T"'

                  do j=1,natoms
                     write(48,'(a,a,3f14.8,a,3f14.8)') at_names(j),"  ",xyz(:,j),"  0  ",grad(:,j)
                  end do
            end if
         else
            read(11,'(a)',iostat=readstat)
         end if
      end if
   end do
   close(11)
   write(*,*) "... done!"
   write(*,*) "Total number of AIMD frames in file: ",frame_act
   write(*,*) "Number of frames converted to training set: ",frames_read
   write(*,*)
end if


if (md_mode .eq. "setup") then
   write(*,*)
   write(*,*) " md_traj=setup: extract the frames of a XDATCAR file and generate folders"
   write(*,*) "  for subsequent VASP calculations."
   write(*,*)

   if (header .eqv. .true.) then
      write(*,*) " Before each frame a header with 8 lines containing information about"
      write(*,*) " the unit cell will be assumed (typical for NpT ensembles with a varying"
      write(*,*) " unit cell)."
      write(*,*)
   else
      write(*,*) " No header before each frame will be assumed (typical for NVT ensembles)."
      write(*,*)
   endif

!
!    Open the XDATCAR file and read in the frames
!

   open(unit=14,file="XDATCAR",status="old",iostat=openstat)
   if (openstat .ne. 0) then
      stop "ERROR! The file 'XDATCAR' could not been found!"
   end if
   read(14,*)
   read(14,*) factor   ! the global geometry conversion factor
   close(14)
!
!    Check if the XDATCAR has the format of NpT trajectory with the
!    full header for each frame, then, skip the headers in each read in
!
   open(unit=14,file="XDATCAR",status="old")
   do i=1,6
      read(14,'(a)') cdum
   end do
   el_nums=0
   read(14,'(a)') cdum
   read(cdum,*,iostat=readstat) el_nums
   read(14,*) cdum
   natoms=sum(el_nums)

   do i=1,natoms
      read(14,'(a)') cdum
   end do
   read(14,*)
   read(14,*) scale_dum
   if (abs(scale_dum-factor) .lt. 1E-10) then
      header=.true.
   else
      header=.false.
   end if
   close(14)

   if (header .eqv. .true.) then
      write(*,*) " The XDATCAR file comes from a NpT ensemble MD simulation."
      write(*,*)
   else
      write(*,*) " The XDATCAR files comes from a NVT ensemble MD simulation."
      write(*,*)
   endif

!
!    Read in the content of the XDATCAR file
!
   open(unit=56,file="XDATCAR",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*)
      write(*,*) "The XDATCAR file is not there!"
      write(*,*)
      stop
   end if
   write(*,*) "  Determine number of atoms and frames in XDATCAR..."
   call system("wc -l XDATCAR > xdat_length")
   open(unit=45,file="xdat_length",status="old")
   read(45,*) xdat_lines
   close(45)
   close(56)
   open(unit=56,file="XDATCAR",status="old")
   do i=1,6
      read(56,*)
   end do
   el_nums=0
   read(56,*,iostat=readstat) el_nums
   close(56)

   natoms=sum(el_nums)

   if (header .eqv. .false.) then
      nframes = (xdat_lines-7)/(natoms+1)
   else
      nframes = (xdat_lines)/(natoms+8)
   end if

   write(*,*) "  Number of atoms: ",natoms
   write(*,*) "  Number of frames: ",nframes
   write(*,*) "  ... completed!"
   write(*,*) " Read XDATCAR content ..."
   allocate(xdat_content(xdat_lines),stat=alloc_stat,errmsg=alloc_err)
   if (alloc_stat .ne. 0) then
      write(*,*) "Allocation of xdat_content array failed:",trim(alloc_err)
      write(*,*) "Please read in a smaller XDATCAR or increase memory!"
      stop
   end if


   open(unit=56,file="XDATCAR",status="old")
   do i=1,xdat_lines
      read(56,'(a)') xdat_content(i)
   end do
   close(56)
   write(*,*) " ... completed!"

!
!    Generate folders for each (selected) frame and write POSCAR file in it
!

   write(*,*) " Write folders for VASP single points ..."
   do k=1,nframes/traj_freq
      i=k*traj_freq
      if (k .lt. 10) then
         write(foldername,'(a,i1)') "frame",k
      else if (k .lt. 100) then
         write(foldername,'(a,i2)') "frame",k
      else if (k .lt. 1000) then
         write(foldername,'(a,i3)') "frame",k
      else
         write(foldername,'(a,i4)') "frame",k
      end if
 
      call system("mkdir "//trim(foldername))
      call chdir(foldername)
      open(unit=57,file="POSCAR",status="replace")

      if (header .eqv. .true.) then
         do j=(i-1)*(natoms+8)+1,i*(natoms+8)
            write(57,*) adjustl(trim(xdat_content(j)))
         end do

      else
         do j = 1,7
            write(57,*) adjustl(trim(xdat_content(j)))
         end do
         do j = 8 + (i-1)*(natoms+1),i*(natoms+1) + 7
            write(57,*) adjustl(trim(xdat_content(j)))
         end do
      end if

      close(57)
      call chdir("..")
   end do
!
!    Write a shell script which enables the simple copy of the other 
!      input files to the respective folders
!

   open(unit=58,file="copy_input.sh",status="replace")
   if (nframes/traj_freq .lt. 10) then
      write(58,'(a,i1,a)') "for i in {1..",nframes/traj_freq,"}"
   else if (nframes/traj_freq .lt. 100) then
      write(58,'(a,i2,a)') "for i in {1..",nframes/traj_freq,"}"
   else if (nframes/traj_freq .lt. 1000) then
      write(58,'(a,i3,a)') "for i in {1..",nframes/traj_freq,"}"
   else
      write(58,'(a,i4,a)') "for i in {1..",nframes/traj_freq,"}"
   end if
   write(58,'(a)') "do"
   write(58,*) "   cp POTCAR frame$i"
   write(58,*) "   cp KPOINTS  frame$i"
   write(58,*) "   cp INCAR  frame$i"
   write(58,*) "   cp slurm_script  frame$i"
   write(58,'(a)') "done"
   close(58)
   write(*,*) " ... completed!"
   write(*,*) " Modify and use script copy_input.sh to start calculations."
   write(*,*)
end if
!
!    If the prepared calculations of all frames are finished, read in the 
!    results and generate the training set files 
!
if (md_mode .eq. "eval") then
!
!    First, determine the number of frames
!
   i=1
   do
      if (i .lt. 10) then
         write(foldername,'(a,i1)') "frame",i
      else if (i .lt. 100) then
         write(foldername,'(a,i2)') "frame",i
      else if (i .lt. 1000) then
         write(foldername,'(a,i3)') "frame",i
      else
         write(foldername,'(a,i4)') "frame",i
      end if
 
      call chdir(trim(foldername),readstat)
      if (readstat .eq. 0) then
         call chdir("..")
         i=i+1
      else
         nframes=i-1
         exit
      end if        
   end do  
   write(*,*) 
   write(*,*) "Number of calculated frames: ",nframes

   allocate(skip_frame(nframes))
   write(*,*)
   write(*,*) "Read in all OUTCAR files and write training set file(s) ..."
!
!    Now, read in the OUTCAR files of all frames 
!
   loop_frames: do i=1,nframes
      if (i .lt. 10) then
         write(foldername,'(a,i1)') "frame",i
      else if (i .lt. 100) then
         write(foldername,'(a,i2)') "frame",i
      else if (i .lt. 1000) then
         write(foldername,'(a,i3)') "frame",i
      else
         write(foldername,'(a,i4)') "frame",i
      end if
      call chdir(trim(foldername))
      open(unit=78,file="OUTCAR",status="old",iostat=readstat)
      if (readstat .ne. 0) then
         write(*,*) "The OUTCAR file for frame ",i," could not been found!"
         write(*,*) "The frame will be skipped in the training set."
         call chdir("..")
         cycle loop_frames
      end if
      natoms=0
      el_nums=0
      inc3=1
      skip_frame(i)=.false.
      do  
         read(78,'(a)',iostat=readstat) a160
         if (readstat .ne. 0) exit

         if (index(a160,"ions per type =") .ne. 0) then
            read(a160,*,iostat=readstat) a80,a80,a80,a80,el_nums
            do j=1,10
               if (el_nums(j) .ne. 0) then
                  nelems=j
               end if
            end do
            natoms=sum(el_nums(1:nelems))
            if (allocated(xyz)) deallocate(xyz)
            allocate(xyz(3,natoms))
            if (allocated(grad)) deallocate(grad)
            allocate(grad(3,natoms))
            if (allocated(at_names)) deallocate(at_names)
            allocate(at_names(natoms))
!
!     Determine element symbols for all atoms
!        
            at_names="XX"
            inc2=1
            do j=1,nelems
               do k=1,el_nums(j)
                  at_names(inc2)=el_syms(j)
                  inc2=inc2+1
               end do
            end do
         end if

         if (index(a160,"TITEL") .ne. 0) then
            read(a160,*) a80,a80,a80,el_syms(inc3)
            inc3=inc3+1
         end if

         if (index(a160,"energy  without entropy=") .ne. 0) then
            read(a160,*) a80,a80,a80,a80,a80,a80,energy
         end if        

       !  if (index(a160,"FORCE on cell =-STRESS in cart") .ne. 0) then
       !     do j=1,13
       !        read(78,*)
       !     end do     
       !     read(78,*) stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(2,3),stress(1,3)
       !     stress(2,1)=stress(1,2)
       !     stress(3,2)=stress(2,3)
       !     stress(3,1)=stress(1,3)
       !  end if
!
!     Directly abort if an error has been writen by the VASP calculation
!
         if (index(a160,"Error") .ne. 0) then
            write(*,*) "The VASP calculation gave an error in frame ",i
            write(*,*) "The frame will be skipped in the training set."
            call chdir("..")
            cycle loop_frames
         end if

         if (index(a160,"POSITION") .ne. 0) then
            read(78,*)

            do j=1,natoms
               read(78,*,iostat=readstat) xyz(:,j),grad(:,j)
               if (readstat .ne. 0) then
                  write(*,*) "Something seems to be wrong with the OUTCAR file of frame ",i
                  write(*,*) "The frame will be skipped in the training set."
                  call chdir("..")
                  cycle loop_frames
               end if
               if ((grad(1,j) .gt. 50.0) .or. (grad(2,j) .gt. 50.0) .or. &
                        & (grad(3,j) .gt. 50.0)) then
                  write(*,*) "WARNING! A force component of atom ",j," in structure ",i," is huge!"
                  if (skip_large) then
                     skip_frame(i)=.false.
                  end if
               end if
            end do
!
!     Read in volume and basis vectors from POSCAR to avoid problem with missing spaces
!      for large numbers
!
            open(unit=127,file="POSCAR",status="old")
            read(127,*)
            read(127,*)
            read(127,*) cell_vecs(1,:)
            read(127,*) cell_vecs(2,:)
            read(127,*) cell_vecs(3,:)
            close(127)
         end if

!         if(index(a160,'VOLUME and BASIS-vectors').ne.0)then
!            read(78,*)
!            read(78,*)
!            read(78,*)
!            read(78,*)
!            read(78,*) cell_vecs(1,:)
!            read(78,*) cell_vecs(2,:)
!            read(78,*) cell_vecs(3,:)
!         end if
  
      end do
      if (natoms .lt. 1) then
         write(*,*) "Something seems to be wrong with the OUTCAR file of frame ",i
         write(*,*) "The frame will be skipped in the training set."
         call chdir("..")
         cycle loop_frames
      end if        
      close(78)

      call chdir("..")
!
!     Write output for current frame (aenet)
!
      if (write_aenet) then
         if (.not. skip_frame(i)) then
            if (i .lt. 10) then
               write(out_name,'(a,a,i1)') trim(basename),"_frame",i
            else if (i .lt. 100) then
               write(out_name,'(a,a,i2)') trim(basename),"_frame",i
            else if (i .lt. 1000) then
               write(out_name,'(a,a,i3)') trim(basename),"_frame",i
            else 
               write(out_name,'(a,a,i4)') trim(basename),"_frame",i
            end if   
            open(unit=30,file="xsf_files/" // trim(out_name) // ".xsf")
            write(30,'(a,f14.8,a)') "# total energy = ",energy," eV"
            write(30,*)
            write(30,'(a)') "CRYSTAL"
            write(30,'(a)') "PRIMVEC"
            do j=1,3
               write(30,'(3f14.8)') cell_vecs(j,:)
            end do
            write(30,'(a)') "PRIMCOORD"
            write(30,'(i6,i1)') natoms,1
            do j=1,natoms
               write(30,'(a,a,6f14.8)') at_names(j)," ",xyz(:,j),grad(:,j)
            end do
            close(30)
!
!     Write name of output file to list of filenames
!
            write(45,'(a,a,a)') "./xsf_files/", trim(out_name),".xsf"
         else
            write(*,*) "Frame ",i," is skipped in the output due to large force components."
         end if

      end if

!
!     Write output for current frame (MACE)
!
      if (write_mace) then
         if (.not. skip_frame(i)) then
            write(48,*) natoms
  !       write(48,'(a,9f13.8,a,9f13.8,a,a,f20.10,a)') 'REF_stress="',stress(1,:),stress(2,:),&
  !               & stress(3,:),'" Lattice="',cell_vecs(1,:),cell_vecs(2,:),&
  !               & cell_vecs(3,:),'" Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3',&
  !               & ' Nmols=1 REF_energy=',energy,' pbc="T T T"'
            write(48,'(a,9f13.8,a,a,f20.10,a)') ' Lattice="',cell_vecs(1,:),cell_vecs(2,:),&
                    & cell_vecs(3,:),'" Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3',&
                    & ' Nmols=1 REF_energy=',energy,' pbc="T T T"'

            if (xyz(3,65) .lt. 0.001) then
               write(*,*) xyz(3,j),j,i
            end if

            do j=1,natoms
               write(48,'(a,a,3f14.8,a,3f14.8)') at_names(j),"  ",xyz(:,j),"  0  ",grad(:,j)
            end do
         else
            write(*,*) "Frame ",i," is skipped in the output due to large force components."
         end if
      end if
   end do loop_frames
   write(*,*) "... completed!"


end if        

if (md_mode .ne. "setup") then
   if (write_aenet) then
      close(45)
      write(*,*)
      write(*,*) "XSF files written to folder 'xsf_files'"
      write(*,*) "List of file names written to 'xsf_list.dat'"
      write(*,*)

   end if

   if (write_mace) then
      close(48)
      write(*,*)
      write(*,*) "MACE training set file "//trim(basename)//".xyz written"
   end if
end if 

write(*,*) " vasp2trainset exited normally."
write(*,*) 

end program vasp2trainset    

