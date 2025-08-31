!
!    analyze_md: Analyze VASP or MLIP molecular dynamics simulations 
!      of bulk or surface slab systems from trajectories, calculates several
!      useful measures like averaged element densities orthogonal
!      to the surface, diffusion coefficients, radial distribution
!      functions and picks example structures for further calculations
!      like core level shifts or partial charges.
!    Part of VASP4CLINT
!     Julien Steffen, 2025 (julien.steffen@fau.de)
!

!
!    Contains all global variables, i.e., those that are used by the main
!     program and one or several subroutines and not transmitted explicitly
!     when they are called
!
module analyze_md_mod 
implicit none 
logical::write_traj  ! if the trajectory shall be written as xyz file
logical::dens_elems  ! if element density profiles shall be calculated
logical::track_atoms  ! if single atoms shall be tracked
logical::calc_rdf  ! if RDF of all element combinations shall be calculated
logical::surf_tension ! if the surface tension of a slab shall be calculated 
logical::calc_diff  ! if the element self-diffusion coefficients shall be calc.
logical::diff_2d   ! if diffusion along a 2d plane shall be considered only
logical::diff_collect   ! if the collective diffusion of one element 
logical::read_dt    ! control if the MD time step has been read in
logical::calc_stable  ! determines the stability of a system
logical::skip_xdat ! skip the read in of XDATCAR for surface tension
logical::dens_cube ! if spacially-resolved average element densities are given

real(kind=8)::pi  ! the pi
integer::analyze_parts ! number of parts of the trajectory to be evaluated
integer::nbins   ! number of density profile bins
integer,allocatable::track_list(:)  ! list of tracked atoms indices
character(len=5),allocatable::track_list_read(:) ! read in list
integer::track_num  ! number of tracked atoms
integer::rdf_bins   ! number of RDF bins
real(kind=8)::rdf_range ! maximum distance for RDF (Angstroms)
real(kind=8)::rdf_binsize  ! width of a single RDF bin (Angstroms)
integer::frame_first  ! first trajectory frame for evaluation
integer::frame_last   ! last trajectory frame for evaluation
real(kind=8)::z_shift   ! shift of z-coordinates for density plots
real(kind=8)::z_val_max  ! maximum z value for density plot, before image shift
character(len=2)::cls_element ! element symbol for which CLS are calculated
integer::cls_rounds ! Number of individual CLS evaluations
integer::atom_slices  ! number of CLS slices along axis
real(kind=8)::time_step  ! the MD time step in fs
character(len=2),allocatable::el_names(:) ! symbols of the elements
integer,allocatable::el_nums(:)  ! numbers of the elements
logical::eval_stat(10)  ! for print out of status for long processes
integer::all_tasks  ! for print out of status
integer::cls_elem_ind  ! The number of the element whose structures are picked
real(kind=8)::xlen,ylen,zlen,zmax  ! sizes of the unit cell
character(len=20)::dens_mode ! which kind of element density profile 
real(kind=8)::dens_depth  ! depth of the slab surface region (concentrations)
real(kind=8)::box_volume  ! volume of the unit cell
real(kind=8)::factor  ! the global POSCAR scaling factor 
real(kind=8)::bond_tol  ! tolerance factor for bond breaking
character(len=2),allocatable::at_names(:)  ! the element symbols
integer::natoms ! number of atoms in system
integer::nelems  ! number of elements in system

end module analyze_md_mod

program analyze_md
use analyze_md_mod

implicit none
integer::i,j,k
integer::xdat_lines
real(kind=8)::rdum1,rdum2
integer::frame_shift
integer::nframes  ! total number of frames
real(kind=8),allocatable::xyz(:,:,:) ! all coordinates of the trajectory
real(kind=8),allocatable::xyz_part(:,:,:)  ! coordinates of the current traj. part
character(len=2),allocatable::el_names_read(:) ! symbols of the elements
character(len=80)::folder_name  ! names of evaluatio/results folders
integer,allocatable::frames_part(:)
real(kind=8)::act_num(3)
character(len=120)::a120
character(len=220)::a220
character(len=1)::atest
character(len=150)::arg,cdum
logical::ana_present ! if the current analysis folder exists
real(kind=8),allocatable::xlens(:),ylens(:),zlens(:) ! NpT unit cell sizes
real(kind=8),allocatable::xlens_p(:),ylens_p(:),zlens_p(:) ! ... for call of subroutines
character(len=2)::el_eval_list(20)  ! the elements to be set to COM (density mode)
integer::el_eval_num  ! number of elements for the COM (density mode)
integer::readstat,openstat
integer::counter,endl 
real(kind=8)::scale_dum
logical::npt_format



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
write(*,*) "PROGRAM analyze_md: Evaluation of VASP DFT/MLIP trajectories"
write(*,*) " trajectories of bulk systems and surface slabs."
write(*,*) "A trajectory can either be evaluates as a whole (results in folder"
write(*,*) " 'analysis') or in several parts (results in folders analysis_part[i])"
write(*,*) "The file XDATCAR must be present!"
write(*,*) "The following command line arguments can/must be given (with - sign!):"
write(*,*)
write(*,*) "General settings and evaluation modes:"
write(*,*) " -overview:  print an overview of all scripts and programs in utils4VASP"
write(*,*) " -part_number=[number]: In how many parts the trajectory shall be divided"
write(*,*) "     and doing all analysis separately for each part."
write(*,*) " -rdf : The radial distribution functions for all element combinations  "
write(*,*) "     shall be calculated. Then, also the total number of neighbors "
write(*,*) "     up to a certain distance will be calculated."
write(*,*) " -dens_elems : The element densities along a axis or around a sphere are calculated."
write(*,*) " -diffusion : calculates the diffusion coefficient, for each element in the slab"
write(*,*) "     separately, via the mean square displacement (MSD)."
write(*,*) " -dt=[time step in fs]: The time step used for MD simulation, in fs."
write(*,*) " -track_atoms=[list of numbers] : Write time-dependent positions of chosen"
write(*,*) "     atoms to file. Example: track_atoms=1,78,178"
write(*,*) " -dens_cube: Write Gaussian cube file with spatially resolved atomic densities."
write(*,*) " -dft_element=[element]: DFT calculation templates (POSCAR) will be generated for"
write(*,*) "     the chosen element for representative parts of the trajectory/the system."
write(*,*) " -tension : calculates the surface tension averaged over all MD frames."
write(*,*) "     For this, the OUTCAR file needs to be present (MD with ISIF=2)"
write(*,*) " -stability : Determines if all bonds in the system are stable, and, if not,"
write(*,*) "     which bond broke first after which time."
write(*,*)
write(*,*) "List of additional options that can/must be added for certain evaluation modes:"
write(*,*) " -frame_first=[number] : First trajectory frame that shall be evaluated by "
write(*,*) "     the script (e.g., in order to skip equilibration parts) (default: 1)"
write(*,*) " -frame_last=[number] : Last trajectory frame that shall be evaluated "
write(*,*) "     (default: last frame in XDATCAR)."
write(*,*) " -write_traj : The file 'trajectory.xyz' containing all frames of XDATCAR"
write(*,*) "     shall be written during the analysis."
write(*,*) " -rdf_bins=[number] : Number of bins for RDF evaluation (default: 201)"
write(*,*) " -rdf_cutoff=[value]: Cutoff for RDF evalulation (default: 8 Angstrom)"
write(*,*) " -diff_collect : calculates the diffusion coefficients of each element, treating all"
write(*,*) "     atoms of them as one effective sum particle."
write(*,*) " -diff_2d : calculates the 2D-diffusion coefficient along x and y, for each element"
write(*,*) "     in the slab separately, via the mean square displacement (MSD)."
write(*,*) " -dens_bins=[number] : Number of bins for element densities (default: 501)"
write(*,*) " -dens_mode=[character] : Decides, along which topology the element densities "
write(*,*) "     shall be calculated, possible are x, y or z for the coordinate axes or "
write(*,*) "     sphere:[el1,...], where a sphere around the center of mass of the atoms of the "
write(*,*) "     chosen element(s) is formed, for example -dens_mode=sphere:Pt,H (default: z)."
write(*,*) " -dens_depth=[value] : How many Angstroms below the slab surface the surface "
write(*,*) "     concentration of elements shall be determined (default: 5.0)"
write(*,*) " -dft_slices=[number]: How many different slices along the z-coordinate where "
write(*,*) "     the atom for which DFT POSCAR templates are written (near or far from the "
write(*,*) "     surface of the slab) (default: 100)."
write(*,*) " -dft_rounds=[number]: In how many parts the trajectory shall be divided for"
write(*,*) "     DFT template generation in each part (default: 1)."
!
!    The pi
!
pi=3.141592653589793238

!
!    Read in the command line arguments and define the evaluation settings
!
call read_arguments()
!
!    First, determine the number of lines in the XDATCAR file
!
open(unit=45,file="XDATCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*)
   write(*,*) "The file XDATCAR with the trajectory is not there!"
   write(*,*)
   stop
end if   
call system("wc -l XDATCAR > xdat_length")
open(unit=45,file="xdat_length",status="old")
read(45,*) xdat_lines
close(45)
!
!    Open the XDATCAR file and read in the frames 
!
        
open(unit=14,file="XDATCAR",status="old",iostat=openstat)
if (openstat .ne. 0) then 
   stop "ERROR! The file 'XDATCAR' could not been found!"
end if        
read(14,*)
read(14,*) factor   ! the global geometry conversion factor 
!
!    Read in the cell dimensions: Cubic unit cell is always assumed!
!
read(14,*) xlen,rdum1,rdum2
read(14,*) rdum1,ylen,rdum2
read(14,*) rdum1,rdum2,zlen
!
!    Box volume for surface tension calculation
!
box_volume=xlen*ylen*zlen
!
!    Read in the elements (so far, only two different allowed)
!
natoms=0
allocate(el_names_read(10))
el_names_read="XX"
read(14,'(a)') cdum
read(cdum,*,iostat=readstat) el_names_read
nelems=0
do i=1,10
   if (el_names_read(i) .eq. "XX") exit
      nelems=nelems+1
   if (el_names_read(i) .eq. cls_element) then
      cls_elem_ind=nelems
   end if        
end do

allocate(el_names(nelems),el_nums(nelems))

do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(14,*) el_nums
close(14)
!
!     Number of atoms in the slab
!
natoms = sum(el_nums)
!
!    Check if the XDATCAR has the format of NpT trajectory with the 
!    full header for each frame, then, skip the headers in each read in
!
open(unit=14,file="XDATCAR",status="old")
do i=1,8
   read(14,'(a)') cdum
end do   
do i=1,natoms
   read(14,'(a)') cdum
end do
read(14,*)
read(14,*) scale_dum
if (abs(scale_dum-factor) .lt. 1E-10) then
   npt_format=.true.
else
   npt_format=.false.
end if        
close(14)
!
!    Open the XDATCAR again for a full 
!
open(unit=14,file="XDATCAR",status="old")
do i=1,7
   read(14,'(a)') cdum
end do
!
!    Define the element symbols for all atoms in the system
!
allocate(at_names(natoms))
nframes=0
counter = 1
do i=1,nelems
   do j=1,el_nums(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do
!
!    If it is a NpT trajectory, allocate the arrays for time-dependent 
!    unit cell sizes
!
if (npt_format) then
   nframes = int((xdat_lines)/(natoms+8))
   allocate(xlens(nframes),ylens(nframes),zlens(nframes))
   xlens(1)=xlen
   ylens(1)=ylen
   zlens(1)=zlen
else        
   nframes = int((xdat_lines - 7)/(natoms+1))
end if
!
!    If the XDATCAR file shall not be read in, skip the rest
!
if (skip_xdat) then
   close(14)
   goto 33
end if

if (npt_format) then
   close(14)
   open(unit=14,file="XDATCAR",status="old")
end if        

allocate(xyz(3,natoms,nframes))
!
!    Read in the coordinates of the frames and correct for image flags
!
eval_stat = .false.
write(*,*)
write(*,*) "Read in the trajectory from XDATCAR..."
do i=1,nframes 
!
!    Every 10% of the read in, give a status update 
!
   do j=1,10
      if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
         if (.not. eval_stat(j)) then
            write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
            eval_stat(j) = .true.
         end if
      end if
   end do
   read(14,*)
   if (npt_format) then
      read(14,*)
      read(14,*) xlens(i),rdum1,rdum2
      read(14,*) rdum1,ylens(i),rdum2
      read(14,*) rdum1,rdum2,zlens(i)
      read(14,*)
      read(14,*)
      read(14,*)
   end if        
   do j=1,natoms 
      read(14,'(a)') a120
      read(a120,*,iostat=readstat) act_num(:)
!
!   In long trajectories, problems might arise with double-digit negative numbers!
!   (no space between them)   Insert a space again!
!

      if (npt_format) then
         xlen=xlens(i)
         ylen=ylens(i)
         zlen=zlens(i)
      end if         
      if (readstat .ne. 0) then
         a220 = ""
         endl = 1
         do k=1,119
            atest=a120(k:k+1)
            if (index(atest,"-") .ne. 0) then
               write(a220(1:endl+4),'(a,a)') a220(1:endl),"   -"
               endl = endl + 4
            else        
               write(a220(1:endl+1),'(a,a)') a220(1:endl),atest
               endl = endl + 1
            end if   
         end do
    !  write(*,*) "Problem in XDATCAR file!"
         read(a220,*,iostat=readstat) act_num(:)
      end if    
       
      do k=1,3
         do 
            if (act_num(k) > 1d0) then
               act_num(k) = act_num(k) - 1d0
            else 
               exit
            end if   
         end do
         do 
            if (act_num(k) < 0d0) then
               act_num(k) = act_num(k) + 1d0   
            else      
               exit
            end if   
         end do
!
!     Convert to cartesian coordinates
!        
         if (k .eq. 1) then
            xyz(k,j,i) = act_num(k)*xlen
         end if
         if (k .eq. 2) then   
            xyz(k,j,i) = act_num(k)*ylen
         end if   
         if (k .eq. 3) then
            xyz(k,j,i) = act_num(k)*zlen
         end if   
      end do
   end do
end do
33 continue

write(*,*) " completed!"
close(14)

!
!    If frame_last was determined explicitly and is larger than frame_first: 
!     set it to the value
!     else: use the total number of frames for frame_last
!
if (frame_last .le. frame_first) then
   frame_last=nframes
end if
!


write(*,*)
write(*,*) "---------- SETTINGS ---------------------------------------"
write(*,'(a,i7)') " Number of atoms in the system:",natoms
write(*,'(a,i9)') " Number of frames in the trajectory:",nframes
write(*,'(a,i5)') " Number of separately evaluated trajectory parts:",analyze_parts
if (npt_format) then
   write(*,*) "The XDATCAR file is a NpT trajectory (variable volume)"
else
   write(*,*) "The XDATCAR file is a NVT trajectory (fixed volume)"
end if  
write(*,'(a,i9,a)') " Frame No. ",frame_first," is the first to be evaluated."
write(*,'(a,i9,a)') " Frame No. ",frame_last," is the last to be evaluated."
if (write_traj) then
   write(*,*) "The trajectory will be written to 'trajectory.xyz'."
end if

if (calc_rdf) then
   write(*,*) "Radial distribution functions of all element combinations"
   write(*,*) " as well as neighbor numbers will be calculated."
   write(*,'(a,i6)') "  * Number of RDF bins: ",rdf_bins
   write(*,'(a,f12.6)') "  * RDF cutoff in Angstroms: ",rdf_range
end if
if (track_atoms) then
   write(*,*) "The positions of the following atoms will be tracked:"
   write(*,'(a)',advance="no") "    "
   do i=1,track_num-1
      write(*,'(i6,a)',advance="no") track_list(i),", "
   end do
   write(*,'(i6,a)') track_list(track_num)
end if        
if (dens_elems) then
   write(*,*) "The one-dimensional density of elements will be evaluated."
   if (trim(dens_mode) .eq. "x") then
      write(*,*) " * Density is evaluated along the x-axis."
   else if (trim(dens_mode) .eq. "y") then
      write(*,*) " * Density is evaluated along the y-axis."  
   else if (trim(dens_mode) .eq. "z") then
      write(*,*) " * Density is evaluated along the z-axis."  
   else if (trim(dens_mode(:7)) .eq. "sphere:") then   
      el_eval_list="XX"
      el_eval_num=0
      read(dens_mode(8:),*,iostat=readstat) el_eval_list
      do i=1,20
         if (el_eval_list(i) .eq. "XX") exit
         el_eval_num=el_eval_num+1
      end do
      write(*,*) " * Density is evaluated in a sphere around the center of mass"
      write(*,'(a)',advance="no") "    by the elements: "
      do i=1,el_eval_num
         write(*,'(a,a)',advance="no") trim(el_eval_list(i)),"  "
      end do
      write(*,*)
   end if   
   write(*,*) " * Number of bins for evaluating the densities:",nbins
   write(*,*) " * Depth below the surface slab (if existent) for which element"
   write(*,*) "    concentrations will be evaluated:",dens_depth
end if
if (calc_diff) then
   write(*,*) "The element self-diffusion coefficients will be calculated."
   write(*,'(a,f12.6)') "  * The assumed MD time step (fs):",time_step
   if (diff_2d) then
      write(*,*) " * Only diffusions parallel to the x-y plane will be considered."
   end if
   if (diff_collect) then
      write(*,*) " * Atoms of one element are treated as collective particle."
   end if        
end if        
if (cls_element .ne. "XX") then
   write(*,'(a)') " DFT example structures will be printed to POSCAR files,"
   write(*,'(a,a,a)') "  with ",cls_element," atoms at different representative positions "
   write(*,*) " * Number of slices along z-axis for selection:",atom_slices
   write(*,*) " * Number of trajectory (sub-)parts for selection:",cls_rounds
end if
if (calc_stable) then
   write(*,*) "The stability of the system (its bonds) shall be investigated."
   write(*,'(a,f12.6,a)') "  * Maximum bong elongation before breaking: ",bond_tol," times"
end if        
write(*,*) "------------------------------------------------------------"
write(*,*)

!
!    Determine the number of frames in each part
!
allocate(frames_part(analyze_parts))
frames_part=0
do i=1,analyze_parts-1
   frames_part(i)=ceiling(real((frame_last-frame_first))/real(analyze_parts))
end do
if (analyze_parts .gt. 1) then
   frames_part(analyze_parts)=nframes-sum(frames_part(1:analyze_parts-1))
else 
   frames_part(1)=frame_last-frame_first+1
end if
!
!    For surface concentrations: write them into global file
!
if (dens_elems .and. (analyze_parts .gt. 1)) then
   open(unit=47,file="surf_concs_layer.dat",status="replace")
   write(47,*) "# This file contains the surface concentrations of elements, "
   write(47,*) "# measured in the first two layers of the surface (w.r.t. the "
   write(47,*) "# total element density), in percent."
   write(47,'(a)',advance="no") " # Part    "
   do i=1,nelems
      write(47,'(a,a,a,a)',advance="no") trim(el_names(i))," (layer1)   ",trim(el_names(i))," (layer2)   "
   end do
   write(47,*)
   open(unit=48,file="surf_concs_depth.dat",status="replace")
   write(48,*) "# This file contains the surface concentrations of elements, "
   write(48,'(a,f12.6,a)') " # measured in the region up to ",dens_depth," Ang. below the "
   write(48,*) "# slab surfaces, in percent."
   write(48,'(a)',advance="no") " # Part    "
   do i=1,nelems
      write(48,'(a,a)',advance="no") trim(el_names(i)),"                "
   end do
   write(48,*)
end if        

!
!    Now loop over all parts of the trajectory and do all
!    assigned evaluations for each part! 
!    Call all subroutines and only give the current part 
!    of the xyz file (and the number of frames) to it
!
frame_shift=0
do i=1,analyze_parts
   if (allocated(xyz_part)) deallocate(xyz_part)
   if (allocated(xlens_p)) deallocate(xlens_p)
   if (allocated(ylens_p)) deallocate(ylens_p) 
   if (allocated(zlens_p)) deallocate(zlens_p)
   allocate(xyz_part(3,natoms,frames_part(i)))
   allocate(xlens_p(frames_part(i)))
   allocate(ylens_p(frames_part(i)))
   allocate(zlens_p(frames_part(i)))
   do j=1,frames_part(i)
      xyz_part(:,:,j)=xyz(:,:,frame_first-1+j+frame_shift)
      if (npt_format) then
!
!    For NpT trajectories: transfer the time-dependent unit cell size
!    to the evaluation parts as well!
!
         xlens_p(j)=xlens(frame_first-1+j+frame_shift)
         ylens_p(j)=ylens(frame_first-1+j+frame_shift)
         zlens_p(j)=zlens(frame_first-1+j+frame_shift)
      else
         xlens_p(j)=xlen
         ylens_p(j)=ylen
         zlens_p(j)=zlen
      end if   
            
   end do
   frame_shift=frame_shift+frames_part(i)
!
!    Generate (if not existent) and move into sub-calculation folder
!
   if (analyze_parts .eq. 1) then
      write(folder_name,'(a)') "analysis"
   else
      if (i .lt. 10) then
         write(folder_name,'(a,i1)') "analysis_part",i
      else if (i .lt. 100) then
         write(folder_name,'(a,i2)') "analysis_part",i
      else if (i .lt. 1000) then
         write(folder_name,'(a,i3)') "analysis_part",i
      else
         write(folder_name,'(a,i4)') "analysis_part",i
      end if
      write(*,*) "-------------------------------------------------"
      write(*,'(a,a,a,i4,a)') "  Analysis for part ",&
                   & trim(folder_name(14:))," of ",analyze_parts,":"
      write(*,*) "-------------------------------------------------"
   end if
   inquire(file=trim(folder_name),exist=ana_present)
   if (ana_present) call system ("rm -r "//trim(folder_name))
   call system ("mkdir "//trim(folder_name))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Now call the calculation parts if they are ordered
! 
   call chdir(trim(folder_name)) 
!
!    Calculate the element densities along the chosen axis
!
   if (dens_elems) then
      write(47,'(i5)',advance="no") i
      write(48,'(i5)',advance="no") i
      call calculate_densities(frames_part(i),xyz_part)
   end if
!
!    Perform the tracking of chosen atoms  
!
   if (track_atoms) then
      call track_atom_pos(frames_part(i),xyz_part)
   end if  
!
!    Calculate radial distribution functions and number of neighbors
!
   if (calc_rdf) then
      call calculate_rdf(frames_part(i),xyz_part,xlens_p,ylens_p,zlens_p)
   end if
!
!    Calculate the surface tension
!
   if (surf_tension) then
      call surface_tension(frames_part(i))
   end if
!
!    Calculate the self-diffusion coefficients for all elements
!
   if (calc_diff) then
      call calculate_diffusion(frames_part(i),xyz_part,xlens_p,ylens_p,zlens_p)
   end if
!
!    Pick structures with placed atoms of chosen element (CLS)
!
   if (cls_element .ne. "XX") then
      call pick_structures(frames_part(i),xyz_part)
   end if
!
!    Calculate the spacial densities of all elements
!
   if (dens_cube) then
      call spacial_density(frames_part(i),xyz_part,xlens_p,ylens_p,zlens_p)
   end if
!
!    Determine if the system remained stable during the MD and if not, which
!    bond broke first at what time
!
   if (calc_stable) then
      call stability(frames_part(i),xyz_part,xlens_p,ylens_p,zlens_p)
   end if        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call chdir("..")
!
!    If only dens_elems is active and more than 100 parts shall be 
!     evaluated: remove the single folders!
!
   if (dens_elems .and. (.not. track_atoms) .and. (.not. calc_rdf) .and. &
                 & (.not. surf_tension) .and. (.not. calc_diff)  .and. &
                 & (.not. dens_cube) .and. (cls_element .ne. "XX") .and. &
                 & (analyze_parts .gt. 100)) then
      call system ("rm "//trim(folder_name))
   end if
end do

if (dens_elems .and. (analyze_parts .gt. 1)) then
   close(47)
   close(48)
end if    
!
!    Write the trajectory in xyz format to file 
!
if (write_traj) then
   write(*,*)
   write(*,*) "Write the trajectory of the system to 'trajectory.xyz'..."     
   open(unit=15,file="trajectory.xyz",status="replace")
   eval_stat=.false.
   do i=frame_first,frame_last 
      do j=1,10
         if (real(i)/real(frame_last-frame_first) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      write(15,*) natoms 
      write(15,*) "Frame No.",i
      do j=1,natoms 
         write(15,*) at_names(j), xyz(:,j,i)
      end do
   end do
   close(15)
   write(*,*) " completed!"
   write(*,*)
end if
write(*,*)
if (dens_elems .and. (.not. track_atoms) .and. (.not. calc_rdf) .and. &
              & (.not. surf_tension) .and. (.not. calc_diff)  .and. &
              & (.not. dens_cube) .and. (cls_element .ne. "XX") .and. & 
              & (analyze_parts .gt. 100)) then
   write(*,*) "Only -dens_elems was ordered and more than 100 trajectory parts"
   write(*,*) " were evaluated. Therefore, no analysis folders are written."
else
   if (analyze_parts .eq. 1) then
      write(*,*) "Results of the analysis are written to folder analysis."
   else
      write(*,'(a)',advance="no") " Results of the analysis are written to folders analysis_part1 to"
      if (analyze_parts .lt. 10) then
         write(*,'(a,i1)') " analysis_part",analyze_parts    
      else if (analyze_parts .lt. 100) then
         write(*,'(a,i2)') " analysis_part",analyze_parts   
      else if (analyze_parts .lt. 1000) then
         write(*,'(a,i3)') " analysis_part",analyze_parts
      else
         write(*,'(a,i4)') " analysis_part",analyze_parts
      end if   
   end if        
end if        

write(*,*)
write(*,*) "analyze_md ended normally..."
write(*,*)
end program analyze_md










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE READ_ARGUMENTS
!    read in all command line arguments and check them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_arguments ()
use analyze_md_mod
implicit none 
integer::i,readstat
character(len=150)::arg
!
!    The number of separate parts in which a trajectory shall
!    be evaluated
!
analyze_parts=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-part_number=") then
      read(arg(14:),*) analyze_parts
   end if
end do

!
!    If the trajectory shall not be written out as xyz
!

write_traj = .true.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-notraj") then
      write_traj = .false.
   end if
end do
!
!    If the densities of elements along the z axis shall be calculated
!
dens_elems=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-dens_elems") then
      dens_elems = .true.
   end if
end do
!
!    For bin packing of abundancies/densities along z-axis
!
nbins=501
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-dens_bins=") then
      read(arg(12:),*) nbins
   end if
end do

!
!    Along which axis/sphere the element density shall be evaluated
!
!
!    The element for which the CLS shall be calculated
!
dens_mode="z"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-dens_mode=") then
      read(arg(12:),*) dens_mode
   end if
end do

!
!    Read in the time step in fs
!
dens_depth = 5.0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-dens_depth=") then
      read(arg(13:),*) dens_depth
      write(*,*)
   end if
end do
!
!    Track the time-dependent positions of one or several atoms, given
!    by their indices/numbers in the system
!    Up to 50 atoms can be chosen
!
track_atoms=.false.
allocate(track_list_read(50))
track_list_read="XXXXX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-track_atoms=") then
      read(arg(14:),*,iostat=readstat) track_list_read
      track_atoms=.true.
   end if
end do
allocate(track_list(50))
track_list=0
track_num=0
do i=1,50
   if (track_list_read(i) .eq. "XXXXX") exit
   read(track_list_read(i),*,iostat=readstat) track_list(i)
   if (readstat .ne. 0) then
      write(*,*) "The format of given atom indices in -track_atoms is wrong!"
      stop
   end if
   track_num=track_num+1
end do
!
!    Look if Radial Distribution functions shall be calculated
!
calc_rdf = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-rdf") then
      calc_rdf = .true.
   end if
end do
!
!    Look if the stability of the system shall be investigated
!
calc_stable = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-stability") then
      calc_stable = .true.
   end if
end do


!
!    For further settings of rdf calculations
!
!    For RDF calculation
!
rdf_bins = 200  ! number of RDF bins
rdf_range = 8.d0  ! maximum distance for RDF

do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-rdf_bins=") then
      read(arg(11:),*) rdf_bins
   end if
end do

do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-rdf_cutoff=") then
      read(arg(13:),*) rdf_range
   end if
end do

rdf_binsize = rdf_range/rdf_bins

!
!    For (optional) diffusion coefficient calculation
!
!
!    If the evaluation shall start not at the first frame but after
!     an equilibration phase
!
frame_first = 1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-frame_first=") then
      read(arg(14:),*) frame_first
   end if
end do
!
!    If the evaluation shall only go to a certain frame number and skip
!     the frames thereafter
!
frame_last = 0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-frame_last=") then
      read(arg(13:),*) frame_first
   end if
end do


!
!    The element for which the POSCAR files usable for further DFT calculations
!    shall be written
!
cls_element="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-dft_element=") then
      read(arg(14:),*) cls_element
   end if
end do
!
!    The number of individual core level shift evaluations
!
cls_rounds = 1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-dft_rounds=") then
      read(arg(13:),*) cls_rounds
   end if
end do

!
!    The number of slices along the z-axis for which CLS values shall be evaluated
!
atom_slices = 100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-dft_slices=") then
      read(arg(13:),*) atom_slices
   end if
end do

!
!    Activates the calculation of surface tension if desired
!
surf_tension=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-tension") then
      surf_tension=.true.
   end if
end do

!
!    Activates the calculation of diffusion coefficients
!

calc_diff=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-diffusion") then
      calc_diff=.true.
   end if
end do

!
!    Activates the calculation of two-dimensional diffusion coefficients
!

diff_2d=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-diff_2d") then
      calc_diff=.true.
      diff_2d=.true.
   end if
end do

!
!    Activates the calculation of collective diffusions, treating all atoms
!     of an element together as an effective particle
!

diff_collect=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-diff_collect") then
      calc_diff=.true.
      diff_collect=.true.
   end if
end do

!
!    Activates the calculation of spacially-resolved element densities,
!     written out as Gaussian Cube files
!
dens_cube = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-dens_cube") then
      dens_cube = .true.
   end if
end do
!
!    Read in the time step in fs
!
read_dt = .false.
if (calc_diff) then
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:4))  .eq. "-dt=") then
         read_dt = .true.
         read(arg(5:32),*) time_step
      end if
   end do
end if

if (calc_diff) then
   if (.not. read_dt) then
      stop "Please set for diffusion a time step with the -dt=... flag!"
   end if
end if

!
!    Look if the xyz trajectory file shall be written
!
write_traj = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-write_traj") then
      write_traj = .true.
   end if
end do
!
!    Look if Radial Distribution functions shall be calculated
!
calc_rdf = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-rdf") then
      calc_rdf = .true.
   end if
end do

!
!    Read in the tolerance factor for bond breakings (stability)
!
bond_tol = 1.7d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-bond_tol=") then
      read(arg(11:32),*) bond_tol
   end if
end do


!
!    If a surface tension calculation is requested and the flag -notraj
!    is turned on (but no other jobs), skip the read in of the XDATCAR
!    file alltogether!
!
skip_xdat = .false.
if ((.not. write_traj) .and. (.not. calc_diff) .and. (.not. calc_rdf) .and. surf_tension) then
   write(*,*) "The calculation of surface tension was required as only"
   write(*,*) " job and no trajectory shall be written out."
   write(*,*) "Therefore, the read in of the XDATCAR file will be skipped!"
   skip_xdat = .true.
end if

!
!    Foe diffusion calculations, the time step must be set!
!
if (calc_diff .and. (time_step .le. 0.0001)) then
   stop "Please give the MD time step in fs!"
end if

!
!     If no job has been chosen, abort the program
!
if ((.not. calc_rdf) .and. (.not. write_traj) .and. (.not. dens_elems) &
          & .and. (.not. track_atoms) .and. (.not. calc_diff) &
          & .and. (cls_element .eq. "XX") .and. (.not. calc_stable)) then
   write(*,*)
   write(*,*) "Please choose at least one of the evaluation options!"
   write(*,*)
   stop
end if

end subroutine read_arguments

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE TRACK_ATOM_POS
!    Print out time-dependent positions of selected atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine track_atom_pos(nframes,xyz)
use analyze_md_mod
implicit none
integer::i,j
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
character(len=60)::track_name
!
!    Write time-dependent positions of selected atoms to files, one file
!    for each selected atom!
!
write(*,*)
write(*,*) "Print time-dependent coordinates of chosen atoms to files..."

do i=1,track_num
   if (track_list(i) .le. 9) then
      write(track_name,'(a,i1,a)') "track_atom",track_list(i),".dat"
   else if (track_list(i) .le. 99) then
      write(track_name,'(a,i2,a)') "track_atom",track_list(i),".dat"
   else if (track_list(i) .le. 999) then
      write(track_name,'(a,i3,a)') "track_atom",track_list(i),".dat"
   else if (track_list(i) .le. 9999) then
      write(track_name,'(a,i4,a)') "track_atom",track_list(i),".dat"
   else
      write(track_name,'(a,i5,a)') "track_atom",track_list(i),".dat"
   end if
   write(*,'(a,i5,a,a,a)') "  * Track atom ",track_list(i)," (",trim(track_name),") ..."
   open(unit=86,file=track_name,status="replace")
   write(86,'(a,i5)') " # In this file, the time-dependent position of atom ",i
   write(86,*) "# is tracked, all coordinates are in Angstroms."
   if (read_dt) then
      write(86,*) "# time step (fs)            x-coordinate              y-coordinate     &
           &         z-coordinate"
   else
      write(86,*) "# frame No.                 x-coordinate              y-coordinate     &
           &         z-coordinate"
   end if
   do j=1,nframes
      if (read_dt) then
         write(86,*) time_step*j,xyz(:,track_list(i),j)
      else
         write(86,*) j,xyz(:,track_list(i),j)
      end if
   end do
   close(86)
end do
write(*,*) " completed!"
write(*,*)

end subroutine track_atom_pos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE CALCULATE_RDF
!    Calcuate radial distribution functions and integrated 
!    number of neighbors for all element combinations in 
!    the system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_rdf(nframes,xyz,xlens,ylens,zlens)
use analyze_md_mod
implicit none
integer::i,j,k,l,m
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
real(kind=8)::xlens(nframes),ylens(nframes),zlens(nframes)  ! the NpT unit cells
integer::task_act

integer::ig,ngr,npart1,npart2
real(kind=8)::nid,r_act,rho,vb

real(kind=8),allocatable::rdf_plot(:,:,:)
real(kind=8),allocatable::neighnum(:,:,:)
real(kind=8)::pos1(3),pos2(3),diff_vec(3)
real(kind=8),allocatable::rdf_sum(:,:,:)
real(kind=8)::dist
real(kind=8)::box_volumes(nframes)  ! the time-dependent unit cell volume
!
!    Calculate the RDFs of all elements if desired
!
write(*,*) "Calculate the RDFs of all element combinations..."
eval_stat=.false.
allocate(rdf_plot(rdf_bins,nelems,nelems))
allocate(neighnum(rdf_bins,nelems,nelems))
allocate(rdf_sum(rdf_bins,nelems,nelems))
rdf_plot=0.d0
rdf_sum=0.d0
neighnum=0.d0
task_act=0
!
!    Determine unit cell volumes for each MD step
!
do i=1,nframes
   box_volumes(i)=xlens(i)*ylens(i)*zlens(i)
end do

all_tasks=((nelems**2-nelems)/2+nelems)*(nframes)
do l=1,nelems
   do m=l,nelems
      do i=1,nframes
         task_act=task_act+1
!
!    Every 10% of the process, give a status update
!
         do j=1,10
            if (real(task_act)/real(all_tasks) .gt. real(j)*0.1d0) then
               if (.not. eval_stat(j)) then
                  write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
                  eval_stat(j) = .true.
               end if
            end if
         end do

         do j=1,el_nums(l)  ! Ni atoms
            do k=1,el_nums(m)  ! Ga atoms
               if (l .gt. 1) then
                  pos1 = xyz(:,sum(el_nums(1:l-1))+j,i)
               else
                  pos1 = xyz(:,j,i)
               end if
               if (m .gt. 1) then
                  pos2 = xyz(:,sum(el_nums(1:m-1))+k,i)
               else
                  pos2 = xyz(:,k,i)
               end if

               diff_vec=pos1-pos2
!
!     Correct the x component
!
               do while (abs(diff_vec(1)) .gt. xlens(i)/2.d0)
                  diff_vec(1)=diff_vec(1)-sign(xlens(i),diff_vec(1))
               end do

!
!     Correct the y component
!
               do while (abs(diff_vec(2)) .gt. ylens(i)/2.d0)
                  diff_vec(2)=diff_vec(2)-sign(ylens(i),diff_vec(2))
               end do

!
!     Correct the z component
!
               do while (abs(diff_vec(3)) .gt. zlens(i)/2.d0)
                  diff_vec(3)=diff_vec(3)-sign(zlens(i),diff_vec(3))
               end do

               dist = sqrt((diff_vec(1))**2 + &
                     & (diff_vec(2))**2 + (diff_vec(3))**2)

!
!     Remainder of the calculation taken from Frenkel Smit, page 86
!
               ig = int(dist/rdf_binsize)
               if ((ig .le. rdf_bins) .and. (ig .ge. 1)) then
                  rdf_plot(ig,l,m) = rdf_plot(ig,l,m) + 2.d0
                  rdf_sum(ig,l,m)=rdf_sum(ig,l,m)+ 1.d0
               end if
            end do
         end do
      end do
      do j=1,rdf_bins
         ngr=frame_last-frame_first
         npart1=el_nums(l)  ! which of both elements?
         npart2=el_nums(m)
         r_act=rdf_binsize*(real(j)+0.5d0)
         vb=((real(j) + 1.d0)**3-real(j)**3)*rdf_binsize**3
         rho=1.d0/(abs(box_volumes(i)))
         nid=4.d0/3.d0*pi*vb*rho*2.0d0
         rdf_plot(j,l,m)=rdf_plot(j,l,m)/(real(ngr)*real(npart1)*real(npart2)*real(nid))
         rdf_plot(j,m,l)=rdf_plot(j,l,m)
         rdf_sum(j,l,m)=rdf_sum(j,l,m)/ngr
         rdf_sum(j,m,l)=rdf_sum(j,l,m)
      end do
   end do
end do
rdf_plot(1,:,:)=0.d0
call system("mkdir RDFs")
call chdir("RDFs")
do i=1,nelems
   do j=1,nelems
!
!    Write the RDF for the current element combination
!
      open(unit=13,file="RDF_"//trim(el_names(i))//"-"//trim(el_names(j))//".dat", &
                    & status="replace")
      write(13,'(a,a,a,a,a)') "#      Distance (A)        RDF ("&
                    &//trim(el_names(i))//"-"//trim(el_names(j))//") "
      do k=1,rdf_bins
         write(13,'(2f19.10)') k*rdf_binsize,rdf_plot(k,i,j)
      end do
!
!    Write the integrated RDF (i.e., the number of surrounding atoms) for the
!    current element combination)
!
      open(unit=13,file="RDF_int_"//trim(el_names(i))//"-"//trim(el_names(j))//".dat", &
                    & status="replace")
      write(13,'(a,a,a,a,a)') "#      Distance (A)       integrated RDF ("&
                    &//trim(el_names(i))//"-"//trim(el_names(j))//") "

      do k=1,rdf_bins
         write(13,'(2f19.10)') k*rdf_binsize,sum(rdf_sum(1:k,i,j))/el_nums(i)
      end do

   end do
end do
call chdir("..")
write(*,*) "Plots of RDFs and integrated RDFs (number of neighbors) for "
write(*,*) " each element combination written in folder RDFs!"

end subroutine calculate_rdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE CALCULATE_DENSITIES
!    Calcuate elemental densities along a given axis (x,y,z)
!    and calculate the surface and subsurface concentrations 
!    of elements (only for slabs!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_densities(nframes,xyz)
use analyze_md_mod
implicit none
integer::i,j,k,l
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
integer::nmins
integer::counter
integer::eval_dim
real(kind=8)::normalize
real(kind=8),allocatable::min_pos(:)
real(kind=8),allocatable::z_vals(:)
real(kind=8),allocatable::int_side(:,:),tot_side(:)
real(kind=8)::z_min_lower1,z_min_lower2,z_min_upper1,z_min_upper2
real(kind=8),allocatable::z_dens(:,:),z_dens_tot(:)
real(kind=8)::zlo,zhi,zdiff,zstep  ! borders of z-density bins 
integer::slab_edge_lo,slab_edge_hi  ! the total elongation of the slab
integer::slab_int_depth  ! the width of dens_depth in bin multiples
real(kind=8),allocatable::concs_depth(:)  ! concentrations in dens_depth
character(len=2)::el_eval_list(20)  ! the elements to be set to COM
integer::el_index(10) ! reference elements for spherical analysis
integer::el_eval_num  ! number of reference elements
real(kind=8)::com_act(3)  ! the center of mass in the current frame
real(kind=8)::diff_vec(3)  ! distance between COM and any atom
real(kind=8)::diff_len  ! absolute value of the distance vector
integer::readstat ! status of read in files
!
!     Allocate arrays for total and element-wise densities
! 
allocate(z_dens(nbins,nelems))
allocate(z_dens_tot(nbins))
!
!     Do evaluations for the case that a coordinate axis has been chosen
!      as dimension/topology along which the density shall be evaluated
!
if (trim(dens_mode) .eq. "x" .or. trim(dens_mode) .eq. "y" .or. &
          & trim(dens_mode) .eq. "z") then
!  
!    Determine lowest and highest value of the coordinates along the 
!     chosen axis
!    MOD: Now o to zmax from POSCAR header
!  
   if (trim(dens_mode) .eq. "x") then
      zlo = 0.d0 
      zhi = xlen
      normalize=ylen*zlen
      eval_dim=1 
   else if (trim(dens_mode) .eq. "y") then
      zlo = 0.d0 
      zhi = ylen
      normalize=xlen*zlen
      eval_dim=2 
   else if (trim(dens_mode) .eq. "z") then
      zlo = 0.d0 
      zhi = zlen
      normalize=xlen*ylen
      eval_dim=3 
   end if
   zdiff = zhi-zlo
   zstep = zdiff/(nbins-1)

!
!    Evaluate the distribution of atoms along the chosen coordinate axis
!
   write(*,*) "Calculate element density distributions along ",trim(dens_mode),"-axis..."
   z_dens = 0.d0
   allocate(z_vals(nbins))
   z_vals=0.d0
   do i=1,nframes
      counter = 1
      do j=1,nelems
         do k=1,el_nums(j)
            do l=1,nbins-1
               if ((xyz(eval_dim,counter,i) .ge. zlo+(l-1)*zstep) .and. (xyz(eval_dim,counter,i) &
                            & .le. zlo+l*zstep)) then
                  z_dens(l,j) = z_dens(l,j) + 1.d0
               end if
            end do
            counter = counter + 1
         end do
      end do
   end do
   do i=1,nbins-1
      z_dens_tot(i)=sum(z_dens(i,:))
   end do
   write(*,*) " completed!"
!
!    Write the density profile to file
!
   open(unit=16,file="dens_elems.dat",status="replace")
   write(16,'(a)',advance="no") " # z-coordinate    "
   do i=1,nelems
      write(16,'(a,a)',advance="no") el_names(i),"      "
   end do
   write(16,*) "    total  "
   do i=1,nbins-1
      z_vals(i)=zlo+(i-0.5d0)*zstep
      z_dens_tot(i)=z_dens_tot(i)/((nframes)*normalize*zstep)
      z_dens(i,:)=z_dens(i,:)/((nframes)*normalize*zstep)
      write(16,*) z_vals(i),z_dens(i,:),z_dens_tot(i)

   end do
   close(16)
!
!    Determine surface concentrations of elements: The parts of the density profile
!     between the last local minimum and the asymptotics as well as of the penultimate
!     local minimum and the asymptotics are integrated
!     This is only done for binary and ternary systems
!
!    First, determine the local minima of the total density profile
!

   if (nelems .gt. 1) then
      write(*,*)
      if (nelems .eq. 2) then
         write(*,*) "Surface concentration of the second element will be determined..."
      else if (nelems .eq. 3) then
         write(*,*) "Surface concentrations of the second and third element will be determined..."
      end if
      allocate(min_pos(nbins))
      allocate(int_side(2,nelems))
      allocate(tot_side(2))
      int_side=0.d0
      tot_side=0.d0
      nmins=0
      do i=14,nbins-13
!
!     Only evaluate parts of the profile that are high enough to event processing of
!     numerical noise or detached atoms
!
         if (z_dens_tot(i) .gt. 0.01d0) then
            if ((z_dens_tot(i) .lt. z_dens_tot(i+1)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-1)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+2)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-2)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+3)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-3)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+4)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-4)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+5)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-5)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+6)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-6)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+7)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-7)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+8)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-8)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+9)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-9)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+10)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-10)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+11)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-11)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+12)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-12)) &
             &   .and. (z_dens_tot(i) .lt. z_dens_tot(i+13)) .and. (z_dens_tot(i) .lt. z_dens_tot(i-13))) then
                nmins=nmins+1
                min_pos(nmins)=z_vals(i)
            end if
         end if
      end do


      open(unit=39,file="surf_concs.dat",status="replace")
      write(39,'(a)') "# This file contains the concentration of different elements in the "
      write(39,'(a)') "# surface region of the slab system analyzed with analyzed_md"
      write(39,*)

      if (nelems .eq. 2) then
         write(39,'(a)') "# The concentration of the second element is calculated."
      else if (nelems .eq. 3) then
         write(39,'(a)') "# The concentrations of the second and the third element are calculated."
      end if
      z_min_lower1 = min_pos(1)
      z_min_lower2 = min_pos(2)
      z_min_upper1 = min_pos(nmins)
      z_min_upper2 = min_pos(nmins-1)
!
!     Now calculate the integrated densities
!
      do i=1,nbins
!
!     Outside the outer minima
!
         if (nelems .eq. 2) then
            if (z_vals(i) .lt. z_min_lower1)then ! .or. z_vals(i) .gt. z_min_upper1) then
               int_side(1,1)=int_side(1,1)+z_dens(i,2)
               tot_side(1)=tot_side(1)+z_dens_tot(i)
            end if
         else if (nelems .eq. 3) then
            if (z_vals(i) .lt. z_min_lower1 .or. z_vals(i) .gt. z_min_upper1) then
               int_side(1,1)=int_side(1,1)+z_dens(i,2)
               int_side(1,2)=int_side(1,2)+z_dens(i,3)
               tot_side(1)=tot_side(1)+z_dens_tot(i)
            end if
         end if
!
!     Outside the penultimate minima
!
         if (nelems .eq. 2) then
            if (z_vals(i) .lt. z_min_lower2 .and. z_vals(i) .gt. z_min_lower1) then
               int_side(2,1)=int_side(2,1)+z_dens(i,2)
               tot_side(2)=tot_side(2)+z_dens_tot(i)
            end if
            if (z_vals(i) .gt. z_min_upper2 .and. z_vals(i) .lt. z_min_upper1) then
               int_side(2,1)=int_side(2,1)+z_dens(i,2)
               tot_side(2)=tot_side(2)+z_dens_tot(i)
            end if
         else if (nelems .eq. 3) then
            if (z_vals(i) .lt. z_min_lower2 .and. z_vals(i) .gt. z_min_lower1) then
               int_side(2,1)=int_side(2,1)+z_dens(i,2)
               int_side(2,2)=int_side(2,2)+z_dens(i,3)
               tot_side(2)=tot_side(2)+z_dens_tot(i)
            end if
            if (z_vals(i) .gt. z_min_upper2 .and. z_vals(i) .lt. z_min_upper1) then
               int_side(2,1)=int_side(2,1)+z_dens(i,2)
               int_side(2,2)=int_side(2,2)+z_dens(i,3)
               tot_side(2)=tot_side(2)+z_dens_tot(i)
            end if
         end if
      end do
!
!     Second method for surface concentration determination: integrate regions
!     on both sides of the slab that are up to a given distance below the 
!     surface
!
!     First determine total elongation of the slab
!
      slab_int_depth=int(dens_depth/zdiff*real(nbins))
!
!     The lower side of the slab
!
      do i=1,nbins
         if (z_dens_tot(i) .gt. 1E-10) then
            slab_edge_lo = i
            exit
         end if         
      end do
      if (slab_edge_lo .eq. 1) then
         write(*,*) "The slab has no vacuum below it!"    
         write(*,*) "Please shift the slab first with modify_xdatcar!"
         stop 
      end if   
!
!     The upper side of the slab
!
      do i=nbins,1,-1
         if (z_dens_tot(i) .gt. 1E-10) then
            slab_edge_hi = i
            exit
         end if
      end do
      if (slab_edge_lo .eq. nbins) then
         write(*,*) "The slab has no vacuum above it!"   
         write(*,*) "Please shift the slab first with modify_xdatcar!"
         stop
      end if   
!
!     Abort if the integration range is too deep (lower and upper edges overlap)
!
      if ((slab_edge_lo+slab_int_depth) .ge. (slab_edge_hi-slab_int_depth)) then
         write(*,*) "Lower and upper edge regions of the slab overlap!"
         write(*,*) "Please choose a smaller -dens_depth parameter!"
         stop
      end if        

!
!     Now sum up the total and relative concentrations in the side region
!
      if (allocated(concs_depth)) deallocate(concs_depth)
      allocate(concs_depth(nelems+1))
      concs_depth=0.d0
!
!     Upper side of the slab
!
      do i=slab_edge_lo,slab_edge_lo+slab_int_depth
         concs_depth(1)=concs_depth(1)+z_dens_tot(i)
         do j=1,nelems
            concs_depth(j+1)=concs_depth(j+1)+z_dens(i,j)
         end do
      end do
!
!     Lower side of the slab
!
      do i=slab_edge_hi,slab_edge_hi-slab_int_depth,-1
         concs_depth(1)=concs_depth(1)+z_dens_tot(i)
         do j=1,nelems
            concs_depth(j+1)=concs_depth(j+1)+z_dens(i,j)
         end do
      end do


!
!     Print concentrations to file
!
      write(39,'(a,f14.8,a,f14.8,a)') "# Second element, below z=",z_min_lower1, &
                  & " A and above z=",z_min_upper1," A (%):"
      write(39,*) int_side(1,1)/(tot_side(1))*100d0
      write(39,'(a,f14.8,a,f14.8,a)') "# Second element, below z=",z_min_lower2, &
                  & " A and above z=",z_min_upper2," A (%):"
      write(39,*) int_side(2,1)/(tot_side(2))*100d0
      if (nelems .eq. 3) then
         write(39,'(a,f14.8,a,f14.8,a)') "# Third element, below z=",z_min_lower1, &
                  & " A and above z=",z_min_upper1," A (%):"
         write(39,*) int_side(1,2)/(tot_side(1))*100d0
         write(39,'(a,f14.8,a,f14.8,a)') "# Third element, below z=",z_min_lower2, &
                  & " A and above z=",z_min_upper2," A (%):"
         write(39,*) int_side(2,2)/(tot_side(2))*100d0
      end if
      write(39,*)
      write(39,*)
      write(39,*) "# In the following, the concentrations of all elements in the "
      write(39,'(a,f12.6,a)') " # region near the slab surface (defined as ",dens_depth," Angstrom"
      write(39,*) "# below the outer edges of the slab (both sides)) are given."
      do i=1,nelems
         write(39,'(a,i1,a,a,a,f12.6,a)') " * Element ",i, "(",trim(el_names(i)),&
                      & "): ",concs_depth(i+1)/concs_depth(1)*100.0," %"
         write(48,'(f15.7,a)',advance="no") concs_depth(i+1)/concs_depth(1)*100.0," "
      end do 
      write(48,*)
      close(39)
      write(*,*) "done!"
      write(*,*) "File 'surf_concs.dat' with concentrations was written."
   end if
!
!     Do evaluations if spherical densities around center of masses of elements
!     shall be calculated
!
else if (trim(dens_mode(:7)) .eq. "sphere:") then
!
!     Determine the width of a density bin: half elongation of the simulation cell
!     in the shortest dimension!
!
   zdiff = min(xlen,ylen,zlen)/2.d0
   zstep = zdiff/(nbins-1)
!
!     Extract the reference elements
!
   el_eval_list="XX"
   el_eval_num=0
   read(dens_mode(8:),*,iostat=readstat) el_eval_list
   do i=1,20
      if (el_eval_list(i) .eq. "XX") exit
      el_eval_num=el_eval_num+1
   end do

!
!     Check if the requested element is part of the system
!
   el_index=0
   do j=1,el_eval_num
      do i=1,nelems
         if (trim(el_eval_list(j)) .eq. trim(el_names(i))) then
            el_index(j)=i
         end if        
      end do
      if (el_index(j) .eq. 0) then
         write(*,*) "The element ",trim(el_eval_list(j)), &
                   & "assigned by -dens_mode is not part of the system!"
         stop
      end if    
   end do

   write(*,*) "Calculate element density distributions along the center of mass"
   write(*,'(a)',advance="no") " of the elements: "
   do i=1,el_eval_num
      write(*,'(a,a)',advance="no") el_eval_list(i)," "
   end do
   write(*,*)
!
!     Now loop through all MD steps, calculate the center of mass of the 
!      selected element atoms and calculate their surroundings
!
   z_dens=0.d0
   do i=1,nframes
      com_act=0.d0
      do j=1,el_eval_num
         do k=1,el_nums(el_index(j))
            com_act=com_act+xyz(:,sum(el_nums(:el_index(j)-1))+k,i)
         end do 
      end do
      do j=1,el_eval_num
         com_act=com_act/el_nums(el_index(j))
      end do         
!
!     Go through all atoms in the system (including the COM element), calculate
!     distances of atoms to COM (corrected by image flags) and sort them into
!     element-resolved bins
!
      counter=1
      do j=1,nelems
        
         do k=1,el_nums(j)
            diff_vec=xyz(:,counter,i)-com_act
!
!     Correct the x component
!
            do while (abs(diff_vec(1)) .gt. xlen/2.d0)
               diff_vec(1)=diff_vec(1)-sign(xlen,diff_vec(1))
            end do

!
!     Correct the y component
!
            do while (abs(diff_vec(2)) .gt. ylen/2.d0)
               diff_vec(2)=diff_vec(2)-sign(ylen,diff_vec(2))
            end do

!
!     Correct the z component
!
            do while (abs(diff_vec(3)) .gt. zlen/2.d0)
               diff_vec(3)=diff_vec(3)-sign(zlen,diff_vec(3))
            end do
!
!     Now assign each atom to one bin, depending of its distance to the COM
!
            diff_len=sqrt(dot_product(diff_vec,diff_vec))
            do l=1,nbins-1
               if ((diff_len .ge. (l-1)*zstep) .and. (diff_len .le. l*zstep)) then
                  z_dens(l,j) = z_dens(l,j) + 1.d0
               end if            
            end do
            counter = counter+1
         end do
      end do
   end do
!
!     Normalize the plot according to the volume of the sphere layers
!
   do i=1,nbins
      z_dens(i,:)=z_dens(i,:)/(nframes*4.d0/3.d0*pi*((i*zstep)**3-((i-1)*zstep)**3))
   end do
!
!     Calculate the sum density of elements
!
   do i=1,nbins-1
      z_dens_tot(i)=sum(z_dens(i,:))
   end do

!
!    Write the density profile to file
!
   open(unit=16,file="dens_elems_sphere.dat",status="replace")
   write(16,'(a)',advance="no") " # radius (A)    "
   do i=1,nelems
      write(16,'(a,a)',advance="no") el_names(i),"      "
   end do
   write(16,*) "    total  "
   do i=1,nbins-1
      write(16,*) (i-0.5d0)*zstep,z_dens(i,:),z_dens_tot(i)
   end do
   close(16)

   write(*,*) "done!"
   write(*,*) "File 'dens_elems_sphere.dat' with sphere densities was written."
else
   write(*,*) "Please give one of the valid -dens_mode commands for -dens_elems!"
   stop
end if
        
end subroutine calculate_densities


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE SURFACE_TENSION
!    Calcuate the surface tension from the pressure tensors
!    given in OUTCAR if a surface slab is investigated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine surface_tension(nframes)
use analyze_md_mod
implicit none
integer::j
integer::nframes  ! current local number of frames
real(kind=8)::pres_tensor(3,3),pres_tensor_total(3,3)  ! the pressure tensor for surface tension
real(kind=8)::pres_xx,pres_yy,pres_zz,pres_xy,pres_yz,pres_zx,tension
character(len=120)::a120
integer::openstat,readstat
integer::task_act
!
!     Calculates the averaged surface tension of the system if desired
!     Opens the OUTCAR file, reads the components of the pressure tensor,
!     averages them and finally calculates the pressure along the z axis
!
write(*,*)
write(*,*) "Calculate surface tension ..."
task_act=0
all_tasks=(nframes)
eval_stat=.false.

open(unit=18,file="OUTCAR",status="old",iostat=openstat)
if (openstat .ne. 0) then
   stop "ERROR! The file 'OUTCAR' could not been found!"
end if
tension=0.d0
pres_tensor_total=0.d0
do

   pres_tensor=0.d0

   read(18,'(a)',iostat=readstat) a120
   if (readstat .ne. 0) exit
   if (index(a120,"stress matrix after NEB project") .ne. 0) then
!
!    Every 10% of the process, give a status update
!
      task_act=task_act+1
      do j=1,10
         if (real(task_act)/real(all_tasks) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      read(18,*,iostat=readstat) pres_xx,pres_xy,pres_zx
      if (readstat .ne. 0) then
         pres_xx=0.d0
         pres_xy=0.d0
         pres_zx=0.d0
      end if
      read(18,*,iostat=readstat) pres_xy,pres_yy,pres_yz
      if (readstat .ne. 0) then
         pres_xy=0.d0
         pres_yy=0.d0
         pres_yz=0.d0
      end if
      read(18,*,iostat=readstat) pres_zx,pres_yz,pres_zz
      if (readstat .ne. 0) then
         pres_zx=0.d0
         pres_yz=0.d0
         pres_zz=0.d0
      end if
      pres_tensor(1,1)=pres_xx
      pres_tensor(2,2)=pres_yy
      pres_tensor(3,3)=pres_zz
      pres_tensor(1,2)=pres_xy
      pres_tensor(2,3)=pres_yz
      pres_tensor(1,3)=pres_zx
      pres_tensor(3,1)=pres_tensor(1,3)
      pres_tensor(3,2)=pres_tensor(2,3)
      pres_tensor(2,1)=pres_tensor(1,2)
 !     pres_tensor=-pres_tensor
   end if
   if (task_act .gt. frame_first) then
      pres_tensor_total=pres_tensor_total+pres_tensor
   end if
end do

!
!     Now calculate surface tension from formula
!
pres_tensor_total=pres_tensor_total/(nframes-frame_first)
tension=(pres_tensor_total(3,3)-0.5d0*(pres_tensor_total(1,1)+&
            & pres_tensor_total(2,2)))*1d8*1602.1766d0/box_volume
tension=0.5d0*zlen*1d-10*tension
write(*,*) " completed!"
write(*,*) "The surface tension is (N/m) or (J/m^2)",tension
open(unit=45,file="tension_act.dat",status="replace")
write(45,*) "# Surface tension of the system (N/m) or (J/m^2)"
write(45,*) tension
close(45)
write(*,*) "File 'tension_act.dat' written."
write(*,*)

end subroutine surface_tension


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE CALCULATE_DIFFUSION
!    Calcuate the self-diffusion coefficients of all elements
!    in the system from mean-displacements of atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_diffusion(nframes,xyz,xlens,ylens,zlens)
use analyze_md_mod
implicit none
integer::i,j,k,l
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
real(kind=8)::xlens(nframes),ylens(nframes),zlens(nframes)  ! the NpT unit cells
real(kind=8),allocatable::vector1(:),vector2(:),vector3(:),pos_diff(:)
real(kind=8),allocatable::vector4(:),vector5(:)
real(kind=8),allocatable::diff(:),times(:),msd_func(:,:)
real(kind=8)::diff_vec(3)
real(kind=8)::avg_diff
!
!     Calculates the diffusion coeffients for each element (self-diffusion)
!
if (diff_2d) then
   write(*,*) "Calculate the diffusion coefficients of all elements..."
else
   write(*,*) "Calculate the 2D diffusion coefficients of all elements..."
end if
allocate(pos_diff(natoms*3))
allocate(msd_func(nframes,nelems))
allocate(times(nframes))
allocate(diff(nframes))

!
!    Correct for box images: reintroduce the images for diffusion coefficients!
!
write(*,*) "Part1: Correct for box images ..."
!
!    Distinguish between NpT and NVT trajectories! (unit cell shapes)
!
do i=1,nframes-1
   do j=1,natoms
      diff_vec(:) = xyz(:,j,i+1) - xyz(:,j,i)
!
!     Correct the x component
!
      do while (abs(diff_vec(1)) .gt. xlens(i)/2d0)
         diff_vec(1)=diff_vec(1)-sign(xlens(i),diff_vec(1))
!         write(*,*) "correct x"
      end do

!
!     Correct the y component
!
      do while (abs(diff_vec(2)) .gt. ylens(i)/2d0)
         diff_vec(2)=diff_vec(2)-sign(ylens(i),diff_vec(2))
!         write(*,*) "correct y"
      end do

!
!     Correct the z component
!
      do while (abs(diff_vec(3)) .gt. zlens(i)/2d0)
         diff_vec(3)=diff_vec(3)-sign(zlens(i),diff_vec(3))
!         write(*,*) "correct z"
      end do

      xyz(:,j,i+1)=xyz(:,j,i)+diff_vec

   end do
end do

write(*,*) " ... done"

write(*,*) "Part 2: Calculate mean square displacements..."
if (nelems .eq. 1) then
   allocate(vector1(el_nums(1)*3))
else if (nelems .eq. 2) then
   allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3))
else if (nelems .eq. 3) then
   allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3),vector3(el_nums(3)*3))
else
   stop "Currently only up to three elements possible for diffusion coefficients!"
end if
do i=1,nframes
   times(i)=(i-1)*time_step
   pos_diff=0.d0
   do j=1,natoms
      if (diff_2d) then
         do k=1,2
            pos_diff((j-1)*3+k)=xyz(k,j,i)-xyz(k,j,1)
         end do
      else
         do k=1,3
            pos_diff((j-1)*3+k)=xyz(k,j,i)-xyz(k,j,1)
         end do
      end if
   end do
   if (nelems .eq. 1) then
!
!     If the collective diffusion coefficient is requested: calculate the net shift of
!      all atoms of each element
!
      if (diff_collect) then
         diff_vec=0.d0
         do k=1,natoms
            do l=1,3
               diff_vec(l)=diff_vec(l)+pos_diff((k-1)*3+l)
            end do
         end do
         msd_func(i,1)=dot_product(diff_vec,diff_vec)/natoms
      else
         msd_func(i,1)=dot_product(pos_diff,pos_diff)/natoms
      end if
   end if
   if (nelems .eq. 2) then
!
!     If the collective diffusion coefficient is requested: calculate the net shift of
!      all atoms of each element
!
      if (diff_collect) then
         diff_vec=0.d0
         do k=1,el_nums(1)
            do l=1,3
               diff_vec(l)=diff_vec(l)+pos_diff((k-1)*3+l)
            end do
         end do
         msd_func(i,1)=dot_product(diff_vec,diff_vec)/el_nums(1)
         diff_vec=0.d0
         do k=el_nums(1)+1,natoms
            do l=1,3
               diff_vec(l)=diff_vec(l)+pos_diff((k-1)*3+l)
            end do
         end do
         msd_func(i,2)=dot_product(diff_vec,diff_vec)/el_nums(2)
      else
         vector1 = pos_diff(1:el_nums(1)*3)
         msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
         vector2 = pos_diff(el_nums(1)*3+1:natoms*3)
         msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
      end if
   end if
   if (nelems .eq. 3) then
!
!     If the collective diffusion coefficient is requested: calculate the net shift of
!      all atoms of each element
!
      if (diff_collect) then
         diff_vec=0.d0
         do k=1,el_nums(1)
            do l=1,3
               diff_vec(l)=diff_vec(l)+pos_diff((k-1)*3+l)
            end do
         end do
         msd_func(i,1)=dot_product(diff_vec,diff_vec)/el_nums(1)
         diff_vec=0.d0
         do k=el_nums(1)+1,el_nums(1)+el_nums(2)
            do l=1,3
               diff_vec(l)=diff_vec(l)+pos_diff((k-1)*3+l)
            end do
         end do
         msd_func(i,2)=dot_product(diff_vec,diff_vec)/el_nums(2)
         diff_vec=0.d0
         do k=el_nums(1)+el_nums(2)+1,natoms
            do l=1,3
               diff_vec(l)=diff_vec(l)+pos_diff((k-1)*3+l)
            end do
         end do
         msd_func(i,3)=dot_product(diff_vec,diff_vec)/el_nums(3)
      else
         vector1 = pos_diff(1:el_nums(1)*3)
         msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
         vector2 = pos_diff(el_nums(1)*3+1:(el_nums(1)+el_nums(2))*3)
         msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
         vector3 = pos_diff((el_nums(1)+el_nums(2))*3+1:natoms*3)
         msd_func(i,3)=dot_product(vector3,vector3)/el_nums(3)
      end if
   end if
   if (nelems .eq. 4) then
      vector1 = pos_diff(1:el_nums(1)*3)
      msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
      vector2 = pos_diff(el_nums(1)*3+1:(el_nums(1)+el_nums(2))*3)
      msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
      vector3 = pos_diff(sum(el_nums(1:2))*3+1:sum(el_nums(1:3))*3)
      msd_func(i,3)=dot_product(vector3,vector3)/el_nums(3)
      vector4 = pos_diff(sum(el_nums(1:3))*3+1:natoms*3)
      msd_func(i,4)=dot_product(vector4,vector4)/el_nums(4)
   end if
   if (nelems .eq. 5) then
      vector1 = pos_diff(1:el_nums(1)*3)
      msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
      vector2 = pos_diff(el_nums(1)*3+1:(el_nums(1)+el_nums(2))*3)
      msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
      vector3 = pos_diff(sum(el_nums(1:2))*3+1:sum(el_nums(1:3))*3)
      msd_func(i,3)=dot_product(vector3,vector3)/el_nums(3)
      vector4 = pos_diff(sum(el_nums(1:3))*3+1:sum(el_nums(1:4))*3)
      msd_func(i,4)=dot_product(vector4,vector4)/el_nums(4)
      vector5 = pos_diff(sum(el_nums(1:4))*3+1:natoms*3)
      msd_func(i,5)=dot_product(vector5,vector5)/el_nums(5)
   end if
end do
write(*,*) " ... done"

open(unit=17,file="msd_plot.dat",status="replace")
write(17,*) "# time (fs),   MSD (A^2)"
do i=1,nframes
   if (nelems .eq. 1) then
      write(17,*) times(i),msd_func(i,1)
   end if
   if (nelems .eq. 2) then
      write(17,*) times(i),msd_func(i,1),msd_func(i,2)
   end if
   if (nelems .eq. 3) then
      write(17,*) times(i),msd_func(i,1),msd_func(i,2),msd_func(i,3)
   end if        
   if (nelems .eq. 4) then
      write(17,*) times(i),msd_func(i,1),msd_func(i,2),msd_func(i,3), &
                        & msd_func(i,4)
   end if
   if (nelems .eq. 5) then
      write(17,*) times(i),msd_func(i,1),msd_func(i,2),msd_func(i,3), &
                        & msd_func(i,4),msd_func(i,5)
   end if
   
end do

close(17)
!
!     Calculate the diffusion coefficient by calculating the MSD
!     based on the last time step!
!
do k=1,nelems
   do i=1,nframes
      if (diff_2d) then
         diff(i)=msd_func(i,k)/(4.d0*times(i))
      else
         diff(i)=msd_func(i,k)/(6.d0*times(i))
      end if
   end do
   diff = diff*(1E-10)**2/(1E-15)
   avg_diff=diff(nframes)

   write(*,*) "calculated_diffusion coefficient, ",trim(el_names(k))," (m^2/s):",avg_diff
   open(unit=19,file="diff_const_MSD_"//trim(el_names(k))//".dat",status="replace")
   write(19,*) avg_diff,"m^2/s"
   close(19)
end do

write(*,*) "Diffusion coefficient calculation finished!"
write(*,*) "Coefficients written to diff_const_MSD_*.dat"

end subroutine calculate_diffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE PICK_STRUCTURES
!    Pick exemplary structures with certain properties, e.g., 
!    position of an atom of an element, from a number of slices
!    of a trajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pick_structures(nframes,xyz)
use analyze_md_mod
implicit none
integer::i,j,k,l,r
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
integer::slice_size,frame_round_first,frame_round_last
integer::atom_first,atom_last
real(kind=8)::slice_step
character(len=60)::filename,roundname
logical,allocatable::struc_used(:)
real(kind=8)::zlo,zhi,zdiff
!
!     Additionally, write POSCAR files for Core Level energy calculations
!     of the chosen active atom species at different positions if wanted
!


zlo = 0.d0 ! minval(xyz(3,:,:))
zhi = zlen !maxval(xyz(3,:,:))
zdiff=zhi-zlo

slice_step = zdiff/(atom_slices)

if (cls_element .ne. "XX") then
   write(*,*) "Write POSCARs for DFT example calculatons..."
   atom_first=sum(el_nums(1:cls_elem_ind-1))+1
   atom_last=sum(el_nums(1:cls_elem_ind))
   call system("mkdir dft_input")
   slice_size=int(nframes/cls_rounds)
   rounds: do r=1,cls_rounds
      if (r .le. 9) then
          write(roundname,'(a,i1)') "round",r
      else if (r .le. 99) then
          write(roundname,'(a,i2)') "round",r
      else
          write(roundname,'(a,i3)') "round",r
      end if
      call system("mkdir dft_input/" //roundname)

      open(unit=18,file="dft_input/" // trim(roundname) // "/active_examples.xyz",status="replace")
      open(unit=19,file="dft_input/" // trim(roundname) // "/active_examples.XDATCAR",status="replace")
!
!     Write header of XDATCAR file
!
      write(19,*) "Example positions of active atoms with different z-coordinates"
      write(19,*) factor
      write(19,*) xlen,0.0,0.0
      write(19,*) 0.0,ylen,0.0
      write(19,*) 0.0,0.0,zlen
      do i=1,nelems
         write(19,'(a,a)',advance="no")"  ", el_names(i)
      end do
      write(19,*)
      write(19,*) el_nums



      frame_round_first=1+(r-1)*slice_size
      frame_round_last=1+r*slice_size
      if (allocated(struc_used)) deallocate(struc_used)
      allocate(struc_used(frame_round_last-frame_round_first))
      struc_used=.false.
      write(*,*) "DFT round ",r,": From frame ", frame_round_first," to ",frame_round_last
      slices: do i=1,atom_slices
         if (i .le. 9) then
            write(filename,'(a,a,a,i1)') "dft_input/",trim(roundname),"/POSCAR_slice",i
         else if (i .le. 99) then
            write(filename,'(a,a,a,i2)') "dft_input/",trim(roundname),"/POSCAR_slice",i
         else
            write(filename,'(a,a,a,i3)') "dft_input/",trim(roundname),"/POSCAR_slice",i
         end if

         frames: do j=frame_round_first,frame_round_last
            if (struc_used(j-frame_round_first+1)) cycle frames
            do k=atom_first,atom_last
               if ((xyz(3,k,j) .ge. zlo+(i-1)*slice_step) .and. (xyz(3,k,j) &
                           & .le. zlo+i*slice_step)) then
                  struc_used(j-frame_round_first+1)=.true.
                  write(18,*) natoms
                  write(18,'(a,i4,a,i4,a,f13.6,a,f13.6)') "Layer ",i," of ", &
                           & atom_slices,": z=",zlo+(i-1)*slice_step," to ", &
                           & zlo+i*slice_step
                  write(19,'(a,i4,a,i4,a,f13.6,a,f13.6)') "Direct Layer ",i," of ", &
                           & atom_slices,": z=",zlo+(i-1)*slice_step," to ", &
                           & zlo+i*slice_step
!
!     Write header of POSCAR file for Core Level shift
!
                  open(unit=20,file=filename)
                  write(20,*) "Active atom at z = ",xyz(3,k,j),", input for DFT calculation (e.g., CLS)."
                  write(20,*) factor
                  write(20,*) xlen,0.0,0.0
                  write(20,*) 0.0,ylen,0.0
                  write(20,*) 0.0,0.0,zlen
                  do l=1,nelems
                     write(20,'(a,a)',advance="no")"  ", el_names(l)
                  end do
                  write(20,*) "   ",el_names(cls_elem_ind)
                  do l=1,nelems-1
                     write(20,'(i10)',advance="no") el_nums(l)
                  end do
                  write(20,'(i10)',advance="no") el_nums(nelems)-1
                  write(20,*) 1
                  write(20,*) "Direct"
                  do l=1,natoms
                     write(18,*) at_names(l), xyz(:,l,j)
                     write(19,*) xyz(1,l,j)/xlen,xyz(2,l,j)/ylen,xyz(3,l,j)/zlen
                     if (l .ne. k) then
                        write(20,*) xyz(1,l,j)/xlen,xyz(2,l,j)/ylen,xyz(3,l,j)/zlen
                     end if
                  end do
                  write(20,*) xyz(1,k,j)/xlen,xyz(2,k,j)/ylen,xyz(3,k,j)/zlen
                  close(20)
                  cycle slices
               end if
            end do
         end do frames
  !       write(*,*) "Warning: No Active atom found in slice No.",i
      end do slices
      close(18)
      close(19)
   end do rounds
   write(*,*) " completed!"
   write(*,*)
end if
if (cls_element .ne. "XX") then
   write(*,*) "POSCAR files for DFT single points written to folder 'dft_input/round...'."
   write(*,*) "Files 'active_examples.xyz/.XDATCAR' with active atom positions written "
   write(*,*) "  to dft_input/round...."
end if
end subroutine pick_structures


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE SPACIAL_DENSITY
!    Calculate the averaged three-dimensional elemental 
!    density of all elements and write them into a 
!    Gaussian cube file, usable for visualization with VMD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
subroutine spacial_density(nframes,xyz,xlens,ylens,zlens)
use analyze_md_mod
implicit none 
integer::i,j,k,l
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
real(kind=8)::xlens(nframes),ylens(nframes),zlens(nframes) ! unit cells
integer::i1,j1,k1,i_new,j_new,k_new
integer::gridx_act,gridy_act,gridz_act
integer::gridx,gridy,gridz
integer::stepx,stepy,stepz
real(kind=8)::distance
integer::inc,elnumber
real(kind=8)::a2bohr
character(len=30)::cube_name
real(kind=8),allocatable::dens_3d(:,:,:,:)

a2bohr=1.8897259886d0

write(*,*) "Calculate the spatially resolved element densities!"
!
!     Determine number grid points per dimension
!

gridx=100
gridy=100
gridz=100
allocate(dens_3d(nelems,gridx,gridy,gridz))
dens_3d=0.d0

!
!    For 3D element densities, determine number of steps around current position
!    One Angstrom around
!
stepx=int(gridx/xlens(nframes))
stepy=int(gridy/ylens(nframes))
stepz=int(gridz/zlens(nframes))
!
!    Convert back to direct coordinates 
!
do i=1,nframes
   do j=1,natoms
      xyz(1,j,i) = xyz(1,j,i)/xlens(i)
      xyz(2,j,i) = xyz(2,j,i)/ylens(i)
      xyz(3,j,i) = xyz(3,j,i)/zlens(i)
   end do
end do

!
!     Remove corrections of box images, project all atoms into central unit cell
!
do i=1,nframes
   do j=1,natoms
!
!     Correct the x component
!
      do while(xyz(1,j,i) .gt. 1.d0)
         xyz(1,j,i)=xyz(1,j,i)-1.d0
      end do
      do while(xyz(1,j,i) .lt. 0.d0)
         xyz(1,j,i)=xyz(1,j,i)+1.d0
      end do
!
!     Correct the y component
!
      do while(xyz(2,j,i) .gt. 1.d0)
         xyz(2,j,i)=xyz(2,j,i)-1.d0
      end do
      do while(xyz(2,j,i) .lt. 0.d0)
         xyz(2,j,i)=xyz(2,j,i)+1.d0
      end do
!
!     Correct the z component
!
      do while(xyz(3,j,i) .gt. 1.d0)
         xyz(3,j,i)=xyz(3,j,i)-1.d0
      end do
      do while(xyz(3,j,i) .lt. 0.d0)
         xyz(3,j,i)=xyz(3,j,i)+1.d0
      end do
   end do
end do
!     Now assign all atoms in each frame to one of the 3D bins
!
do i=frame_first,nframes
   inc=0
   do j=1,nelems
      do k=1,el_nums(j)
         inc=inc+1
         gridx_act=ceiling(xyz(1,inc,i)*gridx)
         gridy_act=ceiling(xyz(2,inc,i)*gridy)
         gridz_act=ceiling(xyz(3,inc,i)*gridz)-1
         do i1=-stepx,stepx
            do j1=-stepy,stepy
               do k1=-stepz,stepz
                  distance=i1**2+j1**2+k1**2
                  i_new=gridx_act+i1
                  if (i_new .gt. gridx) then
                     i_new=i_new-gridx
                  else if (i_new .lt. 1) then
                     i_new=i_new+gridx
                  end if
                  j_new=gridy_act+j1
                  if (j_new .gt. gridy) then
                     j_new=j_new-gridy
                  else if (j_new .lt. 1) then
                     j_new=j_new+gridy
                  end if
                  k_new=gridz_act+k1
                  if (k_new .gt. gridz) then
                     k_new=k_new-gridz
                  else if (k_new .lt. 1) then
                     k_new=k_new+gridz
                  end if

                  dens_3d(j,i_new,j_new,k_new)= &
                          & dens_3d(j,i_new,j_new,k_new)+exp(-(distance)/2.d0)
               end do
            end do
         end do
      end do
   end do
end do
!
!     Write the header of the cube file
!
do i=1,nelems
   write(cube_name,'(a,a,a)') "density_",trim(el_names(i)),".cube"
   write(*,*) "Write density of ",trim(el_names(i))," to file ",cube_name
   open(unit=45,file=cube_name,status="replace")
   write(45,'(a)') "utils4VASP CUBE FILE"
   write(45,'(a,a)') "Contains spatially resolved density of element:",el_names(i)
   write(45,*) natoms,0d0,0.d0,0d0
!   if (.not. npt_traj) then
      write(45,*) gridx,xlen/gridx*a2bohr,0.d0,0.d0
      write(45,*) gridy,0.d0,ylen/gridy*a2bohr,0.d0
      write(45,*) gridz,0.d0,0.d0,zlen/gridz*a2bohr
!   else
!      write(45,*) gridx,a_lens(1)/gridx*a2bohr,0.d0,0.d0
!      write(45,*) gridy,0.d0,b_lens(1)/gridy*a2bohr,0.d0
!      write(45,*) gridz,0.d0,0.d0,c_lens(1)/gridz*a2bohr
!   end if
   inc=0
   do j=1,nelems
      do k=1,el_nums(j)
         inc=inc+1
         call elem(el_names(j),elnumber)
         write(45,*) elnumber,real(elnumber),xyz(1,inc,nframes)*xlens(1)*a2bohr,xyz(2,inc,nframes)* &
                        & ylens(1)*a2bohr,xyz(3,inc,nframes)*zlens(1)*a2bohr
      end do
   end do
   write(45,*) 1,48
!
!    Write the volumetric data in the cube file
!
   do j=1,gridx
      do k=1,gridy
         do l=1,gridz
            write(45,'(f20.10)',advance="no") dens_3d(i,j,k,l)
            if (modulo(l,6) .eq. 5) then
               write(45,*) " "
            end if
         end do
         write(45,*) " "
      end do
   end do
   close(45)
end do

end subroutine spacial_density


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE STABILITY
!    Determine if all bonds remain stable during the part of 
!    the MD. First determine all active bonds in the first 
!    frame, then check in the following frames if any of those
!    bonds got more than stable_crit longer than initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stability(nframes,xyz,xlens,ylens,zlens)
use analyze_md_mod
implicit none
integer::i,j,k,l
integer::nframes  ! current local number of frames
real(kind=8)::xyz(3,natoms,nframes)  ! local part of the trajectory
real(kind=8)::xlens(nframes),ylens(nframes),zlens(nframes) ! unit cells
real(kind=8)::atom_radii(natoms)  ! covalent radii for all atoms 
real(kind=8)::pos1(3),pos2(3),diff_vec(3),dist
integer,allocatable::bond_list(:,:)  ! list with bonds
integer::bond_num  ! number of bonds
real(kind=8),allocatable::bond_lengths(:)   ! The bond lengths at the beginning


write(*,*) "Check the stability of the systems (bond breakings)!"

!
!    Determine covalent radii of all atoms in the system
!

do i=1,natoms
   if (at_names(i) .eq. "H") then
      atom_radii(i)=0.32d0
   else if (at_names(i) .eq. "He") then
      atom_radii(i)=0.46d0
   else if (at_names(i) .eq. "Li") then
      atom_radii(i)=1.33d0
   else if (at_names(i) .eq. "Be") then
      atom_radii(i)=1.02d0
   else if (at_names(i) .eq. "B") then
      atom_radii(i)=0.85d0
   else if (at_names(i) .eq. "C") then
      atom_radii(i)=0.75d0
   else if (at_names(i) .eq. "N") then
      atom_radii(i)=0.71d0
   else if (at_names(i) .eq. "O") then
      atom_radii(i)=0.63d0
   else if (at_names(i) .eq. "F") then
      atom_radii(i)=0.64d0
   else if (at_names(i) .eq. "Ne") then
      atom_radii(i)=0.67d0
   else if (at_names(i) .eq. "Na") then
      atom_radii(i)=1.55d0
   else if (at_names(i) .eq. "Mg") then
      atom_radii(i)=1.39d0
   else if (at_names(i) .eq. "Al") then
      atom_radii(i)=1.26d0
   else if (at_names(i) .eq. "Si") then
      atom_radii(i)=1.16d0
   else if (at_names(i) .eq. "P") then
      atom_radii(i)=1.11d0
   else if (at_names(i) .eq. "S") then
      atom_radii(i)=1.03d0
   else if (at_names(i) .eq. "Cl") then
      atom_radii(i)=0.99d0
   else if (at_names(i) .eq. "Ar") then
      atom_radii(i)=0.96d0
   else if (at_names(i) .eq. "Sc") then
      atom_radii(i)=1.48d0
   else if (at_names(i) .eq. "Ti") then
      atom_radii(i)=1.36d0
   else if (at_names(i) .eq. "V") then
     atom_radii(i)=1.34d0
   else if (at_names(i) .eq. "Cr") then
      atom_radii(i)=1.22d0
   else if (at_names(i) .eq. "Mn") then
      atom_radii(i)=1.19d0
   else if (at_names(i) .eq. "Fe") then
      atom_radii(i)=1.160
   else if (at_names(i) .eq. "Co") then
      atom_radii(i)=1.11d0
   else if (at_names(i) .eq. "Ni") then
      atom_radii(i)=1.10d0
   else if (at_names(i) .eq. "Cu") then
      atom_radii(i)=1.12d0
   else if (at_names(i) .eq. "Zn") then
      atom_radii(i)=1.18d0
   else if (at_names(i) .eq. "Ga") then
      atom_radii(i)=1.24d0
   else if (at_names(i) .eq. "Ge") then
      atom_radii(i)=1.21d0
   else if (at_names(i) .eq. "As") then
      atom_radii(i)=1.21d0
   else if (at_names(i) .eq. "Se") then
      atom_radii(i)=1.16d0
   else if (at_names(i) .eq. "Br") then
      atom_radii(i)=1.14d0
   else if (at_names(i) .eq. "Kr") then
      atom_radii(i)=1.17d0
   else if (at_names(i) .eq. "Mo") then
      atom_radii(i)=1.38d0
   else if (at_names(i) .eq. "Tc") then
      atom_radii(i)=1.28d0
   else if (at_names(i) .eq. "Ru") then
      atom_radii(i)=1.25d0
   else if (at_names(i) .eq. "Rh") then
      atom_radii(i)=1.25d0
   else if (at_names(i) .eq. "Pd") then
      atom_radii(i)=1.20d0
   else if (at_names(i) .eq. "Ag") then
      atom_radii(i)=1.28d0
   else if (at_names(i) .eq. "Cd") then
      atom_radii(i)=1.36d0
   else if (at_names(i) .eq. "In") then
      atom_radii(i)=1.42d0
   else if (at_names(i) .eq. "Sn") then
      atom_radii(i)=1.40d0
   else if (at_names(i) .eq. "Sb") then
      atom_radii(i)=1.40d0
   else if (at_names(i) .eq. "Te") then
      atom_radii(i)=1.36d0
   else if (at_names(i) .eq. "I") then
      atom_radii(i)=1.33d0
   else if (at_names(i) .eq. "Xe") then
      atom_radii(i)=1.31d0
   else if (at_names(i) .eq. "La") then
      atom_radii(i)=1.80d0
   else if (at_names(i) .eq. "Ce") then
      atom_radii(i)=1.63d0
   else if (at_names(i) .eq. "Pr") then
      atom_radii(i)=1.76d0
   else if (at_names(i) .eq. "Nd") then
      atom_radii(i)=1.74d0
   else if (at_names(i) .eq. "Pm") then
      atom_radii(i)=1.73d0
   else if (at_names(i) .eq. "Sm") then
      atom_radii(i)=1.72d0
   else if (at_names(i) .eq. "Eu") then
      atom_radii(i)=1.68d0
   else if (at_names(i) .eq. "Gd") then
      atom_radii(i)=1.69d0
   else if (at_names(i) .eq. "Tb") then
      atom_radii(i)=1.68d0
   else if (at_names(i) .eq. "Dy") then
      atom_radii(i)=1.67d0
   else if (at_names(i) .eq. "Ho") then
      atom_radii(i)=1.66d0
   else if (at_names(i) .eq. "Er") then
      atom_radii(i)=1.65d0
   else if (at_names(i) .eq. "Tm") then
      atom_radii(i)=1.64d0
   else if (at_names(i) .eq. "Yb") then
      atom_radii(i)=1.70d0
   else if (at_names(i) .eq. "Lu") then
      atom_radii(i)=1.62d0
   else if (at_names(i) .eq. "Hf") then
      atom_radii(i)=1.52d0
   else if (at_names(i) .eq. "Ta") then
      atom_radii(i)=1.46d0
   else if (at_names(i) .eq. "W") then
      atom_radii(i)=1.37d0
   else if (at_names(i) .eq. "Re") then
      atom_radii(i)=1.31d0
   else if (at_names(i) .eq. "Os") then
      atom_radii(i)=1.29d0
   else if (at_names(i) .eq. "Ir") then
      atom_radii(i)=1.22d0
   else if (at_names(i) .eq. "Pt") then
      atom_radii(i)=1.23d0
   else if (at_names(i) .eq. "Au") then
      atom_radii(i)=1.24d0
   else if (at_names(i) .eq. "Hg") then
      atom_radii(i)=1.33d0
   else if (at_names(i) .eq. "Tl") then
      atom_radii(i)=1.44d0
   else if (at_names(i) .eq. "Pb") then
      atom_radii(i)=1.44d0
   else if (at_names(i) .eq. "Bi") then
      atom_radii(i)=1.51d0
   else if (at_names(i) .eq. "Po") then
      atom_radii(i)=1.45d0
   else if (at_names(i) .eq. "At") then
      atom_radii(i)=1.47d0
   else if (at_names(i) .eq. "Rn") then
      atom_radii(i)=1.42d0
   else if (at_names(i) .eq. "Ac") then
      atom_radii(i)=1.86d0
   else if (at_names(i) .eq. "Th") then
      atom_radii(i)=1.75d0
   else if (at_names(i) .eq. "Pa") then
      atom_radii(i)=1.69d0
   else if (at_names(i) .eq. "U") then
      atom_radii(i)=1.70d0
   else if (at_names(i) .eq. "Np") then
      atom_radii(i)=1.71d0
   else if (at_names(i) .eq. "Pu") then
      atom_radii(i)=1.72d0
   else if (at_names(i) .eq. "Am") then
      atom_radii(i)=1.66d0
   else if (at_names(i) .eq. "Cm") then
      atom_radii(i)=1.66d0
   end if
end do
!
!    Allocate the bond list with a maximum number of bonds
!
allocate(bond_list(2,int((natoms*natoms-natoms)/2)))
allocate(bond_lengths(int((natoms*natoms-natoms)/2)))
!
!    Determine all bonds in the first MD frame
!
bond_num=0
do i=1,natoms
   do j=i+1,natoms

      pos1 = xyz(:,i,1)
      pos2 = xyz(:,j,1)

      diff_vec=pos1-pos2
!
!     Correct the x component
!
      do while (abs(diff_vec(1)) .gt. xlens(1)/2.d0)
         diff_vec(1)=diff_vec(1)-sign(xlens(1),diff_vec(1))
      end do

!
!     Correct the y component
!
      do while (abs(diff_vec(2)) .gt. ylens(1)/2.d0)
         diff_vec(2)=diff_vec(2)-sign(ylens(1),diff_vec(2))
      end do

!
!     Correct the z component
!
      do while (abs(diff_vec(3)) .gt. zlens(1)/2.d0)
         diff_vec(3)=diff_vec(3)-sign(zlens(1),diff_vec(3))
      end do

      dist = sqrt((diff_vec(1))**2 + &
            & (diff_vec(2))**2 + (diff_vec(3))**2)
!
!     If the distance is smaller than the sum of radii times a factor,
!      mark it as a bond
!
      if (dist .lt. (atom_radii(i)+atom_radii(j))*1.4d0) then
         bond_num=bond_num+1
         bond_list(1,bond_num)=i
         bond_list(2,bond_num)=j
         bond_lengths(bond_num)=dist
      end if        

   end do
end do   

!
!   Now loop through all remaining MD frames and check if a 
!      bond is too elongated
!
do i=1,nframes
   do j=1,bond_num
      pos1=xyz(:,bond_list(1,j),i)
      pos2=xyz(:,bond_list(2,j),i)
      diff_vec=pos1-pos2
!
!     Correct the x component
!
      do while (abs(diff_vec(1)) .gt. xlens(1)/2.d0)
         diff_vec(1)=diff_vec(1)-sign(xlens(1),diff_vec(1))
      end do

!
!     Correct the y component
!
      do while (abs(diff_vec(2)) .gt. ylens(1)/2.d0)
         diff_vec(2)=diff_vec(2)-sign(ylens(1),diff_vec(2))
      end do

!
!     Correct the z component
!
      do while (abs(diff_vec(3)) .gt. zlens(1)/2.d0)
         diff_vec(3)=diff_vec(3)-sign(zlens(1),diff_vec(3))
      end do

      dist = sqrt((diff_vec(1))**2 + &
            & (diff_vec(2))**2 + (diff_vec(3))**2)
!
!     If the distance is more than bond_tol longer than the initial
!       bond length, mark the bond as broken and exit the analysis!
!
      if (dist .gt. (bond_lengths(j)*bond_tol)) then
         write(*,*) "A broken bond has been detected!"
         write(*,'(a,i5,a,a,a,i5,a,a,a)') " Bond between atoms ", &
               & bond_list(1,j),"(", &
               &  trim(at_names(bond_list(1,j))),") and ",bond_list(2,j), &
               & "(",trim(at_names(bond_list(2,j))),")"
         write(*,'(a,f12.6,a,f12.6,a)') " Initial bond length:",bond_lengths(j), &
               & " Ang., current bond length:",dist," Ang."
         write(*,'(a,i9)') " This happened at MD step",i
         write(*,*) "The analysis will be aborted..."
!
!     Write the initial frame and the final frame when the breaking is detected 
!      to a xyz file with two frames
!
         open(unit=78,file="bond_break.xyz",status="replace")
         write(78,*) natoms
         write(78,*) "Initial structure of the trajectory."
         do k=1,natoms
            write(78,*) at_names(k)," ",xyz(:,k,1)
         end do
         write(78,*) natoms
         write(78,'(a,i8,a,i6,a,i6)') "Frame with broken bond (No.",i, &
                    &"), bond ",bond_list(1,j),"-",bond_list(2,j)
         do k=1,natoms
            write(78,*) at_names(k)," ",xyz(:,k,i)
         end do
         close(78)
         write(*,*) "File 'bond_break.xyz' with first frame and broken frame written."
         goto 48
      end if        

   end do
end do
write(*,*) "All bonds were stable during the trajectory."
48 continue

end subroutine stability



subroutine replace_text (s,text,rep,outs)
CHARACTER(*)        :: s,text,rep
CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
INTEGER             :: i, nt, nr

outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
DO
   i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
END DO

return 
end subroutine replace_text


!
!     subroutine elem: read in character with element symbol and
!       give out the number
!
!     part of QMDFF
!
subroutine elem(key1, nat)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CHARACTER(len=*)::KEY1
CHARACTER(len=2)::ELEMNT(107),E

DATA ELEMNT/'h ','he', &
 & 'li','be','b ','c ','n ','o ','f ','ne', &
 & 'na','mg','al','si','p ','s ','cl','ar', &
 & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
 & 'zn','ga','ge','as','se','br','kr', &
 & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
 & 'cd','in','sn','sb','te','i ','xe', &
 & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
 & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
 & 'au','hg','tl','pb','bi','po','at','rn', &
 & 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx', &
 & 'fm','md','cb','xx','xx','xx','xx','xx'/

nat=0
e='  '
do i=1,len(key1)
   if (key1(i:i).ne.' ') L=i
end do
k=1
DO J=1,L
   if (k.gt.2) exit
   N=ICHAR(key1(J:J))
   if (n.ge.ichar('A') .and. n.le.ichar('Z') ) then
      e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
      k=k+1
   end if
   if (n.ge.ichar('a') .and. n.le.ichar('z') ) then
      e(k:k)=key1(j:j)
      k=k+1
   end if
end do

DO I=1,107
   if (e.eq.elemnt(i)) then
      NAT=I
      RETURN
   END IF
END DO

return
end subroutine elem

