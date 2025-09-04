!
!    mlff_select: selection of an effective set of 
!      local reference configurations (basis sets) for 
!      a VASP machine-learning force field from a given
!      ML_AB file. The procedure is similar to the VASP 
!      command ML_MODE = select, but is faster and has more
!      options to be adjusted.
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!
program mlff_select
implicit none
integer::i,j,k,l,m   ! loop indices
integer::inc,inc2,inc3,inc4  ! incremented index
integer::inc_atom  ! element-global atom position increment (basis set)
integer::inc_act   ! the current increment
integer::inc_remain  ! remainder index for basis functions
integer::mlab_num   ! number of ML_AB files to be read in
integer,allocatable::incs(:)  ! moredimensional increment
integer::readstat,ret
integer::basis_mode   ! if all elements have same basis or each one enother
integer::nbasis(50)  ! desired number of local refs. in new ML_AB
integer::nbasis_tmp(50)  ! for resorting of nbasis array, if needed
integer::neigh_min_bas(50)  ! minimum number of confs. per neighbor bin
integer::nbas_grad(50)  ! number of basis functions for max.  gradient comps
integer::neigh_classes ! number of neighborhood subdivision classes (calculated)
integer::neigh_maxsize ! maximum number of basis functions per neigh_class (input)
integer::grad_pre(50)    ! number of basis functions for large gradnorms
integer::conf_num  ! number of configurations (full structures)
integer::act_conf   ! the current configuration number
integer::natoms_max  ! maximum number of atoms per configuration
integer::nat_type_max  ! maximum number of atoms per element 
integer::natoms_sum  ! sum of all atoms within all configurations
integer::nelems  ! current number of elements
integer::nelems_all  ! total number of different elements
integer::nelems_local(20)   ! local number of different elements 
integer::nelems_glob  ! number of different elements in all ML_ABs
integer::ind_act ! current element index
integer::max_environ  ! maximum number of atoms in environment
integer::max_around  ! maximum actual number of atoms around
integer::ngrid  ! number of grid points for RDF calculation
integer::spots_avail  ! number of basis functions for final k-means
integer::remain_class  ! index of the remainder class for basis functions
integer::nconfs_out   ! number of configurations to be written out
integer,allocatable::natoms(:)  ! atom numbers of configurations
integer,allocatable::nelems_confs(:)  ! number of elements in all configurations
real(kind=8),allocatable::xyz(:,:,:) ! the geometries
real(kind=8),allocatable::xyz_dir(:,:,:) ! the geometries (direct coords)
real(kind=8),allocatable::cells(:,:,:) ! the unit cells
real(kind=8),allocatable::ctifors(:)  ! the CTIFOR values for all configurations
real(kind=8),allocatable::stress(:,:)  ! stress tensors for all configurations
real(kind=8),allocatable::angles(:)  ! list of angles for current environment
real(kind=8)::time1,time2  ! Time measurement of program execution
real(kind=8)::grad_min,grad_max  ! minimum and maximum gradient norm (global)
real(kind=8)::grad_step  ! distance between two gradient bins
real(kind=8)::normfac   ! normalizaion factor for RDFs and ADFs
real(kind=8)::rdf_adf   ! relative weighting of RDF and ADFs in clustering
real(kind=8)::rdum   ! dummy real number
integer,allocatable::inds(:,:)  ! the element indices (core charges)
integer,allocatable::natoms_el(:)  ! total number of atoms per element
integer::ind_list_loc(50,20)   ! local list of different element indices
integer::ind_list(50)    ! global list of different element indices
character(len=2),allocatable::el_list_confs(:,:)  ! List of elements in all confs.
character(len=2),allocatable::el_list_tmp(:)  ! List of elements in actual conf
character(len=2)::el_list_loc(50,20)  ! local list of elements
character(len=2)::el_list_glob(50)  ! global list of elements
character(len=2)::el_list_bas(50)  ! list of elements for basis init.
character(len=2)::el_act   ! the current element symbol for printing
character(len=300)::arg,adum ! command line argument
character(len=20)::bas_scale   ! linear or sqare root scaling of basis classes
character(len=50)::mlab_list(20)  ! list of ML_AB files to be read in
character(len=50)::basis_read(50)  ! read in for element-specific basis functions
real(kind=8)::at_mass_glob(50)  ! Global list with atomic masses
real(kind=8)::at_mass(50,20)  ! local list with atomic masses
real(kind=8)::train_div   ! desired trainset diversity (0: none, 1: full)
real(kind=8)::grad_frac   ! desired fraction of basis functions for large gradnorm
real(kind=8),allocatable::energies(:)  ! the energies 
real(kind=8),allocatable::grads(:,:,:)  ! the gradient vectors
real(kind=8),allocatable::environ(:,:,:)  ! the local environment of each atom
real(kind=8),allocatable::dist_list(:)  ! list of distances in environment
real(kind=8),allocatable::rdf_all(:,:)  ! radial distribution functions for atoms
real(kind=8),allocatable::adf_all(:,:)  ! angular distribution functions for atoms
real(kind=8),allocatable::atsel_grad_val(:,:)  ! gradnorms of largest gradnorms 
real(kind=8),allocatable::rdf_tmp(:,:)  ! temporary RDF array
real(kind=8),allocatable::mat_overlap(:,:)  ! overlap matrix for RDF combinations
real(kind=8)::cell_inv(3,3)  ! inverted coordinate unit cell
real(kind=8)::dist_act   ! scalar distance between two atoms
real(kind=8)::cutoff    ! the distance cutoff during the ML-FF learning
real(kind=8)::dist_vec(3)  ! distance vector between two atoms
real(kind=8)::dx   ! distance between two gridpoints (RDF)
real(kind=8)::da   ! distance between two gridpoints (ADF)
real(kind=8)::alpha_r   ! exponent for RDF Gaussian functions
real(kind=8)::alpha_a   ! exponent for ADF Gaussian functions
real(kind=8)::x_act   ! current x-value for RDF buildup
real(kind=8)::rdf_overlap  ! the current RDF overlap integral
real(kind=8)::adf_overlap  ! the current ADF overlap integral
real(kind=8)::rdf_weight  ! relative weight of RDF overlap 
real(kind=8)::adf_weight  ! relative weight of ADF overlap
real(kind=8)::pi  ! the Pi
integer,allocatable::el_nums_confs(:,:)  ! element numbers in all configurations
integer,allocatable::el_nums_tmp(:)  ! element numbers in current conf.
integer,allocatable::confnum_all(:)  ! numbers of configurations 
integer,allocatable::mlab_conf(:) !  to which ML_AB file the configuration belongs
integer,allocatable::mlab_atom(:)  ! to which ML_AB file the atom belongs
integer,allocatable::nat_all(:) ! number of atom in configurations
integer,allocatable::ind_all(:)  ! element indices of all atoms
integer,allocatable::sort_map(:)  ! resorting of atoms within structure
integer,allocatable::sort_map2(:)  ! second resorting array (inverted)
integer,allocatable::num_around(:,:) ! number of atoms around the atom
integer,allocatable::ind_env(:,:)  ! core charges of atoms in environments
integer,allocatable::neigh_global(:,:)  ! list with numbers of neighbors (global)
integer,allocatable::neigh_bas(:,:)  ! number of basis functions per environment
integer,allocatable::grad_histo(:)   ! histogram with gradient norm ranges
integer,allocatable::final_choice(:,:)  ! the final selection of basis functions
integer,allocatable::trans_final(:)  ! for local reordering of final selection
integer,allocatable::atsel_grad(:,:)  ! the atom indices of the largest gradnorms
integer,allocatable::neighnum_global(:)  ! sum of all neighbors for all atoms
integer,allocatable::natoms_neigh(:) ! number of atoms per neighborhood number
integer,allocatable::nbasis_neigh(:)  ! number of basis functions per neighborhood number
integer,allocatable::neigh_div1(:,:),neigh_div2(:,:)  ! the neighborhood diversities
integer,allocatable::basis_sizes(:)  ! sizes of all basis classes (per element)
integer,allocatable::bas_neighs(:)  ! total number of neighbors per class
integer,allocatable::basis_classes(:,:)  ! list of atoms in each basis class
integer,allocatable::k_number(:)   ! list of k-means spots for basis classes
integer,allocatable::smallest(:)  ! atom list with smallest mutual overlaps
integer,allocatable::hierarchy(:,:)  ! the time-dependent cluster indices
integer,allocatable::conf_final(:)  ! the final list of configurations
integer,allocatable::trans_conf(:)  ! translation from old to new configuraiton number
integer,allocatable::confbas_final(:,:)  ! the final list of confs. for each basis func
integer,allocatable::bas_list_num(:)  ! the list of basis function numbers in ML_AB (eval)
integer,allocatable::bas_list_all(:,:,:)  ! the list of basis functions in ML_AB (eval)
integer,allocatable::bas_list_histo(:,:)  ! the histogram of basis functions
integer::mat_coord(2)  ! the current matrix indices
integer::index_act   ! counting for basis set histogram
logical::eval_mode ! Switch for just evaluate the content of a single ML_AB
logical::eval_stat(10)  ! the progress for evaluation loops
logical::no_adf  ! if no angular distribution functions shall be calculated
logical,allocatable::atom_used(:)  ! boolean mask for blocking of treated atoms
!   for time measurement
character(8)  :: date
character(10) :: time
character(5)  :: zone
integer,dimension(8) :: values
!   for Lapack call (matrix diagonalization)
character(len=1)::JOBZ,UPLO
integer::Nn,LDA,LDU,LDV,LWORK,INFO
real(kind=8),dimension(:),allocatable::WORK,W
real(kind=8),dimension(:,:),allocatable::A_mat


real(kind=8),allocatable::gradnorm_all(:)  ! gradient forms for all atoms

character(len=130)::a130

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
write(*,*) "PROGRAM mlff_select: selection of an effective set of"
write(*,*) " local reference configurations (basis sets) for a VASP"
write(*,*) " machine-learning force field from a given ML-AB file."
write(*,*) "The procedure is similar to the VASP command ML_MODE = select,"
write(*,*) " but is faster (since purely based on geometrical comparisons)"
write(*,*) " and has more options to be adjusted (global approach)"
write(*,*) "Usage: mlff_select [list of command line arguments]"
write(*,*) " mlff_select -help  prints a list of all commands."
write(*,*) "Give the -overview command to show an overview of utils4VASP."
write(*,*) 

if (command_argument_count() .eq. 0) then
   stop
end if        
!
!     If the list of all commands shall be printed with the help command
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:5))  .eq. "-help") then
      write(*,*) " ------ LIST OF COMMANDS ------ "
      write(*,*) " all settings without default value must be given explicitly by "
      write(*,*) " one of the commands, all others are optional." 
      write(*,*) " -ml_ab=[file1],[file2],... : The file names of the ML_AB files "
      write(*,*) "    that shall be processed by the program. Between 1 and 20 "
      write(*,*) "    can be given in total."
      write(*,*) " -eval : Evaluate the content of a single ML_AB file to gain"
      write(*,*) "    information about its learning process."
      write(*,*) " -nbasis=[number] : the [number] is the desired number of local"
      write(*,*) "    reference configurations (basis functions) per element in the"
      write(*,*) "    newly written ML_AB file. If this is chosen, each element"
      write(*,*) "    will have the same number of local ref. confs."
      write(*,*) "    Very small values (below 1000) might lead to problems with the"
      write(*,*) "    neighborhood allocations (and will give a bad ML_FF anyway)."
      write(*,*) " -nbasis_el=[el1:number1],[el2:number2],... : (Optional) give the"
      write(*,*) "    number of basis functions per element explicitly. If this is "
      write(*,*) "    done, one number must be given for each element appearing "
      write(*,*) "    in any of the ML_AB files! example: nbasis_el=Ga:3000,Pt:1000"
      write(*,*) " -cutoff=[value] : The radial and angular cutoffs used during the"
      write(*,*) "    learning (ML_RCUT1, ML_RCUT2). Only one global value can be "
      write(*,*) "    given for all ML_AB files. It might be reasonable to try "
      write(*,*) "    a value different to the ones during learning. DEFAULT: 5.0"
      write(*,*) " -grad_frac=[value] : percentage of basis functions to be "
      write(*,*) "    allocated for the largest gradient components in the "
      write(*,*) "    references. The [value]*nbasis atoms with the largest "
      write(*,*) "    will be taken all. DEFAULT: 0.1"
      write(*,*) " -train_div=[value] : the training set diversity (maximum: 1.0,"
      write(*,*) "    minimum: 0.0). The larger the value, the higher percentage"
      write(*,*) "    will be chosen based on different number of neighbor atoms,"
      write(*,*) "    such that the tails of the neighbor distribution will be "
      write(*,*) "    weighted larger (larger diversity). DEFAULT: 0.1"
      write(*,*) " -rdf_only : Deactivates the calculation of ADF profiles, only RDF"
      write(*,*) "    is used for environment comparisons. Recommended for large "
      write(*,*) "    cutoff distances, where calculations get a huge speedup."
      write(*,*) " -rdf2adf=[value] : Relative weighting of radial distribution "
      write(*,*) "    function (RDF) and angular distribution function (ADF) overlap"
      write(*,*) "    in constructing the clustering matrices. Example: -rdf2adf=2.0:"
      write(*,*) "    RDF will have 66% weight, ADF 33% (2 to 1). DEFAULT: 1.0"
      write(*,*) " -neigh_basmax=[number] : Maximum number of basis functions allocated"
      write(*,*) "    to a neighborhood diversity class. A larger number will increase"
      write(*,*) "    the amount of k-means clustering and slow down the calculation. "
      write(*,*) "    A smaller number might harm the diversity of the chosen basis "
      write(*,*) "    functions. DEFAULT: 10"
      write(*,*) " -bas_scale=('linear' or 'root') : How the number of basis functions"
      write(*,*) "    allocated to a certain neighborhood class shall scale with"
      write(*,*) "    the number of atoms contained into it."
      write(*,*) "    possibe are: 'linear' (linear scaling) or 'root' (square root"
      write(*,*) "    scaling). The linear scaling will give more spots to classes "
      write(*,*) "    with many atoms, whereas the root scaling will favor those"
      write(*,*) "    with lower number of atoms. DEFAULT: linear"
      write(*,*) " -max_environ=[number] : Maximum number of atoms within the cutoff "
      write(*,*) "    radius of any atom in the given structures. This value might "
      write(*,*) "    be raised if a larger cutoff shall be used. DEFAULT: 100"
      write(*,*) " -s_grid=[number] : Number of grid points for calculation of RDF"
      write(*,*) "    and ADF overlap integrals. DEFAULT: 200"
      write(*,*) " -rdf_exp=[value] : Exponential prefactor for line broadening "
      write(*,*) "    Gaussians in RDF profile calculations. DEFAULT: 20.0"
      write(*,*) " -adf_exp=[value] : Exponential prefactor for line broadening "
      write(*,*) "    Gaussians in ADF profile calculations. DEFAULT: 0.5"
      write(*,*) 
      stop
   end if
end do

call cpu_time(time1)
!
!     Open logfile for alternative printout to file 
!
open(unit=34,file="mlff_select.log",status="replace")
write(34,*) "Program mlff_select, execution started at:"
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
write(34,'(a,i2,a,i2,a,i4,a,i2,a,i2)') "    ",values(2),".",values(3),".",values(1), &
             & ", at ",values(5),":",values(6)
write(34,*)
flush(34)
!
!     Read in names of ML_AB files that shall be processed
!
mlab_num=0
mlab_list="xxx"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-ml_ab=") then
      read(arg(8:),*,iostat=readstat) mlab_list
      if (readstat .ne. 0) then
    !     write(*,*) "The format of the -ml_ab=[file1],[file2],... command seems to be corrupted!"
    !     write(*,*)
    !     stop
      end if
   end if
end do
mlab_num=0
do i=1,20
   if (mlab_list(i) .eq. "xxx") exit
   mlab_num=mlab_num+1
end do

if (mlab_num .lt. 1) then
   write(*,*) "Please give at least one ML_AB filename that can be read in!"
   write(*,*) "The format is -ml_ab=[file1],[file2], ... (up to 20 files possible)"
   write(*,*)
   write(34,*) "Please give at least one ML_AB filename that can be read in!"
   write(34,*) "The format is -ml_ab=[file1],[file2], ... (up to 20 files possible)"
   write(34,*)
   stop
end if

!
!     Read in the content of a single ML_AB file and analyze the learning process, if desired
!
eval_mode=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-eval") then
      eval_mode = .true.
   end if
end do

!
!     Read in command line arguments
!
!     The desired number of reference environments (basis functions)
!     Option 1: global number for all elements in the ML_AB files
!
if (.not. eval_mode) then
   nbasis=0
   basis_mode=0
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:8))  .eq. "-nbasis=") then
         read(arg(9:),*,iostat=readstat) inc
         if (readstat .ne. 0) then
            write(*,*) "The format of the -nbasis=[number] command seems to be corrupted!"
            write(*,*)
            write(34,*) "The format of the -nbasis=[number] command seems to be corrupted!"
            write(34,*)
            stop
         end if
         nbasis=inc
         basis_mode=1
      end if
   end do
else
   nbasis=1000
end if        
!
!     The desired number of reference environments (basis functions)
!     Option 2: one individual number for each element in the ML_AB file
!
basis_read="xxxx"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-nbasis_el=") then
      read(arg(12:),*,iostat=readstat) basis_read
      if (readstat .ne. 0) then
     !    write(*,*) "The format of the -nbasis=[number] command seems to be corrupted!"
     !    write(*,*)
     !    stop
      end if
      basis_mode=2
   end if
end do
el_list_bas="XX"
do i=1,50
   if (basis_read(i) .eq. "xxxx") exit
   do j=1,50
      if (basis_read(i)(j:j) .eq. ":") then
         el_list_bas(i)=basis_read(i)(1:j-1)
         read(basis_read(i)(j+1:),*) nbasis(i)   
      end if
   end do
end do

if (sum(nbasis) .lt. 1) then
   write(*,*) "Please give the number of local reference configurations after selection"
   write(*,*) " with the keyword -nbasis=[number] or -nbasis_el=[el1:number1],[el2:number2],..."
   write(*,*) "Recommended are values between 2000 and 8000"
   write(*,*)
   write(34,*) "Please give the number of local reference configurations after selection"
   write(34,*) " with the keyword -nbasis=[number] or -nbasis_el=[el1:number1],[el2:number2],..."
   write(34,*) "Recommended are values between 2000 and 8000"
   write(34,*)
   stop
end if

!
!     The cutoff used during the learning process
!
!
!     Default value for the cutoff: 5 Ang
!
cutoff=5.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-cutoff=") then
      read(arg(9:),*,iostat=readstat) cutoff
      if (readstat .ne. 0) then
         write(*,*) "The format of the -cutoff=[value] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -cutoff=[value] command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (cutoff .lt. 0.0d0) then
   write(*,*) "Please give a reasonable value for the radial and angular cutoffs "
   write(*,*) " used during the learning (or shall be used during selection)"
   write(34,*) "Please give a reasonable value for the radial and angular cutoffs "
   write(34,*) " used during the learning (or shall be used during selection)"
   stop
end if

!
!     The fraction of basis functions based on the largest gradient norms
!
grad_frac=0.1d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-grad_frac=") then
      read(arg(12:),*,iostat=readstat) grad_frac
      if (readstat .ne. 0) then
         write(*,*) "The format of the -grad_frac=[value] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -grad_frac=[value] command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (grad_frac .lt. 0.0 .or. grad_frac .gt. 1.0) then
   write(*,*) "Please give the desired fraction of basis functions for the largest "
   write(*,*) " gradient norm atoms with the keyword -grad_frac=[value]"
   write(*,*) " where the maximum fraction is 1.0 (all) and the minimum is 0.0 (none)"
   write(*,*)
   write(34,*) "Please give the desired fraction of basis functions for the largest "
   write(34,*) " gradient norm atoms with the keyword -grad_frac=[value]"
   write(34,*) " where the maximum fraction is 1.0 (all) and the minimum is 0.0 (none)"
   write(34,*)
   stop
end if

!
!     The training set diversity (1: max, 0 min)
!
train_div=0.1d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-train_div=") then
      read(arg(12:),*,iostat=readstat) train_div
      if (readstat .ne. 0) then
         write(*,*) "The format of the -train_div=[value] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -train_div=[value] command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (train_div .lt. 0.0 .or. train_div .gt. 1.0) then
   write(*,*) "Please give the desired trainset diversity with the keyword -train_div=[value]"
   write(*,*) " where the maximum diversity is 1.0 and the minimum is 0.0"
   write(*,*)
   write(34,*) "Please give the desired trainset diversity with the keyword -train_div=[value]"
   write(34,*) " where the maximum diversity is 1.0 and the minimum is 0.0"
   write(34,*)
   stop
end if
!
!     Deactivates the calculation of ADFs in order to accelerate the calculation
!
no_adf=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-rdf_only") then
      no_adf = .true.      
   end if
end do

!
!     Relative weight of RDF and ADF overlaps for the construction of 
!      the clustering matrices
!
rdf_adf=1.0d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-rdf2adf=") then
      read(arg(10:),*,iostat=readstat) rdf_adf
      if (readstat .ne. 0) then
         write(*,*) "The format of the -rdf2adf=[value] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -rdf2adf=[value] command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (rdf_adf .lt. 0.0 ) then
   write(*,*) "Please give a value between zero and infinity for the -rdf2adf command!"
   write(*,*)
   write(34,*) "Please give a value between zero and infinity for the -rdf2adf command!"
   write(34,*)
   stop
end if
rdf_weight=rdf_adf/(rdf_adf+1.d0)
adf_weight=1.d0/(rdf_adf+1.d0)

!
!     The number of neighborhood classes
!
neigh_maxsize=10
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:15))  .eq. "-neigh_basmax=") then
      read(arg(16:),*,iostat=readstat) neigh_maxsize
      if (readstat .ne. 0) then
         write(*,*) "The format of the -neigh_baxmax=[number] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -neigh_basmax=[number] command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (neigh_maxsize .lt. 1) then
   write(*,*) "Please give the maximum number of basis functions allocated to a "
   write(*,*) " neighborhood diversity class, must be larger than 2, but should not "
   write(*,*) " be too large (smaller than 20-30)!"
   write(*,*)
   write(34,*) "Please give the maximum number of basis functions allocated to a "
   write(34,*) " neighborhood diversity class, must be larger than 2, but should not "
   write(34,*) " be too large (smaller than 20-30)!"
   write(34,*)
   stop
end if

!
!     If the number of basis functions for the neighborhood classes shall be scaling
!      linear or square root with the number of contained atoms
!
bas_scale="linear"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-bas_scale=") then
      read(arg(12:),*,iostat=readstat) bas_scale
      if (readstat .ne. 0) then
         write(*,*) "The format of the -bas_scale=(linear or root) command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -bas_scale=(linear or root) command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if ((bas_scale .ne. "linear") .and. (bas_scale .ne. "root")) then
   write(*,*) "Please decide if the number of basis functions per neighborhood class shall "
   write(*,*) " scale linearly or square root wise with the number of contained atoms."
   write(*,*) " Linear scaling: -bas_scale=linear"
   write(*,*) " Square root scaling: -bas_scale=root"
   write(*,*)
   write(34,*) "Please decide if the number of basis functions per neighborhood class shall "
   write(34,*) " scale linearly or square root wise with the number of contained atoms."
   write(34,*) " Linear scaling: -bas_scale=linear"
   write(34,*) " Square root scaling: -bas_scale=root"
   write(34,*)
   stop
end if

!
!     Maximum number of atoms in environment (assumed)
!
max_environ=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-max_environ=") then
      read(arg(14:),*,iostat=readstat) max_environ
      if (readstat .ne. 0) then
         write(*,*) "The format of the -max_environ=(number) command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -max_environ=(number) command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (max_environ .le. 1) then
   write(*,*) "Please give a number for the maximum number of atoms within the cutoff "
   write(*,*) "  radius of any other atom! (default: 500)"
   write(*,*)
   write(34,*) "Please give a number for the maximum number of atoms within the cutoff "
   write(34,*) "  radius of any other atom! (default: 500)"
   write(34,*)
   stop
end if


!
!     Number of grid points of numerical radial distribution functions
!
ngrid=200
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-s_grid=") then
      read(arg(9:),*,iostat=readstat) ngrid
      if (readstat .ne. 0) then
         write(*,*) "The format of the -s_grid=(number) command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -s_grid=(number) command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (ngrid .le. 1) then
   write(*,*) "Please give a reasonable value for the number of RDF/ADF integration "
   write(*,*) " grid points! (default: 200)."   
   write(*,*)
   write(34,*) "Please give a reasonable value for the number of RDF/ADF integration "
   write(34,*) " grid points! (default: 200)."
   write(34,*)
   stop
end if

!
!     Exponent for RDF Gaussian
!
alpha_r=20.0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-rdf_exp=") then
      read(arg(9:),*,iostat=readstat) alpha_r
      if (readstat .ne. 0) then
         write(*,*) "The format of the -rdf_exp=[number] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -rdf_exp=[number] command seems to be corrupted!"
         write(34,*)

         stop
      end if
   end if
end do
if (alpha_r .lt. 0.0001d0) then
   write(*,*) "Please give the preexponential factor of the Gaussian used to broaden"
   write(*,*) " the RDF profiles of all atoms."
   write(*,*)
   write(34,*) "Please give the preexponential factor of the Gaussian used to broaden"
   write(34,*) " the RDF profiles of all atoms."
   write(34,*)
   stop
end if
!
!     Exponent for ADF Gaussian
!
alpha_a=0.5d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-adf_exp=") then
      read(arg(9:),*,iostat=readstat) alpha_a
      if (readstat .ne. 0) then
         write(*,*) "The format of the -adf_exp=[number] command seems to be corrupted!"
         write(*,*)
         write(34,*) "The format of the -adf_exp=[number] command seems to be corrupted!"
         write(34,*)
         stop
      end if
   end if
end do
if (alpha_a .lt. 0.0001d0) then
   write(*,*) "Please give the preexponential factor of the Gaussian used to broaden"
   write(*,*) " the ADF profiles of all atoms."
   write(*,*)
   write(34,*) "Please give the preexponential factor of the Gaussian used to broaden"
   write(34,*) " the ADF profiles of all atoms."
   write(34,*)
   stop
end if
!
!     Abort evaluation mode if more than one ML_AB is defined
!
if (eval_mode) then
   if (mlab_num .gt. 1) then
      write(*,*) "Please give only one ML_AB file for the evaluation mode!"
      stop
   end if
end if        
!
!     Define the Pi
!
pi=4d0*atan(1.0d0)
!
!     Default values for local and global elements array
!
el_list_glob="XX"
el_list_loc="XX"
!
!     Default values for global core charge array
!
ind_list=0
nelems_all=0
!
!     Length between two x values on RDF grid (Ang.)
!
dx=cutoff/real(ngrid)
!
!     Length between two angle values on ADF grid (degrees)
!
da=180.d0/real(ngrid)

!
!     Print out all settings for information
!
write(*,*) "------------- CALCULATION SETTINGS ---------------"
write(*,*) "List of ML_AB files to be read in:"
write(34,*) "------------- CALCULATION SETTINGS ---------------"
write(34,*) "List of ML_AB files to be read in:"
do i=1,mlab_num
   if (i .lt. 10) then
      write(*,'(a,i1,a,a)',advance="no") "   (",i,") ",trim(mlab_list(i))
      write(34,'(a,i1,a,a)',advance="no") "   (",i,") ",trim(mlab_list(i))
   else 
      write(*,'(a,i2,a,a)',advance="no") " (",i,") ",trim(mlab_list(i))
      write(34,'(a,i2,a,a)',advance="no") " (",i,") ",trim(mlab_list(i))
   end if        
end do
write(*,*)
if (eval_mode) then        
   write(*,*) "The content of the ML_AB file will be evaluated."
   write(34,*) "The content of the ML_AB file will be evaluated."
else   
   if (basis_mode .eq. 1) then
      write(*,*) "Number of basis functions for all elements: ",nbasis(1)
      write(34,*) "Number of basis functions for all elements: ",nbasis(1) 
   else if (basis_mode .eq. 2) then
      write(*,*) "Number of basis functions for each element, separately:"
      write(34,*) "Number of basis functions for each element, separately:"
      do i=1,50
         if (el_list_bas(i) .eq. "XX") exit
         write(*,'(a,a,a,i8)') "  * ",el_list_bas(i),":  ",nbasis(i)
         write(34,'(a,a,a,i8)') "  * ",el_list_bas(i),":  ",nbasis(i)
      end do
   end if
   write(*,'(a,f10.4,a)') " The ML_FF descriptor cutoff is ",cutoff," Angstrom."
   write(*,'(a,f10.4,a)') " ",grad_frac*100d0," % of basis functions allocated for large gradients."
   write(*,'(a,f10.4,a)') " ",train_div*100d0," % of basis functions allocated uniformly."
   if (no_adf) then
      write(*,*) "No ADF will be calculated! (-rdf_only activated)"
      rdf_weight=1.0
      adf_weight=0.d0
   end if       
   write(*,'(a,f8.3,a,f8.3,a)') " Weighting of overlap matrices: ",rdf_weight*100d0, &
                   & " % RDF, ",adf_weight*100d0," % ADF."
   write(*,'(a,i5)') " Number of neighborhood classes per environment & 
               &size:",neigh_classes 
   write(34,'(a,f10.4,a)') " The ML_FF descriptor cutoff is ",cutoff," Angstrom."
   write(34,'(a,f10.4,a)') " ",grad_frac*100d0," % of basis functions allocated for large gradients."
   write(34,'(a,f10.4,a)') " ",train_div*100d0," % of basis functions allocated uniformly."
   if (no_adf) then
      write(34,*) "No ADF will be calculated! (-rdf_only activated)"
   end if
   write(34,'(a,f8.3,a,f8.3,a)') " Weighting of overlap matrices: ",rdf_weight*100d0, &
                      & " % RDF, ",adf_weight*100d0," % ADF."
   write(34,'(a,i5)') " Number of neighborhood classes per environment & 
               &size:",neigh_classes

   if (bas_scale .eq. "linear") then
      write(*,*) "Number of basis functions per neighborhood class scaled linearly."
      write(34,*) "Number of basis functions per neighborhood class scaled linearly."
   else if (bas_scale .eq. "root") then     
      write(*,*) "Number of basis functions per neighborhood class scaled with square root."
      write(34,*) "Number of basis functions per neighborhood class scaled with square root."
   end if
   write(*,'(a,i6)') " Maximum number of atoms within the cutoff of any atom: ",max_environ
   write(*,'(a,i8)') " Number of grid points for RDF/ADF overlap integrations: ",ngrid
   write(*,'(a,f10.4)') " Exponential prefactor for RDF line broadening: ",alpha_r
   if (.not. no_adf) then
      write(*,'(a,f10.4)') " Exponential prefactor for ADF line broadening: ",alpha_a
   end if
   write(*,*) "--------------------------------------------------"
   write(*,*)
   write(34,'(a,i6)') " Maximum number of atoms within the cutoff of any atom: ",max_environ
   write(34,'(a,i8)') " Number of grid points for RDF/ADF overlap integrations: ",ngrid
   write(34,'(a,f10.4)') " Exponential prefactor for RDF line broadening: ",alpha_r
   if (.not. no_adf) then
      write(34,'(a,f10.4)') " Exponential prefactor for ADF line broadening: ",alpha_a
   end if
end if
write(*,*) "--------------------------------------------------"
write(*,*)
write(34,*) "--------------------------------------------------"
write(34,*)
flush(34)

! ####################################################
!     Loop over all ML_AB files given in the list and read in their content
!
!     PART A: HEADERS (for global array allocation)
!
natoms_max=0
do l=1,mlab_num
!
!     Open the ML_AB file and check if it's there
!
   open(unit=56,file=mlab_list(l),status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file ",trim(mlab_list(l)), " could not been found!"
      stop
   end if
   write(*,*) "Read header of the file ",trim(mlab_list(l))," ..."
   write(34,*) "Read header of the file ",trim(mlab_list(l))," ..."
   flush(34)
!
!     Read in the ML_AB file
!
!     First, read in the header and determine number of atoms and 
!      configurations
   read(56,*,iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file ",trim(mlab_list(l))," seems to be empty!"
      write(34,*) "The file ",trim(mlab_list(l))," seems to be empty!"
      stop
   end if
   read(56,*)
   read(56,*)
   read(56,*)
   read(56,*) ind_act
   conf_num=conf_num+ind_act
   read(56,*)
   read(56,*)
   read(56,*)
   read(56,*) nelems_local(l)
   read(56,*)
   read(56,*)
   read(56,*)
!
!     All elementwise properties are printed with three values per line!
!
   do i=1,int(nelems_local(l)/3) 
      read(56,*,iostat=readstat) el_list_loc((i-1)*3+1:i*3,l)
   end do
   if (int(nelems_local(l)/3)*3 .lt. nelems_local(l)) then
      if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 3) then
         read(56,*,iostat=readstat) el_list_loc(int(nelems_local(l)/3)*3+1:int(nelems_local(l)/3)*3+3,l)
      else if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 2) then     
         read(56,*,iostat=readstat) el_list_loc(int(nelems_local(l)/3)*3+1:int(nelems_local(l)/3)*3+2,l)
      else if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 1) then
         read(56,*,iostat=readstat) el_list_loc(int(nelems_local(l)/3)*3+1:int(nelems_local(l)/3)*3+1,l)     
      end if
   end if 
  
!
!     Determine different element indices
!
   outer: do i=1,50
      if (el_list_loc(i,l) .eq. "XX") exit
      call elem(el_list_loc(i,l),ind_act)
      inner: do j=1,nelems_local(l)
         if (ind_act .eq. ind_list_loc(j,l)) cycle outer
      end do inner 
      ind_list_loc(i,l) = ind_act
  !    nelems_local=nelems_local+1
   end do outer

   read(56,*)
   read(56,*)
   read(56,*)
!
!     Increment the global maximum number of atoms in the system, if needed
!
   read(56,*) ind_act
   if (ind_act .gt. natoms_max) then
      natoms_max=ind_act
   end if
   read(56,*)
   read(56,*)
   read(56,*)
   read(56,*) nat_type_max
   read(56,*)
   read(56,*)
   read(56,*)
   do i=1,int(nelems_local(l)/3)
      read(56,*)
   end do
   if (int(nelems_local(l)/3)*3 .lt. nelems_local(l)) then
      if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 3) then
         read(56,*)
      else if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 2) then
         read(56,*) 
      else if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 1) then
         read(56,*) 
      end if
   end if

   read(56,*) 
   read(56,*)
   read(56,*)
   at_mass(:,l)=0.d0 
   do i=1,int(nelems_local(l)/3)
      read(56,*,iostat=readstat) at_mass((i-1)*3+1:i*3,l)
   end do
   if (int(nelems_local(l)/3)*3 .lt. nelems_local(l)) then
      if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 3) then
         read(56,*,iostat=readstat) at_mass(int(nelems_local(l)/3)*3+1:int(nelems_local(l)/3)*3+3,l)
      else if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 2) then
         read(56,*,iostat=readstat) at_mass(int(nelems_local(l)/3)*3+1:int(nelems_local(l)/3)*3+2,l)
      else if ((nelems_local(l)-int(nelems_local(l)/3)*3) .eq. 1) then
         read(56,*,iostat=readstat) at_mass(int(nelems_local(l)/3)*3+1:int(nelems_local(l)/3)*3+1,l)
      end if
   end if
!
!    If the evaluation mode is requested: Read in also the lists of basis functions!
!
   if (eval_mode) then
      allocate(bas_list_num(nelems_local(1)))   
      read(56,*)
      read(56,*)
      read(56,*)
      read(56,*) bas_list_num
      allocate(bas_list_all(2,maxval(bas_list_num),nelems_local(1)))
      do i=1,nelems_local(1)
         read(56,*)
         read(56,*)
         read(56,*)
         do j=1,bas_list_num(i)
            read(56,*) bas_list_all(:,j,i)
         end do   
      end do
   end if        
   close(56)
end do


!
!     Define global element names and masses array
!

el_list_glob="XX"
inc=0
do i=1,mlab_num
   do_nelems1: do j=1,nelems_local(i)
      do k=1,inc
         if (trim(el_list_glob(k)) .eq. trim(el_list_loc(j,i))) then
            cycle do_nelems1
         end if
      end do   
      inc=inc+1
      el_list_glob(inc)=el_list_loc(j,i)   
   end do do_nelems1  
end do
nelems_glob=inc


at_mass_glob=0.d0
inc=0
do i=1,mlab_num
   do_nelems2: do j=1,nelems_local(i)
      do k=1,inc-1
         if (at_mass_glob(k) .eq. at_mass(j,i)) then
            cycle do_nelems2
         end if
      end do
      inc=inc+1
      at_mass_glob(inc)=at_mass(j,i)
   end do do_nelems2
end do

!
!     Compare global element list with basis functions element list
!     and resort the latter (or throw an error if something is missing)
!
if (basis_mode .eq. 2) then
   do i=1,nelems_glob
      do j=1,nelems_glob
         if (el_list_glob(i) .eq. el_list_bas(j)) then
            nbasis_tmp(j)=nbasis(i)
            cycle
         end if
      end do
   end do
   nbasis=nbasis_tmp
end if

!
!     Number of basis functions allocated to minimum into different neighbor bins
!
neigh_min_bas=0
do i=1,nelems_glob
   neigh_min_bas(i)=int(nbasis(i)*train_div)
end do

!
!     Determine global different element indices
!
do i=1,nelems_glob
   call elem(el_list_glob(i),ind_act)
   ind_list(i) = ind_act
end do 


!
!    PART B: BODY (for filling of global arrays)
!
!    The global configuration increment index
!
act_conf=0
!
!     Allocate global arrays
!
!     The number of atoms in all reference structures
allocate(natoms(conf_num))
!     The number of elements in all reference structures
allocate(nelems_confs(conf_num))
!     The ML_AB file index of the configuration
allocate(mlab_conf(conf_num))
!     The list of element symbols in all configurations
allocate(el_list_confs(nelems_glob,conf_num))
allocate(el_list_tmp(nelems_glob))
!     The list of atom numbers in all configurations
allocate(el_nums_confs(nelems_glob,conf_num))
allocate(el_nums_tmp(nelems_glob))
!     The geometries (cartesian coordinates)
allocate(xyz(3,natoms_max,conf_num))
!     The geometries (direct coordinates
allocate(xyz_dir(3,natoms_max,conf_num))
!     The unit cell shapes
allocate(cells(3,3,conf_num))
!     The CTIFOR values
allocate(ctifors(conf_num))
!     The element indices (core charges), more effective then names
allocate(inds(natoms_max,conf_num))
!     The reference energies
allocate(energies(conf_num))
!     The reference gradients
allocate(grads(3,natoms_max,conf_num))
!     The reference stress tensors
allocate(stress(6,conf_num))


do l=1,mlab_num
!
!     Again open the ML_AB file 
!
   open(unit=56,file=mlab_list(l),status="old")

   write(*,*) "Read body of the file ",trim(mlab_list(l))," ..."
   write(34,*) "Read body of the file ",trim(mlab_list(l))," ..."
   flush(34)
   do 
      read(56,'(a)',iostat=readstat) a130
      if (readstat .ne. 0) then
         exit
      end if
!
!     If a new configuration has been found, read it in
!
      if (index(a130,"Configuration") .ne. 0) then
         act_conf=act_conf+1
         do i=1,7
            read(56,*)
         end do 
!
!     Read the list and numbers of elements to fulle the inds array
!
         read(56,*) nelems_confs(act_conf)
         do i=1,3
            read(56,*)
         end do
!
!     To which ML_AB file the current configuration belongs
!
         mlab_conf(act_conf)=l
!
!     The current number of atoms
!
         read(56,*) natoms(act_conf)
         read(56,*)
         read(56,*)
         read(56,*)
         do i=1,nelems_confs(act_conf)
            read(56,*) el_list_tmp(i),el_nums_tmp(i) 
         end do
!
!     Sort the element list according to the global element list!
!
!     Further generate mapping array for resorting of falsely sorted
!     structures
!
         if (allocated(sort_map)) deallocate(sort_map)
         if (allocated(sort_map2)) deallocate(sort_map2)
         allocate(sort_map(natoms(act_conf)))
         allocate(sort_map2(natoms(act_conf)))
         inc2=0
         inc3=0
         do i=1,nelems_glob
            do j=1,nelems_confs(act_conf)
               if (el_list_tmp(j) .eq. el_list_glob(i)) then
                  inc2=inc2+1
                  el_list_confs(inc2,act_conf)=el_list_tmp(j)
                  el_nums_confs(inc2,act_conf)=el_nums_tmp(j)
                  do k=1,el_nums_tmp(j)
                     inc3=inc3+1
                     sort_map(inc3)=sum(el_nums_tmp(1:j-1))+k
                  end do
                  exit    
               end if
            end do
         end do
         do i=1,natoms(act_conf)
            sort_map2(sort_map(i))=i
         end do
!
!     Fill element index array
!
         inc=0
         do i=1,nelems_confs(act_conf)
            do j=1,el_nums_confs(i,act_conf)
               inc=inc+1
               call elem(el_list_confs(i,act_conf),ind_act)
               inds(inc,act_conf)=ind_act
            end do
         end do   
         read(56,*)
         read(56,*)
         read(56,*)
         read(56,*) ctifors(act_conf)
         read(56,*)
         read(56,*)
         read(56,*)
         
!
!     Read in the current unit cell vectors
!
         do i=1,3
            read(56,*) cells(:,i,act_conf)
         end do
         read(56,*)
         read(56,*)
         read(56,*)
!
!     Read in the current geometry
!
         do i=1,natoms(act_conf)
            read(56,*) xyz(:,sort_map2(i),act_conf)
         end do
!
!     Invert the unit cell matrix and calculate the geometry 
!        in direct coordinates
!
         call matinv3(cells(:,:,act_conf),cell_inv)
         do i=1,natoms(act_conf)
            xyz_dir(:,i,act_conf)=xyz(1,i,act_conf)*cell_inv(1,:)+ &
                 & xyz(2,i,act_conf)*cell_inv(2,:)+xyz(3,i,act_conf)*cell_inv(3,:)
         end do 
         read(56,*)
         read(56,*)
         read(56,*)
!
!     Read in the current energy
!
         read(56,*) energies(act_conf)
         read(56,*)
         read(56,*)
         read(56,*)
!
!     Read in the current gradient
!
         do i=1,natoms(act_conf)
            read(56,*) grads(:,sort_map2(i),act_conf)
         end do
         read(56,*)
         read(56,*)
         read(56,*)
         read(56,*)
         read(56,*)
!
!     Read in the current stress tensors
!
         read(56,*) stress(1:3,act_conf)
         read(56,*)
         read(56,*)
         read(56,*)
         read(56,*) stress(4:6,act_conf)
      end if
   end do
   close(56)
end do


write(*,*) " ... done!"
write(34,*) " ... done!"
flush(34)

!
!     If only evaluation is requested, do it here and stop afterwards
!
if (eval_mode) then
!
!     Print out the energy of each structure
!
   open(unit=112,file="mlab_energies.dat",status="replace")
   write(112,*) "# This file contains the energies of all structures stored "
   write(112,*) "#  into the ML_AB file" ,mlab_list(1)
   write(112,*) "# Evaluated with mlff_select -eval"
   write(112,*) "# Structure No.     energy(eV)"
   do i=1,conf_num
      write(112,*) i,energies(i)
   end do
   close(112)
   write(*,*) "File 'mlab_energies.dat' with total energies written."
!
!     Print out the CTIFOR value of each structure
!
   open(unit=112,file="mlab_CTIFORs.dat",status="replace")
   write(112,*) "# This file contains the CTIFOR values of all structures stored "
   write(112,*) "#  into the ML_AB file" ,mlab_list(1)
   write(112,*) "# Evaluated with mlff_select -eval"
   write(112,*) "# Structure No.     energy(eV)"
   do i=1,conf_num
      write(112,*) i,ctifors(i)
   end do
   close(112)
   write(*,*) "File 'mlab_CTIFORs.dat' with CTIFOR values written."
!
!     Print out the number of collected reference environments during the learning
!     depending on the structure number 
!     To make this histogram better visible, the frames will be packed together
!     to result in a total number of 200 histogram points
!
   allocate(bas_list_histo(300,nelems_local(1)))
   bas_list_histo=0
   do i=1,nelems_local(1)
      do j=1,bas_list_num(i)
         index_act=int((real(bas_list_all(1,j,i))/real(conf_num))*300)
         bas_list_histo(index_act,i) = bas_list_histo(index_act,i) + 1  
      end do
   end do
   bas_list_histo(300,:)=0
   open(unit=113,file="mlab_basis_histos.dat",status="replace")
   write(113,*) "# This file contains the number of local reference "
   write(113,*) "#  configurations for each element along the on the "
   write(113,*) "#  fly learning, packed together into bins, spanning"
   write(113,*) "#  some structures, each."
   write(113,*) "# Each value contains ca. ",real(conf_num)/real(300), & 
                    & " training structures"
   write(113,'(a)',advance="no") " # Structure No. "
   do i=1,nelems_local(1)
      write(113,'(a,a)',advance="no") "        ",el_list_glob(i)
   end do
   write(113,*)
   do i=1,300
      write(113,*) int(real(i)/300.0*conf_num),bas_list_histo(i,:)
   end do
   close(113)
   write(*,*) "File 'mlab_basis_histos.dat' with basis set histograms written."
   write(*,*) "Evaluation finished!"
   write(*,*)
   stop
end if        
!
!     Setup classifier arrays for the atoms
!     All atoms (no matter to which configuation they belong)
!     are stored into one large array for each classifier!
!
!     From now on, the number of elements is always nelems_glob
!
nelems=nelems_glob
!
!     First, allocation of arrays:
!
!     Total number of atoms within any of the configurations
natoms_sum=sum(natoms)
!     The nelems-dimensional increment array
allocate(incs(nelems))
!     The configuration number of the atom
allocate(confnum_all(natoms_sum))
!     The ML_AB file number of the atom
allocate(mlab_atom(natoms_sum))
!     Boolean mask, if atom was already allocated to final basis
allocate(atom_used(natoms_sum))
!     The atom number in the respective configuration
allocate(nat_all(natoms_sum))
!     The element index (core charge) of each atom
allocate(ind_all(natoms_sum))     
!     The gradient norm for each atom
allocate(gradnorm_all(natoms_sum))
!     The number of atoms in the surrounding (sorted by core charges)
allocate(num_around(nelems,natoms_sum))
!     The coordinates of atoms in the surrounding (for xyz printout)
allocate(environ(3,max_environ,natoms_sum))
!     The core charges/elements of atoms in surrounding
allocate(ind_env(max_environ,natoms_sum))
!     The list of distances to atoms within surrounding (local)
allocate(dist_list(max_environ))
!     The radial distribution functions around all atoms
allocate(rdf_all(ngrid,natoms_sum))
!     The angular distribution functions around all atoms
if (.not. no_adf) then        
   allocate(adf_all(ngrid,natoms_sum))
end if
!     The local list of angles in the environment
allocate(angles(max_environ*max_environ))
!     The global histogram with the number of neighbors
allocate(neigh_global(nelems,max_environ))
!     Number of basis functions per environment
allocate(neigh_bas(nelems,max_environ))
!     Number of atoms with gradient components in certain range
allocate(grad_histo(100))
!     The sum of all neighbors around all atoms
allocate(neighnum_global(natoms_sum))


!     The FINAL choice of basis functions!
allocate(final_choice(maxval(nbasis),nelems))




num_around=0
rdf_all=0.d0
if (.not. no_adf) then
   adf_all=0.d0
end if
write(*,*) "Total number of atoms (possible basis functions):",natoms_sum
write(34,*) "Total number of atoms (possible basis functions):",natoms_sum

write(*,*)
write(*,*) "Calculate classifiers for all atoms in the ML_AB file..."
write(34,*)
write(34,*) "Calculate classifiers for all atoms in the ML_AB file..."
flush(34)

eval_stat=.false.
inc=0
do i=1,conf_num
   do j=1,10
      if (real(i)/real(conf_num) .gt. real(j)*0.1d0) then
         if (.not. eval_stat(j)) then
            write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
            write(34,'(a,i4,a)')  "  ... ",j*10,"% done "
            flush(34)
            eval_stat(j) = .true.
         end if
      end if
   end do

   do j=1,natoms(i)
      inc=inc+1
      confnum_all(inc)=i
      nat_all(inc)=j
      ind_all(inc)=inds(j,i)
      mlab_atom(inc)=mlab_conf(i)
!
!     Calculate the gradient norm
!
      gradnorm_all(inc)=sqrt(grads(1,j,i)*grads(1,j,i)+grads(2,j,i)*&
                     & grads(2,j,i)+grads(3,j,i)*grads(3,j,i))
!
!     Calculate the number of atoms in the direct surrounding
!
!     The atom itself is located in the origin for environment
!     If an environment is smaller than the largest one, set all missing
!     atoms to the origin as well, same element as central atom
!
      environ(:,:,inc)=0.d0
      ind_env(:,inc)=ind_all(inc)
!
!     All distances below cutoff
!
      inc2=1
      do k=1,natoms(i)
         if (k .ne. j) then
            dist_vec=xyz_dir(:,j,i)-xyz_dir(:,k,i)
!
!     Apply image vectors for periodicity
!
!     a-component 
            do while (abs(dist_vec(1)) .gt. 0.5d0) 
               dist_vec(1)=dist_vec(1)-sign(1.d0,dist_vec(1))
            end do
!     b-component                 
            do while (abs(dist_vec(2)) .gt. 0.5d0)
               dist_vec(2)=dist_vec(2)-sign(1.d0,dist_vec(2))
            end do   
!     c-component                 
            do while (abs(dist_vec(3)) .gt. 0.5d0)
               dist_vec(3)=dist_vec(3)-sign(1.d0,dist_vec(3))
            end do
!   
!     Finally, convert distance vector back to cartesian coordinates
!
            dist_vec(:)=dist_vec(1)*cells(1,:,act_conf)+dist_vec(2)*cells(2,:,act_conf)+ &
                       & dist_vec(3)*cells(3,:,act_conf)
!
!     Calculate the distance
!
            dist_act=sqrt(dot_product(dist_vec,dist_vec))            
!
!     Determine if atom is within the cutoff. If yes, add it to the list of 
!      surrounding atoms
!      Fill the respective atom (its coordinates and index) into the local 
!      environment cluster structure
!
            if (dist_act .lt. cutoff) then
               do l=1,nelems
                  if (inds(k,i) .eq. ind_list(l)) then                          
                     num_around(l,inc)=num_around(l,inc)+1
                     inc2=inc2+1
                     environ(:,inc2,inc)=-dist_vec(:)
                     ind_env(inc2,inc)=inds(k,i)
                     dist_list(inc2-1)=dist_act
                  end if
               end do
            end if            
         end if
      end do
!
!     Generate radial distribution function for current environment
!      Precompute it for each atom and normalize its integral!
!
      normfac=0.d0
      do k=1,ngrid
         x_act=k*dx
         do l=2,inc2-1
            rdf_all(k,inc)=rdf_all(k,inc)+exp(-alpha_r*(x_act-dist_list(l))**2)
         end do
         normfac=normfac+rdf_all(k,inc)*rdf_all(k,inc)*dx
      end do 
!
!     Normalize its integral to become 1
!
      rdf_all(:,inc)=rdf_all(:,inc)/sqrt(normfac)

!
!     Generate angular distribution function for current environment
!  
      if (.not. no_adf) then
         angles=0.d0
         inc3=1
         do k=2,inc2-1
            do l=k+1,inc2-1
               angles(inc3)=acos(dot_product(environ(:,k,inc),environ(:,l,inc))/&
                               &sqrt(dist_list(k-1)**2*dist_list(l-1)**2))*180d0/pi 
               inc3=inc3+1        
            end do
         end do  
         normfac=0.d0 
         do k=1,ngrid
            x_act=k*da        
            do l=1,inc3-1
               adf_all(k,inc)=adf_all(k,inc)+exp(-alpha_a*(x_act-angles(l))**2)
            end do
            normfac=normfac+adf_all(k,inc)*adf_all(k,inc)*da
         end do
!
!     Normalize its integral to become 1
!
         adf_all(:,inc)=adf_all(:,inc)/sqrt(normfac)
      end if
   end do
end do
write(*,*) " ... finished!"
write(34,*) " ... finished!"

write(*,*)
write(*,*) " Element   No. of atoms       Basis size       atom usage (%)"
write(34,*)
write(34,*) " Element   No. of atoms       Basis size       atom usage (%)"
flush(34)

allocate(natoms_el(nelems))
do i=1,nelems
   inc=0
   do j=1,natoms_sum
      if (ind_all(j) .eq. ind_list(i)) then
         inc=inc+1
      end if        
   end do
   natoms_el(i)=inc
   write(*,'(3a,i10,a,i10,a,f12.6)') "   ", el_list_glob(i)," : ",inc, "         ",nbasis(i), &
                  & "          ",real(real(nbasis(i))/real(inc)*100.d0)
   write(34,'(3a,i10,a,i10,a,f12.6)') "   ", el_list_glob(i)," : ",inc, "         ",nbasis(i), &
                  & "          ",real(real(nbasis(i))/real(inc)*100.d0)
end do
flush(34)
write(*,*)
write(34,*)
!
!     All atoms are still available
!
atom_used=.true.

!
!     Test write out for radial distribution functions
!
open(unit=57,file="rdf_test.dat",status="replace")
do i=1,ngrid
   write(57,*) i*dx,rdf_all(i,1:5)
end do
close(57)

!
!     Test write out angular distribution functions
!
if (.not. no_adf) then
   open(unit=58,file="adf_test.dat",status="replace")
   do i=1,ngrid
      write(58,*) i*da,adf_all(i,1:5)
   end do
   close(58)
end if

write(*,*) "Select reference confs. based on largest gradient norms..."
write(34,*) "Select reference confs. based on largest gradient norms..."
flush(34)
!
!    A: GRADIENT EXTREMA PRESELECTION:
!    Determine histogram of gradient norms for all atoms, allocate them into
!    bins according to their gradient norms
!    Choose the nbasis*grad_frac atoms with the largest gradient norms
!
!    The number of atoms preselected due to large gradient components
!
grad_pre=0
do i=1,nelems
   grad_pre(i)=int(nbasis(i)*grad_frac)
end do
allocate(atsel_grad(maxval(grad_pre),nelems))
allocate(atsel_grad_val(maxval(grad_pre),nelems))


grad_min=minval(gradnorm_all)
grad_max=maxval(gradnorm_all)
grad_step=(grad_max-grad_min)/100.d0
grad_histo=0

outer2: do i=1,natoms_sum
   inner2: do j=1,100
      if (gradnorm_all(i) .lt. j*grad_step) then
         grad_histo(j)=grad_histo(j)+1
         exit inner2
      end if
   end do inner2
end do outer2
!
!    Write distributions of gradient norms to file
!
open(unit=60,file="grad_histo.dat",status="replace")
write(60,*) "#Here, the number of atoms with gradient norms in a certain"
write(60,*) "#range (x-coordinate) are listed."
do i=1,100
   write(60,*) i*grad_step,grad_histo(i)

end do
close(60)

!
!    Sort the grad_pre atoms with the largest gradients into an array with
!    their atom indices, do a max-heap search
!  
!    Stored in array: atsel_grad
!
incs=0
atsel_grad=0
atsel_grad_val=0.d0
do
   inc=maxloc(gradnorm_all,dim=1)
   do i=1,nelems
      if (ind_all(inc) .eq. ind_list(i)) then
         if (incs(i) .lt. grad_pre(i)) then
            incs(i)=incs(i)+1
            atsel_grad(incs(i),i)=inc
            atom_used(inc)=.false. 
!    Store atom indices in final basis set!
            final_choice(incs(i),i)=inc
            atsel_grad_val(incs(i),i)=gradnorm_all(inc)
         end if
      end if
   end do
   if (sum(incs) .eq. sum(grad_pre)) then
      exit
   end if
   gradnorm_all(inc)=0.d0
end do

write(*,*) " ... done!"
write(*,*) "Sort the configs. according to the number of their neighbors..."
write(34,*) " ... done!"
write(34,*) "Sort the configs. according to the number of their neighbors..."
flush(34)
!    
!    B: THE LOCAL NEIGHBORS: number and diversity, minimum number!
!
neigh_global=0
neighnum_global=0
do i=1,natoms_sum
   if (atom_used(i)) then
      do j=1,nelems
         if (ind_all(i) .eq. ind_list(j)) then
            inc=sum(num_around(:,i))
!
!    If an atom has no neighbors, sort it into the one neighbor bin
!
            if (inc .eq. 0) then
               inc=inc+1
            end if        
            neigh_global(j,inc)=neigh_global(j,inc)+1
            neighnum_global(i)=neighnum_global(i)+inc
         end if
      end do     
   end if      
end do

write(*,*) " ... done!"
write(*,*) "Fill diversity-dependent minimum conf. numbers to neighbor bins..."
write(34,*) " ... done!"
write(34,*) "Fill diversity-dependent minimum conf. numbers to neighbor bins..."
flush(34)
!
!    Fill diversity-dependent minimum number of allocated basis functions
!      into different neighbor bins (each bin one, until full or all atoms
!      allocated, line per line)
!
neigh_bas=0
do_elems: do i=1,nelems
   inc=grad_pre(i)
   do_column:  do
      do_environ: do j=1,max_environ
         if (inc .ge. neigh_min_bas(i)+grad_pre(i)) exit do_column
         if (neigh_bas(i,j) .lt. neigh_global(i,j)) then
            neigh_bas(i,j)=neigh_bas(i,j)+1
            do_atoms: do k=1,natoms_sum
!
!      Allocate chosen atom to final basis set
!      Check if it has the correct element!
!
               if (atom_used(k) .and. (ind_all(k) .eq. ind_list(i))) then
                  if (neighnum_global(k) .eq. j) then
                     atom_used(k) = .false.
                     inc=inc+1
        !             neigh_bas(i,j)=neigh_bas(i,j)+1
                     final_choice(inc,i)=k
                     exit do_atoms
                  end if
               end if
            end do do_atoms
         end if
      end do do_environ
   end do do_column
end do do_elems
write(*,*) " ... done!"
write(34,*) " ... done!"
flush(34)
!
!    Write distribution of neighbor numbers to file
!
open(unit=59,file="neighbor_nums.dat",status="replace")
write(59,*) "#Here, the number of neighbors around atoms in ML_AB are"
write(59,*) "#summed and packed into bins, where the number of atoms with "
write(59,*) "#with an environment of certain size (total number) is listed,"
write(59,*) "#sorted by element."
write(59,'(a)',advance="no") "#  No. neighbors      "
do i=1,nelems
   write(59,'(a,a)',advance="no") el_list_glob(i),"          "
end do
write(59,*)

do i=1,max_environ
   write(59,*) i,neigh_global(1:nelems,i)!,neigh_bas(:,i)

end do
close(59)

write(*,*) "File 'neighbor_nums.dat' with neighbor number histograms written..."
write(*,*)
write(*,*) "Select remaining basis functions via hierarchial cluster analysis"
write(*,*) " of neighborhood-diversity specific atom classes, based on RDF "
write(*,*) " and ADF overlap matrices."
write(34,*) "File 'neighbor_nums.dat' with neighbor number histograms written..."
write(34,*)
write(34,*) "Select remaining basis functions via hierarchial cluster analysis"
write(34,*) " of neighborhood-diversity specific atom classes, based on RDF "
write(34,*) " and ADF overlap matrices."
flush(34)
!    
!     C: THE LOCAL NEIGHBORS: number and diversity, final allocation via 
!      hierarchial cluster analysis!
!      Outer loop over elements: separate calculation for each element!
!
do i=1,nelems
   write(*,*)
   write(*,*) "Cluster all ",el_list_glob(i)," atoms ..."
   write(34,*)
   write(34,*) "Cluster all ",el_list_glob(i)," atoms ..."
   flush(34)
!
!     Number of available basis functions:
!

   inc=0
   inc3=0
   do j=1,max_environ
      if ((neigh_global(i,j)-neigh_bas(i,j)) .gt. 0) then
         inc=inc+1
      end if
   end do
!
!    First determine number of neighborhood diversity classes such that the 
!    maximum number of basis functions to be clustered in one of the classes 
!    it not larger than the given threshold in order to keep the calculation effort
!    manageable
! 
!    BEGIN OF NEIGH_CLASSES CALCULATION
!
   if (allocated(basis_sizes)) deallocate(basis_sizes)
   if (allocated(basis_classes)) deallocate(basis_classes)
   if (allocated(neigh_div1)) deallocate(neigh_div1)
   if (allocated(neigh_div2)) deallocate(neigh_div2)
   if (allocated(bas_neighs)) deallocate(bas_neighs)
!
!    Similar to below, but now with neigh_classes=1 (total number per environment)
!
   allocate(basis_sizes(inc*1+1))
   allocate(bas_neighs(inc*1+1))
   allocate(basis_classes(maxval(neigh_global-neigh_bas),inc*1+1))
   allocate(neigh_div1(2,maxval(neigh_global-neigh_bas)))
   allocate(neigh_div2(2,maxval(neigh_global-neigh_bas)))
   remain_class=inc*1+1
   inc_remain=0
   basis_sizes=0
   basis_classes=0

!
!    Calculate basis function numbers per total environment class, without subdivision
!
   inc=0
   do j=1,max_environ
      neigh_div1=1
      neigh_div2=1
      inc=0
      do k=1,natoms_sum
         if (atom_used(k) .and. (ind_all(k) .eq. ind_list(i))) then
            if (neighnum_global(k) .eq. j) then
               inc=inc+1
               neigh_div1(1,inc)=k
               do l=1,nelems
                  if (num_around(l,k) .gt. 0) then
                     neigh_div1(2,inc)=neigh_div1(2,inc)*num_around(l,k)
                  end if
               end do
            end if     
         end if 
      end do
!
!    Second, sort the array with the surroundings
!
      do k=1,inc
         inc2=maxloc(neigh_div1(2,1:inc),dim=1)
         neigh_div2(:,k)=neigh_div1(:,inc2)
         neigh_div1(:,inc2)=1 
      end do
!
!    Third, allocate the atoms into the neighborhood class, and with this 
!    into the final cluster for the k-means clustering!
! 
      if (inc .gt. 0) then
         inc3=inc3+1
         inc4=0
         do k=1,1-1
            if ((inc/1) .lt. 10) then
               basis_sizes(remain_class)=basis_sizes(remain_class)+inc/1
               bas_neighs(remain_class)=max_environ+1
            else
               basis_sizes((inc3-1)*1+k)=inc/1  
               bas_neighs((inc3-1)*1+k)=j
            end if
            if ((inc/1) .lt. 10) then
               do l=1,inc/1
                  inc4=inc4+1
                  inc_remain=inc_remain+1
                  basis_classes(inc_remain,remain_class)=neigh_div2(1,inc4)
               end do
            else 
               do l=1,inc/1
                  inc4=inc4+1
                  basis_classes(l,(inc3-1)*1+k)=neigh_div2(1,inc4)
               end do 
            end if
         end do  
         if ((inc-inc4) .lt. 10) then 
            basis_sizes(remain_class)=basis_sizes(remain_class)+inc-inc4
            bas_neighs(remain_class)=max_environ+1
         else
            basis_sizes(inc3*1)=inc-inc4
            bas_neighs(inc3*1)=j
         end if
         if ((inc-inc4) .lt. 10) then
            do l=1,inc-inc4
               inc4=inc4+1
               inc_remain=inc_remain+1
               basis_classes(inc_remain,remain_class)=neigh_div2(1,inc4)
            end do            
         else
            do l=1,basis_sizes(inc3*1)
               inc4=inc4+1
               basis_classes(l,inc3*1)=neigh_div2(1,inc4)
            end do
         end if
      end if
   end do
!
!    Fourth, determine the number of available spots in the final basis set
!     for the different local environment diversity classes  
!
!     Total number of available spots
!
   spots_avail=nbasis(i)-grad_pre(i)-neigh_min_bas(i)
!
!     The first index to add an atom to the global basis function list
!
   inc_atom=nbasis(i)-spots_avail+1
!
!     Array with number of spots per neighbor subclass
!
   if (allocated(k_number)) deallocate(k_number)
   allocate(k_number(size(basis_sizes)))
!
   do j=1,size(basis_sizes)
      if (bas_scale .eq. "linear") then
         inc=inc+basis_sizes(j)
      else if (bas_scale .eq. "root") then
         inc=inc+int(sqrt(real(basis_sizes(j))))
      end if
   end do

!
!     Determine the numbers
!   
   do j=1,size(basis_sizes)
!
!     Linear scaling with number of atoms
!
      if (bas_scale .eq. "linear") then
         k_number(j)=nint(real(basis_sizes(j))/real(inc)*real(spots_avail))
!
!     Square root scaling with number of atoms
!
      else if (bas_scale .eq. "root") then
         k_number(j)=nint(sqrt(real(basis_sizes(j)))/real(inc)*real(spots_avail))

      end if
   end do
   neigh_classes=ceiling(real(maxval(k_number))/real(neigh_maxsize))
   write(*,'(a,i5,a)') " (Number of neighborhood diversity classes: ",neigh_classes,")"
! 
!    END OF NEIGH_CLASSES CALCULATION
!
!
!     Number of available basis functions:
!

   inc=0
   inc3=0
   do j=1,max_environ
      if ((neigh_global(i,j)-neigh_bas(i,j)) .gt. 0) then
         inc=inc+1
      end if
   end do

!    Allocate array with atom classes based on neighborhood size 
!     and diversity     
!
   if (allocated(basis_sizes)) deallocate(basis_sizes)
   if (allocated(basis_classes)) deallocate(basis_classes)
   if (allocated(neigh_div1)) deallocate(neigh_div1)
   if (allocated(neigh_div2)) deallocate(neigh_div2)
   if (allocated(bas_neighs)) deallocate(bas_neighs)
   if (allocated(natoms_neigh)) deallocate(natoms_neigh) 
   if (allocated(nbasis_neigh)) deallocate(nbasis_neigh)
!
!     Add an extra dimension (+1) as remainder for all classes with less 
!     then 10 atoms
!
   allocate(basis_sizes(inc*neigh_classes+1))
   allocate(bas_neighs(inc*neigh_classes+1))
   allocate(basis_classes(maxval(neigh_global-neigh_bas),inc*neigh_classes+1))
   allocate(neigh_div1(2,maxval(neigh_global-neigh_bas)))
   allocate(neigh_div2(2,maxval(neigh_global-neigh_bas)))
   allocate(natoms_neigh(max_environ+1))
   allocate(nbasis_neigh(max_environ+1))

   remain_class=inc*neigh_classes+1
   inc_remain=0
   basis_sizes=0
   basis_classes=0
   natoms_neigh=0
   nbasis_neigh=0
!
!    Subdivide the number of neighbor bins into regions based on neighborhood
!    diversity (how many of the different elements are within the neighborhood)
!
!    Fill local array with product of all numbers of elements in the surrounding
!
   inc=0
   do j=1,max_environ
      neigh_div1=1
      neigh_div2=1
      inc=0
      do k=1,natoms_sum
         if (atom_used(k) .and. (ind_all(k) .eq. ind_list(i))) then
            if (neighnum_global(k) .eq. j) then
               inc=inc+1
               neigh_div1(1,inc)=k
               do l=1,nelems
                  if (num_around(l,k) .gt. 0) then
                     neigh_div1(2,inc)=neigh_div1(2,inc)*num_around(l,k)
                  end if
               end do
            end if     
         end if 
      end do
      natoms_neigh(j)=inc
!
!    Second, sort the array with the surroundings
!
      do k=1,inc
         inc2=maxloc(neigh_div1(2,1:inc),dim=1)
         neigh_div2(:,k)=neigh_div1(:,inc2)
         neigh_div1(:,inc2)=1 
      end do
!
!    Third, allocate the atoms into the neighborhood class, and with this 
!    into the final cluster for the k-means clustering!
! 
      if (inc .gt. 0) then
         inc3=inc3+1
         inc4=0
         do k=1,neigh_classes-1
            if ((inc/neigh_classes) .lt. 10) then
               basis_sizes(remain_class)=basis_sizes(remain_class)+inc/neigh_classes
               bas_neighs(remain_class)=max_environ+1
            else
               basis_sizes((inc3-1)*neigh_classes+k)=inc/neigh_classes  
               bas_neighs((inc3-1)*neigh_classes+k)=j
            end if
            if ((inc/neigh_classes) .lt. 10) then
               do l=1,inc/neigh_classes
                  inc4=inc4+1
                  inc_remain=inc_remain+1
                  basis_classes(inc_remain,remain_class)=neigh_div2(1,inc4)
               end do
            else 
               do l=1,inc/neigh_classes
                  inc4=inc4+1
                  basis_classes(l,(inc3-1)*neigh_classes+k)=neigh_div2(1,inc4)
               end do 
            end if
         end do  
         if ((inc-inc4) .lt. 10) then 
            basis_sizes(remain_class)=basis_sizes(remain_class)+inc-inc4
            bas_neighs(remain_class)=max_environ+1
         else
            basis_sizes(inc3*neigh_classes)=inc-inc4
            bas_neighs(inc3*neigh_classes)=j
         end if
         if ((inc-inc4) .lt. 10) then
            do l=1,inc-inc4
               inc4=inc4+1
               inc_remain=inc_remain+1
               basis_classes(inc_remain,remain_class)=neigh_div2(1,inc4)
            end do            
         else
            do l=1,basis_sizes(inc3*neigh_classes)
               inc4=inc4+1
               basis_classes(l,inc3*neigh_classes)=neigh_div2(1,inc4)
            end do
         end if
      end if
   end do
!
!    Fourth, determine the number of available spots in the final basis set
!     for the different local environment diversity classes  
!
!     Total number of available spots
!
   spots_avail=nbasis(i)-grad_pre(i)-neigh_min_bas(i)
!
!     The first index to add an atom to the global basis function list
!
   inc_atom=nbasis(i)-spots_avail+1
!
!     Array with number of spots per neighbor subclass
!
   if (allocated(k_number)) deallocate(k_number)
   allocate(k_number(size(basis_sizes)))
!
   do j=1,size(basis_sizes)
      if (bas_scale .eq. "linear") then
         inc=inc+basis_sizes(j)
      else if (bas_scale .eq. "root") then
         inc=inc+int(sqrt(real(basis_sizes(j))))
      end if
   end do

!
!     Determine the numbers
!   
!   write(*,*) basis_sizes
   do j=1,size(basis_sizes)
!
!     Linear scaling with number of atoms
!
      if (bas_scale .eq. "linear") then
         k_number(j)=nint(real(basis_sizes(j))/real(inc)*real(spots_avail))
!
!     Square root scaling with number of atoms
!
      else if (bas_scale .eq. "root") then
         k_number(j)=nint(sqrt(real(basis_sizes(j)))/real(inc)*real(spots_avail))

      end if
   end do

!
!     If the total number of basis functions is too high or too low, increase 
!     or decrease them until the exact number is reached
!
   if (sum(k_number(1:size(basis_sizes))) .gt. spots_avail) then
      inc_remain=sum(k_number(1:size(basis_sizes)))-spots_avail
      decrem: do 
         do j=1,size(basis_sizes)
            if (k_number(j) .ge. 5) then
               k_number(j) = k_number(j)-1
               inc_remain=inc_remain-1
               if (inc_remain .lt. 1) exit decrem
            end if
         end do
      end do decrem
   end if
   if (sum(k_number(1:size(basis_sizes))) .lt. spots_avail) then
      inc_remain=spots_avail-sum(k_number(1:size(basis_sizes)))
      increm: do
         do j=1,size(basis_sizes)
            if (k_number(j) .ge. 5) then
               k_number(j) = k_number(j)+1
               inc_remain=inc_remain-1
               if (inc_remain .lt. 1) exit increm
            end if
         end do
      end do increm
   end if
!
!     Now perform the k-means clustering for each diversity-neighbor 
!     number group and obtain final choice of atoms for basis set
!
   do j=1,size(basis_sizes)
      do k=1,max_environ+1
         if (bas_neighs(j) .eq. k) then
            natoms_neigh(k)=natoms_neigh(k)+basis_sizes(j)
            nbasis_neigh(k)=nbasis_neigh(k)+k_number(j)
         end if        
      end do  
   end do

   do j=1,size(basis_sizes)
      if (allocated(mat_overlap)) deallocate(mat_overlap)
      allocate(mat_overlap(basis_sizes(j),basis_sizes(j)))
      if (allocated(hierarchy)) deallocate(hierarchy)
      allocate(hierarchy(basis_sizes(j),basis_sizes(j)))
      mat_overlap=1.1d0
      if (k_number(j) .gt. 0) then
         if (j .gt. 1) then
            if (bas_neighs(j) .ne. bas_neighs(j-1)) then
               if (j .eq. size(basis_sizes)) then
                  write(*,'(a,i8,a,i5,a)') " Cluster the remaining atoms (", &
                            & natoms_neigh(bas_neighs(j))," atoms to ",nbasis_neigh(bas_neighs(j)), &
                            & " basis functions) ..." 
                  write(34,*) " Cluster the remaining atoms ..."
               else
                  write(*,'(a,i5,a,i8,a,i5,a)') "  Cluster atoms with ",bas_neighs(j)," neighbors (", &
                            & natoms_neigh(bas_neighs(j))," atoms to ",nbasis_neigh(bas_neighs(j)), &
                            & " basis functions) ..."
                  write(34,'(a,i5,a)') "  Cluster atoms with ",bas_neighs(j)," neighbors ..."
               end if         
               flush(34)
            end if
         end if     
                    
         do k=1,basis_sizes(j)
            do l=1,k-1
!
!     Calculate radial distribution function (RDF) and angular distribution
!      function (ADF) overlaps of all atoms in the cluster
!      Calculate 1-overlap, in order to get the lowest values for the 
!      largest overlaps (diagonal: zero)
!
               rdf_overlap=0.d0
               do m=1,ngrid
                  rdf_overlap=rdf_overlap+rdf_all(m,basis_classes(k,j))* &
                              & rdf_all(m,basis_classes(l,j))*dx
               end do
               adf_overlap=0.d0
               if (.not. no_adf) then
                  do m=1,ngrid
                     adf_overlap=adf_overlap+adf_all(m,basis_classes(k,j))* &
                              & adf_all(m,basis_classes(l,j))*da                
                  end do
               end if
!
!     Weight the RDF and ADF overlaps according to the value given within
!      the -rdf2adf command
!
               mat_overlap(k,l)=1.d0-(rdf_overlap*rdf_weight+adf_overlap*adf_weight)
            end do
         end do 
!         if (basis_sizes(j) .gt. 98) then
!         do k=1,basis_sizes(j)
!            do l=1,k-1
!               write(97,*) k,l,mat_overlap(k,l)
!            end do
!            do l=k,basis_sizes(j)
!               write(97,*) k,l,mat_overlap(l,k)
!            end do
!            write(97,*)
!         end do
!         stop "Gupog"
!         end if        

!
!     Perform the hierarchial clustering of the current RDF overlap matrix
!
!     The initial hierarchy: all indices as usual
!
         do l=1,basis_sizes(j)
            hierarchy(l,1)=l
         end do
         inc2=basis_sizes(j)+1
         do k=2,basis_sizes(j)
            do l=1,basis_sizes(j)
               hierarchy(l,k)=hierarchy(l,k-1)
            end do
!
!     Determine indices with lowest current value and update cluster list
!
            mat_coord=minloc(mat_overlap) 
            hierarchy(mat_coord(1),k)=inc2
            hierarchy(mat_coord(2),k)=inc2
            rdum=mat_overlap(mat_coord(1),mat_coord(2))
!
!     If the current index is already a cluster, set all of its atoms to the 
!      new cluster
!
            do l=1,basis_sizes(j)
               if (hierarchy(l,k) .eq. hierarchy(mat_coord(1),k-1)) then
                  hierarchy(l,k) = inc2
               end if
               if (hierarchy(l,k) .eq. hierarchy(mat_coord(2),k-1)) then
                  hierarchy(l,k) = inc2
               end if
            end do
      
!
!     Update all elements in the row and columns of the minloc values to 
!     the larger values of the respective column/row
!
            do l=1,basis_sizes(j)
               if ((mat_overlap(l,mat_coord(1)) .lt. 1.d0) .and. & 
                    & (mat_overlap(l,mat_coord(2)) .lt. 1.d0)) then
                  if (mat_overlap(l,mat_coord(1)) .gt. mat_overlap(l,mat_coord(2))) then
                     mat_overlap(l,mat_coord(2)) = mat_overlap(l,mat_coord(1))
                  else
                     mat_overlap(l,mat_coord(1)) = mat_overlap(l,mat_coord(2))
                  end if 
               end if
            end do

            do l=1,basis_sizes(j)
               if ((mat_overlap(mat_coord(1),l) .lt. 1.d0) .and. &
                    & (mat_overlap(mat_coord(2),l) .lt. 1.d0)) then
                  if (mat_overlap(mat_coord(1),l) .gt. mat_overlap(mat_coord(2),l)) then
                     mat_overlap(mat_coord(2),l) = mat_overlap(mat_coord(1),l)
                  else
                     mat_overlap(mat_coord(1),l) = mat_overlap(mat_coord(2),l)
                  end if
               end if
            end do
!
!     Set all values that appear twice or more in the merged cluster to 1.0
!
            inc3=1
            do l=1,basis_sizes(j)
               if (hierarchy(l,k) .eq. inc2) then
                  if (inc3 .gt. 1) then
                     mat_overlap(:,l) = 1.d0
                  end if
                  inc3=inc3+1
               end if
            end do
            inc3=1
            do l=1,basis_sizes(j)
               if (hierarchy(l,k) .eq. inc2) then
                  if (inc3 .gt. 1) then
                     mat_overlap(l,:) = 1.d0
                  end if
                  inc3=inc3+1
               end if
            end do

         !   do l=1,basis_sizes(j)
         !   write(199,'(12f10.6)') mat_overlap(l,:)
         !   end do
         !   write(199,*)

            inc2=inc2+1       
         end do
!
!     Finally, determine the atom indices used for the basis set by choosing them 
!      from the intermediate clustering, where the number of clusters is equal to
!      the number of needed basis functions, one function is then taken from
!      each cluster
!     The first atom always belongs to a new cluster and is thus always taken
!
         if (inc_atom .gt. nbasis(i)) cycle 
         final_choice(inc_atom,i)=basis_classes(1,j)
         inc_atom=inc_atom+1 

         inc_act=basis_sizes(j)-k_number(j)+1
         if (k_number(j) .gt. basis_sizes(j)) then
            write(*,*) "For one neighborhood class, the number of basis functions to "
            write(*,*) " be allocated is larger than the total number of atoms in this"
            write(*,'(a,i7,a,i7,a)') " spot! (current:",k_number(j),", allowed:",basis_sizes(j)," )"
            write(*,*) "Please reduce the number of neigh_classes (globally) or the "
            write(*,*) "  number of basis functions for ",el_list_glob(i),"!"
            write(34,*) "For one neighborhood class, the number of basis functions to "
            write(34,*) " be allocated is larger than the total number of atoms in this"
            write(34,'(a,i7,a,i7,a)') " spot! (current:",k_number(j),", allowed:",basis_sizes(j)," )"
            write(34,*) "Please reduce the number of neigh_classes (globally) or the "
            write(34,*) "  number of basis functions for ",el_list_glob(i),"!"
            stop 
         end if        
         atom_list: do k=2,basis_sizes(j)
            
!
!     The cluster must be different to that before but also should not appeared earlier at all!
!
            if (hierarchy(k-1,inc_act) .ne. hierarchy(k,inc_act)) then
               do l=1,k-1
                  if (hierarchy(l,inc_act) .eq. hierarchy(k,inc_act)) then 
                     cycle atom_list
                  end if                            
               end do   
               if (inc_atom .gt. nbasis(i)) cycle atom_list   
               final_choice(inc_atom,i)=basis_classes(k,j)
               inc_atom=inc_atom+1              
            end if        
         end do atom_list

      end if
   end do

end do
write(*,*) " ... done!"
write(*,*) 
write(*,*) "Prepare final writeout of selected configurations..."
write(34,*) " ... done!"
write(34,*)
write(34,*) "Prepare final writeout of selected configurations..."
flush(34)

!do i=1,nelems
!   do j=1,nbasis(i)
!      write(23,*) i,j,final_choice(j,i)
!   end do
!end do


!
!     Reorder final atom selection to ML_AB files
!
do i=1,nelems
   inc=0
   if (allocated(trans_final)) deallocate(trans_final)
   allocate(trans_final(nbasis(i)))
   do j=1,mlab_num
      do k=1,nbasis(i)
         if (mlab_atom(final_choice(k,i)) .eq. j) then
            inc=inc+1
            trans_final(inc)=final_choice(k,i)
         end if      
      end do
   end do
   final_choice(:,i)=trans_final(:)
end do
!
!     Finally, write out the newly generated ML_AB file!
!
!     Determine the number and the list of written out full configurations
!
allocate(conf_final(conf_num))
allocate(trans_conf(conf_num))
allocate(confbas_final(maxval(nbasis),nelems))
nconfs_out=0
trans_conf=0
inc3=0
conf_final=0

do l=1,mlab_num
   do i=1,nelems
      do_basis: do j=1,nbasis(i)
         if (mlab_atom(final_choice(j,i)) .ne. l) cycle do_basis
         inc2=confnum_all(final_choice(j,i))
         do k=1,inc3
            if (inc2 .eq. conf_final(k)) then
               cycle do_basis
            end if   
         end do
         inc3=inc3+1
         conf_final(inc3)=inc2
      end do do_basis
   end do
end do
nconfs_out=inc3
!
!    Initialize array for translation from new to old configuration numbers
!    (important if not all old configurations will be part of the new ML_AB)
!
do i=1,conf_num
   do j=1,nconfs_out
      if (conf_final(j) .eq. i) then
         trans_conf(j)=i
      end if     
   end do
end do

!
!     Initialize array for reverse translation from old to new configuration
!     numbers, important for list of basis functions and configurations
!
conf_final=0
do i=1,nconfs_out
   do j=1,conf_num
      if (trans_conf(i) .eq. j) then
         conf_final(j)=i
      end if        
   end do
end do


!open(unit=45,file="trans_conf.dat")
!do i=1,nconfs_out
!   write(45,*) i,trans_conf(i)
!end do
!close(45)

!open(unit=45,file="final_choice.dat")
!do i=1,nbasis(1)
!   write(45,*) i,final_choice(i,:)
!end do
!close(45)

!open(unit=45,file="conf_reverse.dat")
!do i=1,conf_num
!   write(45,*) i,conf_final(i)
!end do
!close(45)

!
!     Write local environments of all chosen basis atoms to a trajectory file
!
max_around=maxval(sum(num_around,dim=1))
open(unit=56,file="environments.xyz",status="replace")
do i=1,nelems
   do j=1,nbasis(i)
      write(56,*) max_around
      write(56,'(a,a,a,i8,a)') " ",el_list_glob(i)," basis atom ",j," (environment)"
      do k=1,max_around
         call atomname(ind_env(k,final_choice(j,i)),el_act)
         write(56,*) el_act,environ(:,k,final_choice(j,i))
      end do
   end do
end do
close(56)

write(*,*)
write(*,*) "File 'environments.xyz' with local environments of basis functions written."
write(*,*)
write(34,*)
write(34,*) "File 'environments.xyz' with local environments of basis functions written."
write(34,*)
flush(34)

!
!    First write the header of the new ML_AB file
!
open(unit=56,file="ML_AB_sel",status="replace")
write(56,'(a)') " 1.0 Version (written by mlff_select of VASP4CLINT)"
write(56,'(a)') "**************************************************"
write(56,'(a)') "     The number of configurations"
write(56,'(a)') "--------------------------------------------------"
write(56,'(i10)') nconfs_out
write(56,'(a)') "**************************************************"
write(56,'(a)') "     The maximum number of atom type"
write(56,'(a)') "--------------------------------------------------"
write(56,'(i8)') nelems
write(56,'(a)') "**************************************************"
write(56,'(a)') "     The atom types in the data file"
write(56,'(a)') "--------------------------------------------------"
do i=1,int(nelems/3)
   write(56,'(a)',advance="no") "    "
   do j=1,3
      write(56,'(a,a)',advance="no") " ",el_list_glob((i-1)*3+j)
   end do   
   write(56,*) 
end do
if (int(nelems/3)*3 .lt. nelems) then
   write(56,'(a)',advance="no") "    "     
   do j=int(nelems/3)*3,nelems-1
      write(56,'(a,a)',advance="no") " ",el_list_glob(1+j)
   end do
   write(56,*)
end if   
write(56,'(a)') "**************************************************"
write(56,'(a)') "     The maximum number of atoms per system "
write(56,'(a)') "--------------------------------------------------"
write(56,'(i15)') natoms_max
write(56,'(a)') "**************************************************"
write(56,'(a)') "     The maximum number of atoms per atom type"
write(56,'(a)') "--------------------------------------------------"
write(56,'(i15)') nat_type_max
write(56,'(a)') "**************************************************"
write(56,'(a)') "     Reference atomic energy (eV)"
write(56,'(a)') "--------------------------------------------------"
do i=1,int(nelems/3)
   write(56,'(3e18.9)') 0.d0, 0.d0, 0.d0
end do
if (int(nelems/3)*3 .lt. nelems) then
   do j=int(nelems/3)*3,nelems 
      write(56,'(e18.10)',advance="no") 0.d0
   end do
   write(56,*)
end if
write(56,'(a)') "**************************************************"
write(56,'(a)') "     Atomic mass"
write(56,'(a)') "--------------------------------------------------"
do i=1,int(nelems/3)
   do j=1,3
      write(56,'(f20.12)',advance="no") at_mass_glob((i-1)*3+j)
   end do   
   write(56,*) 
end do
if (int(nelems/3)*3 .lt. nelems) then
   do j=int(nelems/3)*3,nelems-1 
      write(56,'(f20.12)',advance="no") at_mass_glob(1+j)
   end do
   write(56,*)
end if
write(56,'(a)') "**************************************************"
write(56,'(a)') "     The numbers of basis sets per atom type"
write(56,'(a)') "--------------------------------------------------"
do i=1,int(nelems/3)
   write(56,'(a)',advance="no") "    "
   do j=1,3
      write(56,'(i7)',advance="no") nbasis(j)
   end do   
   write(56,*) 
end do
if (int(nelems/3)*3 .lt. nelems) then
   if ((nelems-int(nelems/3)*3) .eq. 3) then
      write(56,'(3i7)') nbasis(1:3)     
   else if ((nelems-int(nelems/3)*3) .eq. 2) then   
      write(56,'(2i7)') nbasis(1:2)
   else if ((nelems-int(nelems/3)*3) .eq. 1) then
      write(56,'(i7)') nbasis(1)
   end if
end if

!
!    Now write the basis set information for each element, separately
!

do i=1,nelems
   write(56,'(a)') "**************************************************"
   write(56,'(a,a)') "     Basis set for ",el_list_glob(i)
   write(56,'(a)') "--------------------------------------------------"
   do j=1,nbasis(i)
      write(56,'(a,3i7)') "    ",conf_final(confnum_all(final_choice(j,i))), &
                       & nat_all(final_choice(j,i))
   end do
end do
write(56,'(a)') "**************************************************"

!
!     Write out the information for all remaining configurations in the set
!

do i=1,nconfs_out
   write(56,'(a,i7)') "     Configuration num. ",i
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     System name"
   write(56,'(a)') "--------------------------------------------------"
   write(56,'(a)') "     Reselected by mlff_select (VASP4CLINT)"
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     The number of atom types"
   write(56,'(a)') "--------------------------------------------------"
   write(56,'(i8)') nelems_confs(trans_conf(i))
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     The number of atoms"
   write(56,'(a)') "--------------------------------------------------"
   write(56,'(i11)') natoms(trans_conf(i))
   write(56,'(a)') "**************************************************"
   write(56,'(a)') "     Atom types and atom numbers"
   write(56,'(a)') "--------------------------------------------------"
   do j=1,nelems_confs(trans_conf(i))
      write(56,'(a,a,i7)') "     ",el_list_confs(j,trans_conf(i)), &
                           & el_nums_confs(j,trans_conf(i))
   end do
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     CTIFOR"
   write(56,'(a)') "--------------------------------------------------"
   write(56,*) ctifors(trans_conf(i))  
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     Primitive lattice vectors (ang.)"
   write(56,'(a)') "--------------------------------------------------"
   do j=1,3
      write(56,*) cells(:,j,trans_conf(i))
   end do
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     Atomic positions (ang.)"
   write(56,'(a)') "--------------------------------------------------"
   do j=1,natoms(trans_conf(i))
      write(56,*) xyz(:,j,trans_conf(i))
   end do
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     Total energy (eV)"
   write(56,'(a)') "--------------------------------------------------"
   write(56,*) energies(trans_conf(i))
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     Forces (eV ang.^-1)"
   write(56,'(a)') "--------------------------------------------------"
   do j=1,natoms(trans_conf(i))
      write(56,*) grads(:,j,trans_conf(i))
   end do
   write(56,'(a)') "=================================================="
   write(56,'(a)') "     Stress (kbar)"
   write(56,'(a)') "--------------------------------------------------"
   write(56,'(a)') "     XX YY ZZ"
   write(56,'(a)') "--------------------------------------------------"
   write(56,*) stress(1:3,trans_conf(i))
   write(56,'(a)') "--------------------------------------------------"
   write(56,'(a)') "     XY YZ ZX"
   write(56,'(a)') "--------------------------------------------------"
   write(56,*) stress(4:6,trans_conf(i))
   write(56,'(a)') "**************************************************"
end do
close(56)

call cpu_time(time2)

write(*,*) "New ML_AB file written to file 'ML_AB_sel'."
write(*,*)
write(*,'(a,f12.4,a)') " Total execution time: ",time2-time1," seconds"
write(*,*) "mlff_select finished normally."
write(*,*)
write(34,*) "New ML_AB file written to file 'ML_AB_sel'."
write(34,*)
write(34,'(a,f12.4,a)') " Total execution time: ",time2-time1," seconds"
write(34,*) "mlff_select finished normally."
write(34,*)
!
!     Close logfile for alternative output to file
!
close(34)

end program mlff_select

!
!     subroutine elem: read in character with element symbol and
!       give out the number
!
!     part of QMDFF (taken from Caracal)
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

!
!     subroutine atomname: is an atomic charge is given, extract the
!        corresponding element label
!
!     part of EVB
!
subroutine atomname(nat,key1)
implicit none
CHARACTER(len=2)::KEY1  ! the name of the element
CHARACTER(len=2)::ELEMNT(107)  ! array with all names
integer::nat  ! atomic charge

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

!
!     Return the i'th entry of the elements array of i is the given
!     atomic charge of the atom
!
key1=elemnt(nat)


return
end subroutine atomname



!
!    subroutine matinv3: inversion of 3x3 matrices
!    modified from https://fortranwiki.org/fortran/show/Matrix+inversion
!
subroutine matinv3(A,B)
!! Performs a direct calculation of the inverse of a 3×3 matrix.
real(kind=8), intent(in) :: A(3,3)   !! Matrix
real(kind=8), intent(out)   :: B(3,3)   !! Inverse matrix
real(kind=8)             :: detinv

! Calculate the inverse determinant of the matrix
detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

! Calculate the inverse of the matrix
B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

return
end subroutine
