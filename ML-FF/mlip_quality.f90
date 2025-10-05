!
!    mlip_quality: Evaluate the quality of a MLP with respect to the 
!      reproduction of a given training or verification set.
!      Currently, Behler neural networks, message-passing atomic cluster
!      expansion (MACE) and VASP ML-FF production runs can be evaluated.
!
!    Part of utils4VASP
!
!     Julien Steffen,     2025 (julien.steffen@fau.de)
!     Maximilian Bechtel, 2025 (maxi.bechtel@fau.de)
!

program mlp_quality

implicit none
integer::readstat
integer::i,j,k,l
integer::nframes
integer::natoms_max,natoms_act
integer::ngrads
integer,allocatable::natoms_list(:)
integer::frame_act
integer::force_num
integer::nhisto_en,nhisto_grad
integer::nhisto_2d_abs,nhisto_2d_angle
integer,allocatable::histo_en(:),histo_grad(:)
integer,allocatable::histo_grad_2d(:,:)
real(kind=8)::nhisto_en_range,nhisto_grad_range
real(kind=8)::nhisto_2d_abs_range,nhisto_2d_angle_range
real(kind=8)::fdum,ediff,gdiff
real(kind=8)::angle,cos_theta
real(kind=8)::max_histo
real(kind=8)::gradnorm_ann,gradnorm_ref
real(kind=8),allocatable::energies_ann(:),energies_ref(:)
real(kind=8),allocatable::gradients_ann(:,:,:),gradients_ref(:,:,:)
real(kind=8)::mae_energy,rmse_energy,mae_forces,rmse_forces
character(len=120)::a120,a180,adum,arg
logical,allocatable::select_all(:,:,:)
character(len=1)::select_line(3)
character(len=80)::ref_filename,a80
character(len=80),allocatable::filename_list(:)
logical::eval_ann,eval_mace,eval_vasp

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
write(*,*) "PROGRAM mlip_quality: Evaluation of a MLP with respect to the "
write(*,*) " reproduction of a given training or verification set."
write(*,*) "The analysis can be done for Behler neural network potentials "
write(*,*) " optimized with the aenet program (ANN) or for message-passing "
write(*,*) " atomic cluster expansion (MACE) potentials."
write(*,*) "For ANN, the predict.x program must have been executed in the "
write(*,*) " current folder, such that the predict.out file and a folder"
write(*,*) " 'xsf_files' are present."
write(*,*) "For MACE, the file vasp_results.xyz must be there, generated"
write(*,*) " by vasp2trainset and containing the VASP energies and forces"
write(*,*) " Further, the files energies_mace.dat, resulting from the ASE"
write(*,*) " MACE trajectory (without header lines) and gradients_mace.dat"
write(*,*) " also resulting from the ASE calculation, must be there."
write(*,*) "For VASP ML-FF, the file vasp_results.xyz must be there, generated"
write(*,*) " by vasp2trainset containing the VASP energies and forces"
write(*,*) " (the same as for MACE). Further, the OUTCAR from the VASP ML-FF"
write(*,*) " production run must be there."
write(*,*) "The ordering of the structures must be the same in all files!"
write(*,*) "The following keywords can be given (those with default values"
write(*,*) " are optional):"
write(*,*) " -overview:  print an overview of all scripts and programs "
write(*,*) "   in utils4VASP"
write(*,*) " -ann  : ANN        output will be evaluated."
write(*,*) " -mace : MACE       output will be evaluated."
write(*,*) " -vasp : VASP ML-FF output will be evaluated."
write(*,*) " -natoms_max=[number] : Maximum number of atoms per frame of the "
write(*,*) "   training set (default: 1000)"
write(*,*) " -nhisto_en=[number] : Number of histogram bins for plot of "
write(*,*) "   errors in the energy deviation per atom (default: 100) "
write(*,*) " -nhisto_grad=[number] : Number of histogram bins for plot of "
write(*,*) "   errors in the energy deviation per atom (default: 400) "
write(*,*) " -nhisto_2d_abs=[number] : Number of histogram bins for the 2D"
write(*,*) "   plot of angle deviations, the absolute forces (default: 100)"
write(*,*) " -nhisto_2d_angle=[number] : Number of histogram bins for the 2D"
write(*,*) "   plot of angle deviations, the angle deviations (default: 100)"
write(*,*) " -histo_en_range=[value] : Energy range for energy deviation "
write(*,*) "   histogram, in meV per atom (default: 10.0)"
write(*,*) " -histo_grad_range=[value] : Gradient range for gradient deviation"
write(*,*) "   histogram, in meV/Ang (default: 100.0)"
write(*,*) " -histo_2d_abs_range=[value] : Range of absolute gradient values"
write(*,*) "   for the 2D gradient histogram, in eV/Ang (default: 3.0)"
write(*,*) " -histo_2d_angle_range=[value] : Range of direction error for the 2D"
write(*,*) "   2D gradient histogram, in degrees (default: 60.0)"
write(*,*) 

!
!     If ANN data shall be evaluated
!
eval_ann=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-ann") then
      eval_ann=.true.
   end if
end do


!
!     If MACE data shall be evaluated
!
eval_mace=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-mace") then
      eval_mace=.true.
   end if
end do

!
!     If VASP ML-FF data shall be evaluated
!
eval_vasp=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-vasp") then
      eval_vasp=.true.
   end if
end do

if (.not. eval_ann .and. .not. eval_mace .and. .not. eval_vasp) then
   write(*,*) "Please choose either ANN, MACE or VASP ML-FF!"
   stop
end if
if (eval_ann .and. eval_mace) then
   write(*,*) "Please choose either ANN or MACE!"
   stop
end if
if (eval_ann .and. eval_vasp) then
   write(*,*) "Please choose either ANN or VASP ML-FF!"
   stop
end if
if (eval_mace .and. eval_vasp) then
   write(*,*) "Please choose either MACE or VASP ML-FF!"
   stop
end if
if (eval_ann .and. eval_mace .and. eval_vasp) then
   write(*,*) "Please choose either ANN, MACE or VASP ML-FF!"
   stop
end if
if (eval_ann) then
   write(*,*) "Artificial Behler Neural Network (ANN) will be evaluated!"
end if
if (eval_mace) then
   write(*,*) "Message-Passing Atomic Cluster Expansion (MACE) will be evaluated!"
end if   
if (eval_vasp) then
   write(*,*) "Energies and gradients from a VASP ML-FF OUTCAR will be evaluated ..."
end if
write(*,*)
!
!     Maximum number of atoms that can appear in the training set
!
natoms_max=1000
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-natoms_max=") then
      read(arg(13:),*,iostat=readstat) natoms_max
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -natoms_max=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Number of histogram points for energies and gradients
!
nhisto_en=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-nhisto_en=") then
      read(arg(12:),*,iostat=readstat) nhisto_en
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -nhisto_en=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_grad=400
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-nhisto_grad=") then
      read(arg(14:),*,iostat=readstat) nhisto_grad
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -nhisto_grad=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_2d_abs=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:15))  .eq. "-nhisto_2d_abs=") then
      read(arg(16:),*,iostat=readstat) nhisto_2d_abs
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -nhisto_2d_abs=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_2d_angle=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:17))  .eq. "-nhisto_2d_angle=") then
      read(arg(18:),*,iostat=readstat) nhisto_2d_abs
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -nhisto_2d_angle=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Energy range for histogram
!
nhisto_en_range=10.0d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:16))  .eq. "-histo_en_range=") then
      read(arg(17:),*,iostat=readstat) nhisto_en_range
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -histo_en_range=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Gradient range for histogram
!
nhisto_grad_range=100.0d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:18))  .eq. "-histo_grad_range=") then
      read(arg(19:),*,iostat=readstat) nhisto_grad_range
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -histo_grad_range=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

!
!     Ranges for 2D plots of gradient angles and absolute values
!
nhisto_2d_abs_range=3.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:20))  .eq. "-histo_2d_abs_range=") then
      read(arg(21:),*,iostat=readstat) nhisto_2d_abs_range
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -histo_2d_abs_range=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_2d_angle_range=60d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:22))  .eq. "-histo_2d_angle_range=") then
      read(arg(23:),*,iostat=readstat) nhisto_2d_angle_range
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -histo_2d_angle_range=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Histogram arrays
!

allocate(histo_en(nhisto_en))
allocate(histo_grad(nhisto_grad+10))
allocate(histo_grad_2d(nhisto_2d_abs,nhisto_2d_angle))

!
!     Evaluate ANN results 
!
if (eval_ann) then
!
!     Range of histogram plot 
!
   open(unit=56,file="predict.out",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file 'predict.out' is not there!"
      stop
   end if
   do
      read(56,'(a)',iostat=readstat) a120
      if (readstat .ne. 0) exit

      if (index(a120,'Number of structures in the data set:') .ne. 0 ) then
         read(a120,*) adum,adum,adum,adum,adum,adum,adum,nframes
         write(*,*) "Number of structures to evaluate: ",nframes
         exit
      end if 
   end do
   write(*,*) 
   write(*,*) "Start quality evaluation ..."
   close(56)
   allocate(natoms_list(nframes))
   allocate(energies_ann(nframes),energies_ref(nframes))
   allocate(gradients_ann(3,natoms_max,nframes))
   allocate(gradients_ref(3,natoms_max,nframes))
   allocate(select_all(3,natoms_max,nframes))
   select_all=.true.
   allocate(filename_list(nframes))
!
!    Loop through all evaluated training set pieces
!
   open(unit=56,file="predict.out",status="old")
   frame_act=0
   do
      read(56,'(a)',iostat=readstat) a80
      if (readstat .ne. 0) exit
      if (index(a80,'Structure number  :') .ne. 0 ) then
!
!    Read in ANN results
!
         frame_act=frame_act+1
         read(56,'(a)') a120
         ref_filename=a120(22:)
         filename_list(frame_act)=ref_filename
         read(56,*) adum,adum,adum,adum,natoms_act
         natoms_list(frame_act)=natoms_act
         do i=1,13
            read(56,'(a)') adum
         end do
         do i=1,natoms_act
            read(56,'(a)') a180
            read(a180,*,iostat=readstat) adum,fdum,fdum,fdum,gradients_ann(:,i,frame_act)
         end do
         read(56,*)
         read(56,*)
         read(56,*) adum,adum,adum,energies_ann(frame_act)
!
!    Read in reference results
!
         open(unit=57,file=ref_filename,status="old",iostat=readstat)
         if (readstat .ne. 0) then
            write(*,*) "The reference file ",trim(ref_filename)," is not there!"
            stop
         end if
         read(57,*) adum,adum,adum,adum,energies_ref(frame_act)
         do i=1,8
            read(57,'(a)') adum
         end do
         do i=1,natoms_act
            read(57,'(a)') a180
            read(a180,*,iostat=readstat) adum,fdum,fdum,fdum,gradients_ref(:,i,frame_act),select_line
            if (readstat .eq. 0) then
               do j=1,3
                  if (select_line(j) .eq. "T") select_all(j,i,frame_act)=.true.
                  if (select_line(j) .eq. "F") select_all(j,i,frame_act)=.false.
               end do
            end if

         end do    
         close(57)
      end if
   end do
   close(56)
   write(*,*) "... finished!"
end if


!
!     Evaluate MACE results 
!
if (eval_mace) then
!
!     Range of histogram plot 
!
   open(unit=56,file="vasp_results.xyz",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file 'vasp_results.xyz' is not there!"
      stop
   end if
   nframes=0
   do
      read(56,'(a)',iostat=readstat) a120
      if (readstat .ne. 0) exit

      if (index(a120,'Lattice=') .ne. 0 ) then
         nframes=nframes+1
      end if
   end do
   write(*,*) "Number of structures to evaluate: ",nframes
   write(*,*)
   write(*,*) "Start quality evaluation ..."
   close(56)
   allocate(natoms_list(nframes))
   allocate(energies_ann(nframes),energies_ref(nframes))
   allocate(gradients_ann(3,natoms_max,nframes))
   allocate(gradients_ref(3,natoms_max,nframes))
   allocate(select_all(3,natoms_max,nframes))
   allocate(filename_list(nframes))
   filename_list="MACE"
   select_all=.true.
!
!    Loop through all evaluated training set pieces
!
   open(unit=56,file="vasp_results.xyz",status="old")
   open(unit=57,file="gradients_mace.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file 'gradients_mace.dat' is not there!"
      stop
   end if
   open(unit=58,file="energies_mace.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file 'energies_mace.dat' is not there!"
      stop
   end if

   frame_act=0
   do i=1,nframes
      read(56,*) natoms_act
      natoms_list(i)=natoms_act
      read(56,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,energies_ref(i)
      do j=1,natoms_act
         read(56,*,iostat=readstat) adum,fdum,fdum,fdum,fdum,gradients_ref(:,j,i)
      end do
      read(58,*) adum,energies_ann(i)
      read(57,*)
      do j=1,natoms_act
         read(57,*,iostat=readstat) gradients_ann(:,j,i)
      end do
   end do
   close(56)
   close(57)
   close(58)
   write(*,*) "... finished!"
end if
write(*,*)

!
!     Evaluate VASP ML-FF results 
!
if (eval_vasp) then

   open(unit=56,file="vasp_results.xyz",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file 'vasp_results.xyz' is not there!"
      stop
   end if

   !
   ! Determine the number of frames in vasp_results.xyz
   !
   nframes = 0
   do
      read(56,'(a)',iostat=readstat) a120
      if (readstat .ne. 0) exit  ! Break if End of File reached

      if (index(a120,'Lattice=') .ne. 0 ) then  ! "Lattice=" marks begin of a new frame
         nframes = nframes + 1
      end if
   end do

   write(*,*) "Number of structures to evaluate: ", nframes
   write(*,*)
   write(*,*) "Start quality evaluation ..."

   close(56)

   !
   ! Allocate all needed arrays
   !
   allocate(natoms_list(nframes))
   allocate(energies_ann(nframes),energies_ref(nframes))
   allocate(gradients_ann(3,natoms_max,nframes))
   allocate(gradients_ref(3,natoms_max,nframes))
   allocate(select_all(3,natoms_max,nframes))
   allocate(filename_list(nframes))
   filename_list="VASP ML-FF"
   select_all=.true.

   !
   ! Loop through all evaluated training set pieces
   !
   open(unit=56,file="vasp_results.xyz",status="old")
   open(unit=57,file="OUTCAR",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file 'OUTCAR' is not there!"
      stop
   end if

   ! REFERENCE DATA !
   write(*,*) "Collecting reference data from vasp_results.xyz ..."
   frame_act = 0
   do i = 1, nframes
      read(56,*) natoms_act
      natoms_list(i) = natoms_act  ! Get current number of atoms

      ! Get reference energy
      read(56,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,adum,energies_ref(i)

      ! Get reference gradients
      do j = 1, natoms_act
         read(56,*,iostat=readstat) adum,fdum,fdum,fdum,fdum,gradients_ref(:,j,i)
      end do
   end do

   ! VASP ML-FF DATA FROM OUTCAR !
   write(*,*) "Collecting VASP ML-FF data from OUTCAR ..."
   i = 1
   do
      ! Loop through OUTCAR until End of File reached
      read(57,'(a)',iostat=readstat) a120
      if (readstat .ne. 0) exit

      ! Check if line with ML gradients was found
      if ( index(a120, "TOTAL-FORCE (eV/Angst) (ML)") .ne. 0 ) then
        read(57,'(a)') adum

        do j = 1, natoms_list(i)
           read(57,*,iostat=readstat) adum,adum,adum,gradients_ann(:,j,i)
        end do
      end if

      ! Check if line with ML energy was found
      if ( index(a120, "ML energy(sigma->0)") .ne. 0 ) then
         read(a120,*,iostat=readstat) adum,adum,adum,adum,adum,adum,adum,adum,energies_ann(i)
         i = i + 1
      end if
   end do

   close(56)
   close(57)
   write(*,*) "... finished!"

end if
write(*,*)

!
!     Calculate mean absolute errors (MAE) and root mean square errors (RMSE)
!
mae_energy=0d0
rmse_energy=0d0
do i=1,nframes
   mae_energy=mae_energy+abs(energies_ref(i)-energies_ann(i))   
   rmse_energy=rmse_energy+(energies_ref(i)-energies_ann(i))**2
end do
mae_energy=mae_energy/sum(natoms_list)*1000d0
rmse_energy=sqrt(rmse_energy/sum(natoms_list))*1000d0

write(*,*) "---------------------------------------------------------"
write(*,'(a,f13.6,a)') "MAE(energy):",mae_energy," meV/atom"
write(*,'(a,f13.6,a)') "RMSE(energy):",rmse_energy," meV/atom"

!
!     MAE and RMSE for forces, only for selected atoms!
!
mae_forces=0d0
rmse_forces=0d0
force_num=0
do i=1,nframes
   do j=1,natoms_list(i)
      do k=1,3
         if (select_all(1,j,i)) then
            mae_forces=mae_forces+abs(gradients_ref(k,j,i)-gradients_ann(k,j,i))
            rmse_forces=rmse_forces+(gradients_ref(k,j,i)-gradients_ann(k,j,i))**2
            force_num=force_num+1
         end if
      end do
   end do
end do
mae_forces=mae_forces/real(force_num)*1000d0
rmse_forces=sqrt(rmse_forces/real(force_num))*1000d0

write(*,'(a,f13.6,a)') "MAE(force components):",mae_forces," meV/Ang."
write(*,'(a,f13.6,a)') "RMSE(force components):",rmse_forces," meV/Ang."
write(*,*) "---------------------------------------------------------"
write(*,*)

!
!     Generate scattering plot for the energies
!

open(unit=60,file="energies_compare.dat",status="replace")
do i=1,nframes
   write(60,*) energies_ref(i),energies_ann(i)," # ",trim(filename_list(i))
end do
close(60)
write(*,*) "Energy scattering plot written to 'energies_compare.dat'"
!
!     Generate scattering plot for the gradients
!

open(unit=60,file="gradients_compare.dat",status="replace")
do i=1,nframes
   do j=1,natoms_list(i)
      if (select_all(1,j,i)) then
         gradnorm_ann=sqrt(gradients_ann(1,j,i)**2+gradients_ann(2,j,i)**2+gradients_ann(3,j,i)**2)
         gradnorm_ref=sqrt(gradients_ref(1,j,i)**2+gradients_ref(2,j,i)**2+gradients_ref(3,j,i)**2)
         write(60,*) gradnorm_ann,gradnorm_ref
      end if
   end do
end do
close(60)
write(*,*) "Gradient scattering plot written to 'gradients_compare.dat'"

!
!    Generate histogram plots for energies and gradients
!
histo_en=0
do i=1,nframes
   ediff=abs(energies_ref(i)-energies_ann(i))/natoms_list(i)*1000
   j=nint(ediff/nhisto_en_range*nhisto_en)+1
   if (j .gt. nhisto_en) cycle
   histo_en(j)=histo_en(j)+1
end do

open(unit=60,file="energies_histo.dat",status="replace")
do i=1,nhisto_en
   write(60,*) real(i)/real(nhisto_en)*nhisto_en_range,real(histo_en(i))/real(nframes)
end do
close(60)
write(*,*) "Energy error histogram written to 'energies_histo.dat'"

histo_grad=0
ngrads=0
do i=1,nframes
   do j=1,natoms_list(i)
      if (select_all(1,j,i)) then
         gradnorm_ann=sqrt(gradients_ann(1,j,i)**2+gradients_ann(2,j,i)**2+gradients_ann(3,j,i)**2)
         gradnorm_ref=sqrt(gradients_ref(1,j,i)**2+gradients_ref(2,j,i)**2+gradients_ref(3,j,i)**2)
         gdiff=abs(gradnorm_ann-gradnorm_ref)*1000
         k=nint(gdiff/nhisto_grad_range*nhisto_grad)+1
         if (k .gt. nhisto_grad) cycle
         histo_grad(k)=histo_grad(k)+1
         ngrads=ngrads+1
      end if
   end do
end do


open(unit=60,file="gradients_histo.dat",status="replace")
do i=1,nhisto_grad
   write(60,*) real(i)/real(nhisto_grad)*nhisto_grad_range,real(histo_grad(i))/real(ngrads)
end do
close(60)
write(*,*) "Gradient error histogram written to 'gradients_histo.dat'"

!
!    Generate 2D density plot for gradient angles
!
ngrads=0
histo_grad_2d=0
do i=1,nframes
   do j=1,natoms_list(i)
      if (select_all(1,j,i)) then
 !     write(*,*) gradients_ref(:,j,i),gradients_ann(:,j,i)
         gradnorm_ann=sqrt(gradients_ann(1,j,i)**2+gradients_ann(2,j,i)**2+gradients_ann(3,j,i)**2)
         gradnorm_ref=sqrt(gradients_ref(1,j,i)**2+gradients_ref(2,j,i)**2+gradients_ref(3,j,i)**2)
         cos_theta=dot_product(gradients_ref(:,j,i),gradients_ann(:,j,i))/gradnorm_ann/gradnorm_ref
!      write(*,*) cos_theta
         if (cos_theta .lt. -1.d0) cos_theta=-1.d0
         if (cos_theta .gt. 1.d0) cos_theta=1.d0
         angle=acos(cos_theta)
         angle=abs(angle*180.d0/acos(-1.d0))

         k=nint(angle/nhisto_2d_angle_range*nhisto_2d_angle)+1

         l=nint(gradnorm_ref/nhisto_2d_abs_range*nhisto_2d_abs)+1

         if (k .gt. nhisto_2d_angle .or. l .gt. nhisto_2d_abs) cycle
         histo_grad_2d(l,k)=histo_grad_2d(l,k)+1
         ngrads=ngrads+1
      end if
   end do
end do

max_histo=maxval(histo_grad_2d)


open(unit=60,file="gradnorm_vs_angle_histo.dat",status="replace")
do i=1,nhisto_2d_abs
   do j=1,nhisto_2d_angle
      write(60,*) real(i-1)/real(nhisto_2d_abs)*nhisto_2d_abs_range,real(j-1)/real(nhisto_2d_angle)* &
              & nhisto_2d_angle_range,real(histo_grad_2d(i,j))/max_histo
   end do
   write(60,*)
end do
close(60)
write(*,*) "2D gradient direction histogram written to 'gradnorm_vs_angle_histo.dat'"
write(*,*)
!
!     Generate gnuplot files for the different plots and execute the gnuplot commands!
!
!     A) The energy per atom histogram
!
write(*,*) "Execute gnuplot file 'plot_energy_histo.gnu'"
open(unit=61,file="plot_energy_histo.gnu",status="replace")
write(61,*) "set terminal svg lw 1.5 size 400,350"
write(61,*) "set output 'energies_histo.svg'"
write(61,*) "set xlabel 'error in energy (meV/atom)'"
write(61,*) "set ylabel 'relative frequency'"
write(61,*) "set grid"
write(61,*) "set style fill solid"
write(61,*) "plot 'energies_histo.dat' u 1:2 w boxes lc rgb 'dark-red' title 'validation set'"
close(61)
call system("gnuplot plot_energy_histo.gnu")
write(*,*) "Plot picture written to 'energies_histo.svg'!"
!
!     B) The force per atom histogram
!
write(*,*) "Execute gnuplot file 'plot_gradient_histo.gnu'"
open(unit=61,file="plot_gradient_histo.gnu",status="replace")
write(61,*) "set encoding iso_8859_1"
write(61,*) "set terminal svg lw 1.5 size 400,350"
write(61,*) "set output 'gradients_histo.svg'"
write(61,*) "set xlabel 'error in absolute force (meV/{\305})'"
write(61,*) "set ylabel 'relative frequency'"
write(61,*) "set grid"
write(61,*) "set style fill solid"
write(61,*) "plot 'gradients_histo.dat' u 1:2 w boxes lc rgb 'dark-red' title 'validation set'"
close(61)
call system("gnuplot plot_gradient_histo.gnu")
write(*,*) "Plot picture written to 'gradients_histo.svg'!"
!
!    C) The 2D force direction error histogram
!
write(*,*) "Execute gnuplot file 'plot_2d_grad_histo.gnu'"
open(unit=61,file="plot_2d_grad_histo.gnu",status="replace")
write(61,*) "set encoding iso_8859_1"
write(61,*) "set terminal png lw 4.5 size 2000,2000 font 'Helvetica,48'"
write(61,*) "set output '2d_gradient_histo.png'"
write(61,*) "set xlabel 'absolute atomic force (eV/{\305})'"
write(61,*) "set ylabel 'error in direction (°)'"
write(61,*) "set zlabel 'relative frequency'"
write(61,*) "unset key"
write(61,*) "set pm3d map interpolate 5,5"
write(61,*) "set palette defined ( 0 'white',\"
write(61,*) "     0.01 'black',\"
write(61,*) "     0.15 'blue',\"
write(61,*) "     0.4 'green',\"
write(61,*) "     0.7 'yellow',\"
write(61,*) "     1 'red' )"
write(61,*) "splot 'gradnorm_vs_angle_histo.dat' u 1:2:3 with pm3d"
close(61)
call system("gnuplot plot_2d_grad_histo.gnu")
write(*,*) "Plot picture written to '2d_gradient_histo.png'!"

write(*,*) 
write(*,*) "Program mlp_quality exited normally."
write(*,*)

end program mlp_quality

