!
!    modify_xdatcar: Perform different tasks on a XDATCAR file,
!      such as translating or multiplying the unit cell, picking
!      certain structures or write an xyz trajectory file
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!


program modify_xdatcar
implicit none 
integer::i,j,k,l,m
integer::readstat,openstat,counter,endl
integer::natoms,nelems,xdat_lines,nframes
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
real(kind=8),allocatable::a_vecs(:,:),b_vecs(:,:),c_vecs(:,:)
character(len=2),allocatable::el_names_read(:),el_names(:)
character(len=2),allocatable::at_names(:),at_names2(:)
real(kind=8),allocatable::at_masses(:)
real(kind=8),allocatable::vel_act(:,:),vel_old(:,:)
real(kind=8)::totmass
integer,allocatable::el_nums(:),el_nums_tmp(:)
real(kind=8)::act_num(3),factor
real(kind=8),allocatable::xyz(:,:,:),xyz2(:,:,:)
real(kind=8)::xlen,ylen,zlen
real(kind=8)::shift_vec(3),act_val,xyz_print(3)
integer::multiply_vec(3),pick_ind,pos_new,multiply_prod
integer::frame_first,frame_last,line_num
integer::read_freq,status
logical::npt,print_npt
logical::remove_mode
logical::remove_transrot
!  For removal of overall translation and rotation
real(kind=8)::vtot(3),mang(3),eps,weigh,xdel,ydel,zdel
real(kind=8)::xx,xy,xz,yy,yz,zz,xtot,ytot,ztot
real(kind=8)::tensor(3,3),vang(3)
character(len=40)::remove_com
character(len=1)::remove_dim
character(len=2)::remove_sign
real(kind=8)::com_act(3)
real(kind=8)::remove_border
logical,allocatable::keep_atoms(:,:)
logical::smooth_mode
logical::eval_stat(10)
logical::shift_cell,multiply_cell,pick_frame,print_xyz,print_last
character(len=320)::a120,cdum,arg,adum
character(len=220)::a220
character(len=50)::atest
integer::alloc_stat
character(len=100)::alloc_err
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
write(*,*) "PROGRAM modify_xdatcar: Modification of XDATCAR trajectory"
write(*,*) " files for arbitrary systems."
write(*,*) "The file XDATCAR must be present."
write(*,*) "The following calculations can be done, called by one of the"
write(*,*) " listed command line arguments:"
write(*,*) " -overview:  print an overview of all scripts and programs in "
write(*,*) "    utils4VASP"
write(*,*) " -shift=a,b,c : Shift each frame of the trajectory by a vector"
write(*,*) "    given in direct coordinates. Example: -shift=0.1,0.0,0.2"
write(*,*) " -multiply=a,b,c : Multiply the unit cell of each frame by some"
write(*,*) "    replications in each of the coordinate directions. Integers"
write(*,*) "    must be given as arguments. Example: -multiply=2,2,1"
write(*,*) " -remove=[value] : Remove all atoms in all frames that fulfill"
write(*,*) "    the criteria. Example: -remove=xgt8 (all atoms removed that"
write(*,*) "    have an x value of 8 or larger. The trajectory will then be"
write(*,*) "    written in NpT format. Attention with variable atom number!"
write(*,*) " -read_freq=[number]: Only read in and process every nth frame"
write(*,*) "    Example: -read_freq=10 : every 10th frame will be read in and "
write(*,*) "    processed and written out. (DEFAULT: 1)"  
write(*,*) " -pick_frame=[number] : Pick one frame of the XDATCAR file and print"
write(*,*) "    it to a POSCAR file (POSCAR_pick). Example: -pick_frame=283"
write(*,*) " -print_xyz : Print each frame to a xyz trajectory: xdat_mod.xyz"
write(*,*) " -frame_last=[number] : The number of the last frame to be printed "
write(*,*) "   (default: last frame of the trajectory)"
write(*,*) " -frame_first=[number] : The number of the first frame to be printed "
write(*,*) "   (default: frame No. 1)"
write(*,*) " -remove_transrot : Remove all translation and rotation of the "
write(*,*) "   system by placing each frame at the center of the unit cell "
write(*,*) "   (center of mass) and aligned to principial moments of inertia."
write(*,*) "   Should only be applied to isolated molecular systems!"
write(*,*) " -smooth : Add this to remove all sudden jumps of atoms between "
write(*,*) "   frames due to image flags (e.g. for VDOS calculations)."
write(*,*) " -print_npt : Print a NVT trajectory in the NPT format, e.g., for "
write(*,*) "    subsequent analysis with TRAVIS."
write(*,*) "One or several of these commands can be executed. The ordering"
write(*,*) " will be the same as the listing of the commands."
write(*,*) "If one job apart pick_frame is ordered, the trajectory will be "
write(*,*) " written to XDATCAR_mod or xdat_mod.xyz"
!
!    PART A: Read in the command line arguments !!!!!!!!!!!!!!!!
!
!    The coordinate shift vector
!
shift_cell=.false.
shift_vec=0.0d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-shift=") then
      read(arg(8:),*,iostat=readstat) shift_vec
      shift_cell=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -shift=..., something went wrong!"
         write(*,*)
      end if        
   end if
end do
!
!    The multiplication of unit cells
!
multiply_cell=.false.
multiply_vec=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-multiply=") then
      read(arg(11:),*,iostat=readstat) multiply_vec
      multiply_cell=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -multiply=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Remove all atoms in a certain coordinate range
!
remove_mode=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-remove=") then
      read(arg(9:),*,iostat=readstat) remove_com
      remove_mode=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -remove=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (remove_mode) then
   read(remove_com(:1),*) remove_dim
   read(remove_com(2:3),*) remove_sign
   read(remove_com(4:),*) remove_border
end if        

!
!     Read in the read in frequency (every n'th frame will be 
!          read in and processed)
!
read_freq=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-read_freq=") then
      read(arg(12:),*,iostat=readstat) read_freq
      pick_frame=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -read_freq=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (read_freq .lt. 1) then
   write(*,*) "The read-in frequency must be at least 1!"
   stop
end if        
!
!     Pick a certain frame from the trajectory
!
pick_frame=.false.
pick_ind=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-pick_frame=") then
      read(arg(13:),*,iostat=readstat) pick_ind
      pick_frame=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -pick_frame=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Write the trajectory to a xyz file
!
print_xyz=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-print_xyz") then
      print_xyz=.true.
   end if
end do

!
!     Write the trajectory in the NPT format 
!
print_npt=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-print_npt") then
      print_npt=.true.
   end if
end do
!
!     Remove all translation and rotation of the system by placing 
!     its center of mass in the center of the unit cell, aligned 
!     along its principial moments of inertia
!
remove_transrot=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:16))  .eq. "-remove_transrot") then
      remove_transrot=.true.
   end if
end do

!
!     Remove sudden jumps of atoms between frames: smooth mode
!
smooth_mode=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-smooth") then
      smooth_mode=.true.
   end if
end do

!
!    Activate printout in NPT format if the total number of atoms 
!      changes for each frame after atom deletions
!
if (remove_mode) then
   print_npt=.true.
end if
!
!     Only print frames until frame No frame_last
!
print_last=.false.
frame_last=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-frame_last=") then
      read(arg(13:),*,iostat=readstat) frame_last
      print_last=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -frame_last=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

!
!     Start printing from frame No. frame_first
!
frame_first=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-frame_first=") then
      read(arg(14:),*,iostat=readstat) frame_first
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -frame_first=..., something went wrong!"
         write(*,*)
      end if
   end if
end do



if ((.not. shift_cell) .and. (.not. multiply_cell) .and. (.not. pick_frame) .and. &
           &  (.not. print_xyz) .and. (.not. print_last) .and. (read_freq .eq. 1) .and. &
           &  (.not. print_npt) .and. (.not. remove_mode) .and. (.not. remove_transrot)) then
   write(*,*)
   write(*,*) "Please give at least one of the possible commands!"
   write(*,*)
   stop
end if

!
!    PART B: Read in the XDATCAR file !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    First, determine the number of lines in the XDATCAR file
!
write(*,*)
write(*,*) "Determine number of atoms and frames in XDATCAR..."
status = system("wc -l XDATCAR > xdat_length")
if (status .ne. 0) then
   write(*,*) "The number of atoms/frames could not be determined! XDATCAR missing?"
   stop
end if
open(unit=45,file="xdat_length",status="old")
read(45,*) xdat_lines
close(45)
write(*,*) "completed!"

open(unit=14,file="XDATCAR",status="old",iostat=openstat)
if (openstat .ne. 0) then
   stop "ERROR! The file 'XDATCAR' could not been found!"
end if
read(14,*)
read(14,*) factor   ! the global geometry conversion factor
!
!    Read in the cell dimensions
!
read(14,*) a_vec(1),a_vec(2),a_vec(3)
read(14,*) b_vec(1),b_vec(2),b_vec(3)
read(14,*) c_vec(1),c_vec(2),c_vec(3)
!
!    Read in the elements
!
allocate(el_names_read(50))
el_names_read="XX"
read(14,'(a)') cdum
read(cdum,*,iostat=readstat) el_names_read
nelems=0
do i=1,50
   if (el_names_read(i) .eq. "XX") exit
   nelems=nelems+1
end do
allocate(el_names(nelems),el_nums(nelems))
allocate(el_nums_tmp(nelems))

do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(14,*) el_nums
!
!    Define the element symbols for all atoms in the system
!
natoms = sum(el_nums)
allocate(at_names(natoms))

!
!    Check if the trajectory is an NVT or NPT trajectory
!
read(14,*)
do i=1,natoms
   read(14,*)
end do
read(14,*) adum
!
!    How does the first line of the second frame begin? If it begins 
!    with Direct or direct, the trajectory is a NVT trajectory
!
if (trim(adum) .eq. "Direct" .or. trim(adum) .eq. "direct") then
   npt=.false.
   write(*,*) "A NVT trajectory (constant unit cell) has been detected."
else
   npt=.true.
   write(*,*) "A NpT trajectory (varying unit cell) has been detected."
end if

counter = 1
do i=1,nelems
   do j=1,el_nums(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do
if (npt) then
   nframes = (xdat_lines - 7)/(natoms+8)
else
   nframes = (xdat_lines - 7)/(natoms+1)
end if
!write(*,*) xdat_lines,natoms,nframes
!
!    If the read frequency is larger than 1, divide the number of 
!       frames by it!
!
if (read_freq .gt. 1) then
   nframes=int(nframes/read_freq)
end if        
!
!    From the print-last command, determine the first frame to be written
!
if (print_last) then
   frame_last=frame_last
else 
   frame_last=nframes
end if

!
!    Allocate global coordinate arrays
!
allocate(xyz(3,natoms,nframes),stat=alloc_stat,errmsg=alloc_err)
if (alloc_stat .ne. 0) then
   write(*,*) "Allocation of xyz array failed:",trim(alloc_err)
   write(*,*) "Please read in a smaller trajectory or increase memory!"
   stop
end if
if (npt) then
   allocate(a_vecs(3,nframes))
   allocate(b_vecs(3,nframes))
   allocate(c_vecs(3,nframes))
end if
allocate(keep_atoms(natoms,nframes))
!
!    Read in the coordinates of the frames and correct for image flags
!

close(14)
open(unit=14,file="XDATCAR",status="old")

eval_stat = .false.
write(*,*)
write(*,*) "Read in the trajectory from XDATCAR..."
line_num=1
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
!
!     A NPT trajectory
!
   if (npt) then
      read(14,*)
      read(14,*) a_vecs(:,i)
      read(14,*) b_vecs(:,i)
      read(14,*) c_vecs(:,i) 
      read(14,*)
      read(14,*)
      read(14,*)
!
!     A NVT trajectory
!
   else
      if (i .eq. 1) then 
         read(14,*)
         read(14,*) 
         read(14,*) 
         read(14,*) 
         read(14,*)
         read(14,*)
         read(14,*)
      end if
   end if
   
!
!     A NVT trajectory
!
   do j=1,natoms
      read(14,'(a)') a120
      read(a120,*,iostat=readstat) act_num(:)
!
!   In long trajectories, problems might arise with double-digit negative numbers!
!   (no space between them)   Insert a space again!
!
      if (readstat .ne. 0) then
         a220 = ""
         endl = 1
         do k=1,119
            atest=a120(k:k+1)
            if (index(atest,"-") .ne. 0) then
               write(a220(1:endl+4),'(a,a)') a220(1:endl),"   -"
               endl = endl + 4
            else
               write(a220(1:endl+1),'(a,a)',iostat=readstat) a220(1:endl),atest
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
            else if (act_num(k) < 0.d0) then
               act_num(k) = act_num(k) + 1d0
            else 
               exit
            end if
         end do
         xyz(k,j,i)=act_num(k)
      end do
   end do
!
!     If read_freq > 1, skip the next read_freq-1 frames
!
   if (read_freq .gt. 1) then
      do j=1,read_freq-1
         if (npt) then
            do k=1,natoms+8
               read(14,*)
            end do
         else
            do k=1,natoms+1
               read(14,*)
            end do
         end if
      end do
   end if
end do
write(*,*) " completed!"
close(14)

write(*,*)
write(*,*) "---------- SETTINGS ---------------------------"
write(*,*) "Number of atoms in the system:",natoms
write(*,*) "Number of frames in the trajectory:",nframes*read_freq
if (read_freq .gt. 1) then
   write(*,'(a,i5,a)') " Every ",read_freq,"th frame will be evaluated."
end if        
write(*,*) "-----------------------------------------------"
write(*,*)
!
!     C: Shift the unit cell along a vector
!
if (shift_cell) then
   eval_stat = .false.
   write(*,'(a,f9.4,a,f9.4,a,f9.4)') " Shift all frames of the XDATCAR along the vector ", &
             & shift_vec(1),",",shift_vec(2),",",shift_vec(3)
   do i=1,nframes
      do j=1,10
         if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
      do j=1,natoms
         do k=1,3
            act_val=xyz(k,j,i)+shift_vec(k)
            do 
               if (act_val > 1.0d0) then
                  act_val=act_val-1.0d0
               else if (act_val < 0.0d0) then
                  act_val=act_val+1.d0
               else 
                  exit
               end if
            end do
            xyz(k,j,i)=act_val
         end do
      end do
   end do
   write(*,*) "completed!"
   write(*,*)
end if
!
!    D: Multiply the unit cell according to given integer tuple
!
if (multiply_cell) then
   eval_stat = .false.
   write(*,'(a,i4,a,i4,a,i4,a)') " Multipy the unit cell of all frames: ",multiply_vec(1), &
               & " times a, ",multiply_vec(2)," times b, ",multiply_vec(3)," times c"
   multiply_prod=multiply_vec(1)*multiply_vec(2)*multiply_vec(3)
   if (npt) then
      do i=1,nframes
         a_vecs(:,i)=a_vecs(:,i)*multiply_vec(1)
         b_vecs(:,i)=b_vecs(:,i)*multiply_vec(2)
         c_vecs(:,i)=c_vecs(:,i)*multiply_vec(3)
      end do
   else
      a_vec=a_vec*multiply_vec(1)
      b_vec=b_vec*multiply_vec(2)
      c_vec=c_vec*multiply_vec(3)
   end if
   el_nums=el_nums*multiply_prod
   allocate(xyz2(3,natoms*multiply_prod,nframes),at_names2(natoms*multiply_prod))
   do i=1,nframes
      do j=1,10
         if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
      pos_new=0
      do j=1,natoms
         do k=1,multiply_vec(1)
            do l=1,multiply_vec(2)
               do m=1,multiply_vec(3)
                  pos_new=pos_new+1
                  xyz2(1,pos_new,i)=(xyz(1,j,i)+1.0*real(k))/real(multiply_vec(1))
                  xyz2(2,pos_new,i)=(xyz(2,j,i)+1.0*real(l))/real(multiply_vec(2))
                  xyz2(3,pos_new,i)=(xyz(3,j,i)+1.0*real(m))/real(multiply_vec(3))
                  at_names2(pos_new)=at_names(j)
               end do
            end do
         end do
      end do
   end do
   natoms=natoms*multiply_prod

   write(*,*) "completed!"
   write(*,*)
else 
!
!    If no multiplication shall be done, simply copy contents of xyz to xyz2
!
   allocate(xyz2(3,natoms,nframes))
   allocate(at_names2(natoms))
   xyz2=xyz
   at_names2=at_names
end if
!
!    E: Pick a frame of the trajectory and print it as POSCAR or xyz file
!
if (pick_frame) then
   if (print_xyz) then
      write(*,'(a,i10,a)') " Frame ",pick_ind," will be written to file pos_pick.xyz"
      open(unit=33,file="pos_pick.xyz",status="replace")
      write(33,*) natoms
      write(33,'(a,i10,a)') " Frame ",pick_ind," picked from XDATCAR by modify_xdatcar"
      do j=1,natoms
         if (npt) then
            xyz_print(1)=(xyz2(1,j,pick_ind)*a_vecs(1,pick_ind)+xyz2(2,j,pick_ind)* &
                         & b_vecs(1,pick_ind)+xyz2(3,j,pick_ind)*c_vecs(1,pick_ind))*factor
            xyz_print(2)=(xyz2(1,j,pick_ind)*a_vecs(2,pick_ind)+xyz2(2,j,pick_ind)* &
                         & b_vecs(2,pick_ind)+xyz2(3,j,pick_ind)*c_vecs(2,pick_ind))*factor
            xyz_print(3)=(xyz2(1,j,pick_ind)*a_vecs(3,pick_ind)+xyz2(2,j,pick_ind)* &
                         & b_vecs(3,pick_ind)+xyz2(3,j,pick_ind)*c_vecs(3,pick_ind))*factor
         else
            xyz_print(1)=(xyz2(1,j,pick_ind)*a_vec(1)+xyz2(2,j,pick_ind)*b_vec(1)+ &
                         & xyz2(3,j,pick_ind)*c_vec(1))*factor
            xyz_print(2)=(xyz2(1,j,pick_ind)*a_vec(2)+xyz2(2,j,pick_ind)*b_vec(2)+ &
                         & xyz2(3,j,pick_ind)*c_vec(2))*factor
            xyz_print(3)=(xyz2(1,j,pick_ind)*a_vec(3)+xyz2(2,j,pick_ind)*b_vec(3)+ &
                      & xyz2(3,j,pick_ind)*c_vec(3))*factor
         end if
         write(34,*) at_names2(j),xyz_print(:)
      end do
      close(33)  
      write(*,*) "completed!"
      write(*,*)
   else 
      write(*,'(a,i10,a)') " Frame ",pick_ind," will be written to file POSCAR_pick"
      open(unit=33,file="POSCAR_pick",status="replace")
      write(33,'(a,i10,a)') " Frame ",pick_ind," picked from XDATCAR by modify_xdatcar"
      write(33,*) factor
      if (npt) then
         write(33,'(3f15.6)') a_vecs(:,pick_ind)
         write(33,'(3f15.6)') b_vecs(:,pick_ind)
         write(33,'(3f15.6)') c_vecs(:,pick_ind)
      else
         write(33,'(3f15.6)') a_vec
         write(33,'(3f15.6)') b_vec
         write(33,'(3f15.6)') c_vec
      end if
      do j=1,nelems
         write(33,'(a,a)',advance="no") el_names(j),"  "
      end do
      write(33,*)
      do j=1,nelems
         write(33,'(i6,a)',advance="no") el_nums(j)," "
      end do
      write(33,*)

      write(33,*) "Direct configuration=  ",pick_ind
      do j=1,natoms
         write(33,'(3f15.8)') xyz2(:,j,pick_ind)
      end do
      close(33)
      write(*,*) "completed!"
      write(*,*)
   end if
end if

!
!    E: Remove atoms from all frames of the trajectory that are within a defined range
!
keep_atoms=.true.
if (remove_mode) then
   if (.not. npt) then
      npt=.true.
      allocate(a_vecs(3,nframes))
      allocate(b_vecs(3,nframes))
      allocate(c_vecs(3,nframes))
      do i=1,nframes
         a_vecs(:,i)=a_vec
         b_vecs(:,i)=b_vec
         c_vecs(:,i)=c_vec
      end do
   end if

   do i=1,nframes
      do j=1,natoms
         if (remove_dim .eq. "x") then
            if (remove_sign .eq. "gt") then
               if (xyz2(1,j,i)*a_vecs(1,i) .gt. remove_border) then
                  keep_atoms(j,i) = .false.     
               end if   
            else if (remove_sign .eq. "lt") then
               if (xyz2(1,j,i)*a_vecs(1,i) .lt. remove_border) then
                  keep_atoms(j,i) = .false.
               end if   
            end if
         else if (remove_dim .eq. "y") then
            if (remove_sign .eq. "gt") then
               if (xyz2(2,j,i)*b_vecs(2,i) .gt. remove_border) then
                  keep_atoms(j,i) = .false.
               end if
            else if (remove_sign .eq. "lt") then
               if (xyz2(2,j,i)*b_vecs(2,i) .lt. remove_border) then
                  keep_atoms(j,i) = .false.
               end if
            end if
         else if (remove_dim .eq. "z") then
            if (remove_sign .eq. "gt") then
               if (xyz2(3,j,i)*c_vecs(3,i) .gt. remove_border) then
                  keep_atoms(j,i) = .false.
               end if
            else if (remove_sign .eq. "lt") then
               if (xyz2(3,j,i)*c_vecs(3,i) .lt. remove_border) then
                  keep_atoms(j,i) = .false.
               end if
            end if
         end if   
      end do
   end do
end if
!
!    D: Remove all net translation and rotation from the system (should be an
!      isolated molecule). Each frame is placed with its center of mass in the 
!      middle of the unit cell; and the content is aligned with respect to its 
!      principial moments of inertia
!
if (remove_transrot) then
!
!    First, remove all image flags to avoid sudden jumps of COM
!
   do i=frame_first,frame_last
      if (i .gt. frame_first) then
         do j=1,natoms     
            do k=1,3
               if (xyz2(k,j,i) .gt. xyz2(k,j,i-1)+0.5) then
                  xyz2(k,j,i) = xyz2(k,j,i)-1.0
               end if
               if (xyz2(k,j,i) .lt. xyz2(k,j,i-1)-0.5) then
                  xyz2(k,j,i) = xyz2(k,j,i)+1.0
               end if
            end do
         end do
      end if
   end do
!
!    Determine the atomic masses 
!
   allocate(at_masses(natoms))
   allocate(vel_act(3,natoms),vel_old(3,natoms))
   do i=1,natoms
      call atommass(at_names(i),at_masses(i)) 
   end do
   totmass=sum(at_masses)
!
!    Determine the center of mass of the current trajectory frame
!
   do i=frame_first,frame_last
      com_act=0.d0   
      do j=1,natoms
         com_act=com_act+xyz2(:,j,i)*at_masses(j)
      end do
      com_act=com_act/totmass    
      write(*,*) "com",com_act
!      do j=1,natoms
!         do k=1,3
!            xyz2(k,j,i)=xyz2(k,j,i)-com_act(k)+0.5d0
!         end do
!      end do 

!
!    Rotate the current frame to the standard moments of inertia
!
      vel_act=0
      if (i .gt. frame_first) then
         do j=1,natoms
            vel_act(:,j)=xyz2(:,j,i)-xyz2(:,j,i-1)
         end do 
         vel_old=vel_act
!
!     Compute linear velocity of the system center of mass
!
         do j = 1, natoms
            weigh = at_masses(j)
            do k = 1, 3
               vtot(j) = vtot(j) + vel_act(k,j)*weigh
            end do
         end do
         vtot=vtot/totmass
         write(*,*) "vtot",vtot
!
!     Compute angular momentum for overall system
!
         mang = 0.0d0
         xtot=com_act(1)
         ytot=com_act(2)
         ztot=com_act(3)

         do j = 1, natoms
            weigh = at_masses(j)
            mang(1) = mang(1) + (xyz2(2,j,i)*vel_act(3,j)-xyz2(3,j,i)*vel_act(2,j))*weigh
            mang(2) = mang(2) + (xyz2(3,j,i)*vel_act(1,j)-xyz2(1,j,i)*vel_act(3,j))*weigh
            mang(3) = mang(3) + (xyz2(1,j,i)*vel_act(2,j)-xyz2(2,j,i)*vel_act(1,j))*weigh
         end do

         mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
         mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
         mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass
         write(*,*) "mang",mang
!
!     calculate the moment of inertia tensor
!
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0

         do j = 1, natoms
            weigh = at_masses(j)
            xdel = xyz2(1,j,i) - xtot
            ydel = xyz2(2,j,i) - ytot
            zdel = xyz2(3,j,i) - ztot
            xx = xx + xdel*xdel*weigh
            xy = xy + xdel*ydel*weigh
            xz = xz + xdel*zdel*weigh
            yy = yy + ydel*ydel*weigh
            yz = yz + ydel*zdel*weigh
            zz = zz + zdel*zdel*weigh
         end do
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
!
!     avoid bad behavior (singularity) for diatomic systems
!

         if (natoms .le. 2) then
            eps = 0.000001d0
            tensor(1,1) = tensor(1,1) + eps
            tensor(2,2) = tensor(2,2) + eps
            tensor(3,3) = tensor(3,3) + eps
         end if
! 
!     invert the moment of inertia tensor
!
         call invert (3,tensor)
!
!     compute angular velocity 
!
         do k = 1, 3
            vang(k) = 0.0d0
            do l = 1, 3
               vang(k) = vang(k) + tensor(k,l)*mang(l)
            end do
        end do


!
!     eliminate any translation of the overall system
!
         do j = 1, natoms
            do k = 1, 3
               vel_act(k,j) = vel_act(k,j) - vtot(j)
            end do
         end do
!
!     eliminate any rotation of the overall system
!
         do j = 1, natoms
            xdel = xyz2(1,j,i) - xtot
            ydel = xyz2(2,j,i) - ytot
            zdel = xyz2(3,j,i) - ztot
            vel_act(1,j) = vel_act(1,j) - vang(2)*zdel + vang(3)*ydel
            vel_act(2,j) = vel_act(2,j) - vang(3)*xdel + vang(1)*zdel
            vel_act(3,j) = vel_act(3,j) - vang(1)*ydel + vang(2)*xdel
         end do
!
!     Now correct the xyz structure from the velocity: 
!     Subtract the old velocity (difference to previous frame) 
!     and then add the new velocity
! 
         do j=1,natoms
            xyz2(:,j,i)=xyz2(:,j,i)-vel_old(:,j)
            xyz2(:,j,i)=xyz2(:,j,i)+vel_act(:,j)
         end do 
      end if
   end do
    

end if        


!
!    D: Write structure to XYZ trajectory file
!      During this, convert coordinates to cartesian!
!      If not ordered, write a modified XDATCAR file (in direct coordinates) instead
!      Print nothing if the pick_frame option is chosen
!
eval_stat = .false.
if (print_xyz) then
   write(*,*) "Write trajectory in xyz format to file xdat_mod.xyz"
   write(*,'(a,i10,a,i10,a)') " Frame ",frame_first," to frame ",frame_last," will be written."
   open(unit=34,file="xdat_mod.xyz",status="replace")
   do i=frame_first,frame_last
      do j=1,10
         if (real(i)/real(frame_last-frame_first) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
      counter=0
      do j=1,natoms
         if (keep_atoms(j,i)) then
            counter=counter+1
         end if
      end do
      write(34,*) counter
      write(34,*) "Trajectory converted from XDATCAR file via modify_xdatcar"
      do j=1,natoms
!
!     If the smooth mode is activated, remove all sudden jumps between frames
!      due to changed image flags
!
         if (smooth_mode) then
            if (i .gt. frame_first) then
               do k=1,3
                  if (xyz2(k,j,i) .gt. xyz2(k,j,i-1)+0.5) then
                     xyz2(k,j,i) = xyz2(k,j,i)-1.0
                  end if
                  if (xyz2(k,j,i) .lt. xyz2(k,j,i-1)-0.5) then
                     xyz2(k,j,i) = xyz2(k,j,i)+1.0
                  end if
               end do
            end if
         end if
!
!     If removal mode is activated: only print the remaining atoms
!
         if (keep_atoms(j,i)) then
            if (npt) then
               xyz_print(1)=(xyz2(1,j,i)*a_vecs(1,i)+xyz2(2,j,i)*b_vecs(1,i)+ &
                            & xyz2(3,j,i)*c_vecs(1,i))*factor
               xyz_print(2)=(xyz2(1,j,i)*a_vecs(2,i)+xyz2(2,j,i)*b_vecs(2,i)+ &
                            & xyz2(3,j,i)*c_vecs(2,i))*factor
               xyz_print(3)=(xyz2(1,j,i)*a_vecs(3,i)+xyz2(2,j,i)*b_vecs(3,i)+ & 
                            & xyz2(3,j,i)*c_vecs(3,i))*factor
            else 
               xyz_print(1)=(xyz2(1,j,i)*a_vec(1)+xyz2(2,j,i)*b_vec(1)+ &
                            & xyz2(3,j,i)*c_vec(1))*factor
               xyz_print(2)=(xyz2(1,j,i)*a_vec(2)+xyz2(2,j,i)*b_vec(2)+ &
                            & xyz2(3,j,i)*c_vec(2))*factor
               xyz_print(3)=(xyz2(1,j,i)*a_vec(3)+xyz2(2,j,i)*b_vec(3)+ &
                            & xyz2(3,j,i)*c_vec(3))*factor
            end if
            write(34,*) at_names2(j),xyz_print(:)
         end if
      end do
   end do
   close(34)
   write(*,*) "completed!"
else 
   if (print_npt) then
      if (.not. npt) then
         npt=.true.
         allocate(a_vecs(3,nframes))
         allocate(b_vecs(3,nframes))
         allocate(c_vecs(3,nframes))
         do i=1,nframes
            a_vecs(:,i)=a_vec
            b_vecs(:,i)=b_vec
            c_vecs(:,i)=c_vec
         end do
      end if
      write(*,*) "The -print_npt command is given, the trajectory is written as NpT"
   end if
   if (shift_cell .or. multiply_cell .or. print_last .or. (read_freq .gt. 1) .or. &
                   & print_npt .or. remove_transrot) then
      write(*,*) "Write trajectory in VASP format to file XDATCAR_mod"
      write(*,'(a,i10,a,i10,a)') " Frame ",frame_first," to frame ",frame_last," will be written."
      open(unit=34,file="XDATCAR_mod",status="replace")
      do i=frame_first,frame_last
         do j=1,10
            if (real(i)/real(frame_last-frame_first) .gt. real(j)*0.1d0) then
               if (.not. eval_stat(j)) then
                  write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
                  eval_stat(j) = .true.
               end if
            end if
         end do
         if (npt) then
!
!    For removed atoms: decrease the number of total atoms per element each frame
!
            el_nums_tmp=0
            counter = 1
            do j=1,nelems
               do k=1,el_nums(j)
                  if (keep_atoms(counter,j)) then
                     el_nums_tmp(j) = el_nums_tmp(j)+1
                  end if
                  counter = counter +1                  
               end do
            end do

            write(34,*) "NpT Trajectory written by modify_xdatcar"
            write(34,*) factor
            write(34,'(3f15.6)') a_vecs(:,i)
            write(34,'(3f15.6)') b_vecs(:,i)
            write(34,'(3f15.6)') c_vecs(:,i)
            do j=1,nelems
               if (el_nums_tmp(j) .gt. 0) then
                  write(34,'(a,a)',advance="no") el_names(j),"  "
               end if
            end do
            write(34,*)
            do j=1,nelems
               if (el_nums_tmp(j) .gt. 0) then
                  write(34,'(i6,a)',advance="no") el_nums_tmp(j)," "
               end if
            end do
            write(34,*)
         else if (i .eq. frame_first) then
            write(34,*) "NVT Trajectory written by modify_xdatcar"
            write(34,*) factor
            write(34,'(3f15.6)') a_vec
            write(34,'(3f15.6)') b_vec
            write(34,'(3f15.6)') c_vec
            do j=1,nelems
               write(34,'(a,a)',advance="no") el_names(j),"  "
            end do
            write(34,*)
            do j=1,nelems
               write(34,'(i6,a)',advance="no") el_nums(j)," "
            end do
            write(34,*)
         end if
 
         write(34,*) "Direct configuration=  ",i
         do j=1,natoms
!
!     If the smooth mode is activated, remove all sudden jumps between frames
!      due to changed image flags
!
            if (smooth_mode) then
               if (i .gt. frame_first) then
                  do k=1,3
                     if (xyz2(k,j,i) .gt. xyz2(k,j,i-1)+0.5) then
                        xyz2(k,j,i) = xyz2(k,j,i)-1.0
                     end if
                     if (xyz2(k,j,i) .lt. xyz2(k,j,i-1)-0.5) then
                        xyz2(k,j,i) = xyz2(k,j,i)+1.0
                     end if
                  end do
               end if
            end if

            if (keep_atoms(j,i)) then
               write(34,'(3f15.8)') xyz2(:,j,i) 
            end if
         end do
      end do
      close(34)
      write(*,*) "completed!"
   end if
end if

write(*,*)
write(*,*) "modify_xdatcar has finished all tasks, goodbye!"
write(*,*)

end program modify_xdatcar  

!
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine atommass  --  initialize atomic masses        ##
!     ##                                                           ##
!     ###############################################################
!
!     Normally the atomic masses needed for dynamics are read in
!     from the ffield file.
!     But in the case of EVB-QMDFF no masses are specified in the ffield
!     so they need to defined here in this file.
!     The masses from H up to Rn are availiable
!     They are taken from "Das große Tafelwerk, 1. Auflage 2003"
!
!     part of EVB --> Taken from Caracal (J. Steffen et al.)
!
subroutine atommass(adum,mass_act)
implicit none
character(len=2)::adum
real(kind=8)::mass_act
!
!     initialize the element-masses
!     Convert all letters to upcase letters in order to avoid
!     errors with atom symbols!
!
call upcase(adum)
if (adum .eq. "H") then
   mass_act=1.00782503207d0
else if (adum .eq. "D") then
   mass_act=2.0141017778d0
else if (adum .eq. "HE") then
   mass_act=4.00260d0
else if (adum .eq. "LI") then
   mass_act=6.94000d0
else if (adum .eq. "BE") then
   mass_act=9.01218d0
else if (adum .eq. "B") then
   mass_act=10.81000d0
else if (adum .eq. "C") then
   mass_act=12.00000d0
else if (adum .eq. "N") then
   mass_act=14.0030740048d0
else if (adum .eq. "O") then
   mass_act=15.99491461956d0
else if (adum .eq. "F") then
   mass_act=18.99840d0
else if (adum .eq. "NE") then
   mass_act=20.18d0
else if (adum .eq. "NA") then
   mass_act=23.0d0
else if (adum .eq. "MG") then
   mass_act=24.31d0
else if (adum .eq. "AL") then
   mass_act=26.98d0
else if (adum .eq. "SI") then
   mass_act=28.09d0
else if (adum .eq. "P") then
   mass_act=30.97d0
else if (adum .eq. "S") then
   mass_act=32.06000d0
else if (adum .eq. "CL") then
   mass_act=34.96885268d0
else if (adum .eq. "AR") then
   mass_act=39.95d0
else if (adum .eq. "K") then
   mass_act=39.1d0
else if (adum .eq. "CA") then
   mass_act=40.08d0
else if (adum .eq. "SC") then
   mass_act=44.96d0
else if (adum .eq. "TI") then
   mass_act=47.88d0
else if (adum .eq. "V") then
   mass_act=50.94d0
else if (adum .eq. "CR") then
   mass_act=52.0d0
else if (adum .eq. "MN") then
   mass_act=54.94d0
else if (adum .eq. "FE") then
   mass_act=55.85d0
else if (adum .eq. "CO") then
   mass_act=58.93d0
else if (adum .eq. "NI") then
   mass_act=58.69d0
else if (adum .eq. "CU") then
   mass_act=63.55d0
else if (adum .eq. "ZN") then
   mass_act=65.39d0
else if (adum .eq. "GA") then
   mass_act=69.72d0
else if (adum .eq. "GE") then
   mass_act=72.61d0
else if (adum .eq. "AS") then
   mass_act=74.92d0
else if (adum .eq. "SE") then
   mass_act=78.96d0
else if (adum .eq. "BR") then
   mass_act=78.9183371d0
else if (adum .eq. "KR") then
   mass_act=83.80d0
else if (adum .eq. "RB") then
   mass_act=85.47d0
else if (adum .eq. "SR") then
   mass_act=87.62d0
else if (adum .eq. "Y") then
   mass_act=88.91d0
else if (adum .eq. "ZR") then
   mass_act=91.22d0
else if (adum .eq. "NB") then
   mass_act=92.91d0
else if (adum .eq. "MO") then
   mass_act=95.94d0
else if (adum .eq. "TC") then
   mass_act=98d0
else if (adum .eq. "RU") then
   mass_act=101.07d0
else if (adum .eq. "RH") then
   mass_act=102.91d0
else if (adum .eq. "PD") then
   mass_act=106.42d0
else if (adum .eq. "AG") then
   mass_act=107.87d0
else if (adum .eq. "CD") then
   mass_act=112.41d0
else if (adum .eq. "IN") then
   mass_act=114.82d0
else if (adum .eq. "SN") then
   mass_act=118.71d0
else if (adum .eq. "SB") then
   mass_act=121.76d0
else if (adum .eq. "TE") then
   mass_act=127.60d0
else if (adum .eq. "I") then
   mass_act=125.90d0
else if (adum .eq. "XE") then
   mass_act=131.29d0
else if (adum .eq. "CS") then
   mass_act=132.91d0
else if (adum .eq. "BA") then
   mass_act=137.33d0
else if (adum .eq. "LA") then
   mass_act=138.91d0
else if (adum .eq. "HF") then
   mass_act=178.49d0
else if (adum .eq. "TA") then
   mass_act=180.95d0
else if (adum .eq. "W") then
   mass_act=183.95d0
else if (adum .eq. "RE") then
   mass_act=186.21d0
else if (adum .eq. "OS") then
   mass_act=190.23d0
else if (adum .eq. "IR") then
   mass_act=192.22d0
else if (adum .eq. "PT") then
   mass_act=195.08d0
else if (adum .eq. "AU") then
   mass_act=196.97d0
else if (adum .eq. "HG") then
   mass_act=200.59d0
else if (adum .eq. "TL") then
   mass_act=204.38d0
else if (adum .eq. "PB") then
   mass_act=207.2d0
else if (adum .eq. "BI") then
   mass_act=208.98d0
else if (adum .eq. "PO") then
   mass_act=209d0
else if (adum .eq. "AT") then
   mass_act=210d0
else if (adum .eq. "RN") then
   mass_act=222d0
end if

return
end subroutine atommass

!
!     subroutine upcase: converts a text string to all upper
!            case letters
!
!     part of EVB  --> Taken from Caracal (J. Steffen et al.)
!
subroutine upcase (string)
implicit none
integer::i,size,len
integer::code,ichar
character(len=1)::char
character(len=1)::letter
character(len=*)::string
!
!     convert lower case to upper case one letter at a time
!
size = len(string)
do i = 1, size
   letter = string(i:i)
   code = ichar(letter)
   if (letter.ge.'a' .and. letter.le.'z') &
      &   string(i:i) = char(code-32)
end do

return
end subroutine upcase

!
!     subroutine invert: inverts a matrix using the Gauss-Jordan method
!     needed for dynamics
!     variables and parameters:
!     n     dimension of the matrix to be inverted
!     a     matrix to invert; contains inverse on exit
!
!     part of EVB --> Taken from Caracal (J. Steffen et al.)
!
subroutine invert (n,a)
implicit none
integer::i,j,k,n
integer::icol,irow
integer,allocatable::ipivot(:)
integer,allocatable::indxc(:)
integer,allocatable::indxr(:)
real(kind=8)::big,temp
real(kind=8)::pivot
real(kind=8)::a(n,*)
!
!     perform dynamic allocation of some local arrays
!
allocate (ipivot(n))
allocate (indxc(n))
allocate (indxr(n))
!
!     perform matrix inversion via the Gauss-Jordan algorithm
!
do i = 1, n
   ipivot(i) = 0
end do
do i = 1, n
   big = 0.0d0
   do j = 1, n
      if (ipivot(j) .ne. 1) then
         do k = 1, n
            if (ipivot(k) .eq. 0) then
               if (abs(a(j,k)) .ge. big) then
                  big = abs(a(j,k))
                  irow = j
                  icol = k
               end if
            else if (ipivot(k) .gt. 1) then
               write (*,'(/," INVERT  --  Cannot Invert a Singular Matrix")')
               stop
            end if
         end do
      end if
   end do
   ipivot(icol) = ipivot(icol) + 1
   if (irow .ne. icol) then
      do j = 1, n
         temp = a(irow,j)
         a(irow,j) = a(icol,j)
         a(icol,j) = temp
      end do
   end if
   indxr(i) = irow
   indxc(i) = icol
   if (a(icol,icol) .eq. 0.0d0) then
      write (*,'(/," INVERT  --  Cannot Invert a Singular Matrix")')
      stop
   end if
   pivot = a(icol,icol)
   a(icol,icol) = 1.0d0
   do j = 1, n
      a(icol,j) = a(icol,j) / pivot
   end do
   do j = 1, n
      if (j .ne. icol) then
         temp = a(j,icol)
         a(j,icol) = 0.0d0
         do k = 1, n
            a(j,k) = a(j,k) - a(icol,k)*temp
         end do
      end if
   end do
end do
do i = n, 1, -1
   if (indxr(i) .ne. indxc(i)) then
      do k = 1, n
         temp = a(k,indxr(i))
         a(k,indxr(i)) = a(k,indxc(i))
         a(k,indxc(i)) = temp
      end do
   end if
end do
!
!     perform deallocation of some local arrays
!
deallocate (ipivot)
deallocate (indxc)
deallocate (indxr)

return
end subroutine invert
